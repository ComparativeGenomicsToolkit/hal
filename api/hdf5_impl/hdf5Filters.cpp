/*
 * Copyright (C) 2024 by Glenn Hickey (hickey@soe.ucsc.edu)
 * Copyright (C) 2012-2024 by UCSC Computational Genomics Lab
 *
 * Released under the MIT license, see LICENSE.txt
 *
 * LZ4 and Zstd HDF5 compression filters for HAL.
 *
 * Filter IDs follow the HDF5 registered filter convention:
 *   LZ4:  32004 (compatible with nexusformat/HDF5-External-Filter-Plugins)
 *   Zstd: 32015 (compatible with aparamon/HDF5Plugin-Zstandard)
 *
 * The cd_values parameter conventions match the standard implementations
 * so files created by HAL can be read by any HDF5 tool with the
 * corresponding plugin installed, and vice versa.
 */

#include "hdf5Filters.h"
#include <cstdlib>
#include <cstring>
#include <lz4.h>
#include <lz4hc.h>
#include <zstd.h>

/* ========================================================================
 * LZ4 filter (ID 32004)
 *
 * cd_values[0] = blockSize (0 means use LZ4 default: 1<<17 = 128KB)
 *
 * Wire format per chunk (matches nexusformat convention):
 *   [4 bytes big-endian] original (uncompressed) size
 *   [rest]               LZ4 compressed data
 * ======================================================================== */

static size_t H5Z_filter_lz4(unsigned int flags, size_t cd_nelmts,
                              const unsigned int cd_values[], size_t nbytes,
                              size_t *buf_size, void **buf) {
    if (flags & H5Z_FLAG_REVERSE) {
        /* Decompress */
        const char *src = (const char *)*buf;
        if (nbytes < 4) return 0;

        /* Read original size (big-endian 4 bytes) */
        uint32_t origSize = ((uint32_t)(unsigned char)src[0] << 24) |
                            ((uint32_t)(unsigned char)src[1] << 16) |
                            ((uint32_t)(unsigned char)src[2] << 8)  |
                            ((uint32_t)(unsigned char)src[3]);

        void *outbuf = malloc(origSize);
        if (!outbuf) return 0;

        int decompressed = LZ4_decompress_safe(src + 4, (char *)outbuf,
                                               (int)(nbytes - 4), (int)origSize);
        if (decompressed < 0) {
            free(outbuf);
            return 0;
        }

        free(*buf);
        *buf = outbuf;
        *buf_size = origSize;
        return (size_t)origSize;
    } else {
        /* Compress */
        const char *src = (const char *)*buf;
        int srcSize = (int)nbytes;
        int maxDst = LZ4_compressBound(srcSize);
        void *outbuf = malloc(4 + maxDst);
        if (!outbuf) return 0;

        /* Write original size as big-endian 4 bytes */
        char *dst = (char *)outbuf;
        dst[0] = (char)((nbytes >> 24) & 0xFF);
        dst[1] = (char)((nbytes >> 16) & 0xFF);
        dst[2] = (char)((nbytes >> 8)  & 0xFF);
        dst[3] = (char)((nbytes)       & 0xFF);

        int level = (cd_nelmts > 1 && cd_values[1] > 0) ? (int)cd_values[1] : 0;
        int compressed;
        if (level > 0) {
            compressed = LZ4_compress_HC(src, dst + 4, srcSize, maxDst, level);
        } else {
            compressed = LZ4_compress_default(src, dst + 4, srcSize, maxDst);
        }
        if (compressed <= 0) {
            free(outbuf);
            return 0;
        }

        free(*buf);
        *buf = outbuf;
        *buf_size = 4 + compressed;
        return 4 + (size_t)compressed;
    }
}

static const H5Z_class2_t H5Z_LZ4_class = {
    H5Z_CLASS_T_VERS,      /* H5Z_class_t version */
    H5Z_FILTER_LZ4,        /* Filter id number */
    1,                      /* encoder_present flag */
    1,                      /* decoder_present flag */
    "lz4",                  /* Filter name for debugging */
    NULL,                   /* The "can apply" callback */
    NULL,                   /* The "set local" callback */
    H5Z_filter_lz4         /* The actual filter function */
};

/* ========================================================================
 * Zstd filter (ID 32015)
 *
 * cd_values[0] = compression level (default ZSTD_CLEVEL_DEFAULT = 3)
 *
 * Wire format: raw zstd frame (self-describing, includes original size).
 * This matches the aparamon/HDF5Plugin-Zstandard convention.
 * ======================================================================== */

static size_t H5Z_filter_zstd(unsigned int flags, size_t cd_nelmts,
                               const unsigned int cd_values[], size_t nbytes,
                               size_t *buf_size, void **buf) {
    if (flags & H5Z_FLAG_REVERSE) {
        /* Decompress */
        const void *src = *buf;
        unsigned long long origSize = ZSTD_getFrameContentSize(src, nbytes);
        if (origSize == ZSTD_CONTENTSIZE_UNKNOWN || origSize == ZSTD_CONTENTSIZE_ERROR) {
            /* Fall back: allocate a generous buffer */
            origSize = nbytes * 8;
        }

        void *outbuf = malloc((size_t)origSize);
        if (!outbuf) return 0;

        size_t decompressed = ZSTD_decompress(outbuf, (size_t)origSize, src, nbytes);
        if (ZSTD_isError(decompressed)) {
            free(outbuf);
            return 0;
        }

        free(*buf);
        *buf = outbuf;
        *buf_size = decompressed;
        return decompressed;
    } else {
        /* Compress */
        int level = (cd_nelmts > 0) ? (int)cd_values[0] : ZSTD_CLEVEL_DEFAULT;
        size_t maxDst = ZSTD_compressBound(nbytes);
        void *outbuf = malloc(maxDst);
        if (!outbuf) return 0;

        size_t compressed = ZSTD_compress(outbuf, maxDst, *buf, nbytes, level);
        if (ZSTD_isError(compressed)) {
            free(outbuf);
            return 0;
        }

        free(*buf);
        *buf = outbuf;
        *buf_size = compressed;
        return compressed;
    }
}

static const H5Z_class2_t H5Z_ZSTD_class = {
    H5Z_CLASS_T_VERS,      /* H5Z_class_t version */
    H5Z_FILTER_ZSTD,       /* Filter id number */
    1,                      /* encoder_present flag */
    1,                      /* decoder_present flag */
    "zstd",                 /* Filter name for debugging */
    NULL,                   /* The "can apply" callback */
    NULL,                   /* The "set local" callback */
    H5Z_filter_zstd         /* The actual filter function */
};

/* ========================================================================
 * Registration
 * ======================================================================== */

void halRegisterCompressors(void) {
    static bool registered = false;
    if (registered) return;

    H5Zregister(&H5Z_LZ4_class);
    H5Zregister(&H5Z_ZSTD_class);
    registered = true;
}
