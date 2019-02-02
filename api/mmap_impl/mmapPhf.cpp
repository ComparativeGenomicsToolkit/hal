/* ==========================================================================
 * phf.cc - Tiny perfect hash function library.
 * --------------------------------------------------------------------------
 * Copyright (c) 2014-2015  William Ahern
 *
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to permit
 * persons to whom the Software is furnished to do so, subject to the
 * following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
 * NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
 * USE OR OTHER DEALINGS IN THE SOFTWARE.
 * ==========================================================================
 * This code was taken PHF package by William Ahern,
 * https://github.com/wahern/phf for use in HAL.  It has been modified from
 * the original source.  Copied at commit
 * 2b781079946401ce160eb577b4534e94b763242a.
 * ==========================================================================
 */
#include <halDefs.h>
#include <limits.h>   /* CHAR_BIT SIZE_MAX */
#include <inttypes.h> /* PRIu32 PRIu64 PRIx64 */
#include <stdint.h>   /* UINT32_C UINT64_C uint32_t uint64_t */
#include <stdlib.h>   /* abort(3) calloc(3) free(3) qsort(3) */
#include <string.h>   /* memset(3) */
#include <errno.h>    /* errno */
#include <assert.h>   /* assert(3) */
#if !PHF_NO_LIBCXX
#include <string>     /* std::string */
#endif
#include <vector>
#include <algorithm>

#include "mmapPhf.h"

#define PHF_NO_COMPUTED_GOTOS 1

#ifdef __clang__
#pragma clang diagnostic ignored "-Wunused-function"
#if __cplusplus < 201103L
#pragma clang diagnostic ignored "-Wc++11-long-long"
#endif
#elif PHF_GNUC_PREREQ(4, 6)
#pragma GCC diagnostic ignored "-Wunused-function"
#if __cplusplus < 201103L
#pragma GCC diagnostic ignored "-Wlong-long"
#pragma GCC diagnostic ignored "-Wformat" // %zu
#endif
#endif


/*
 * M A C R O  R O U T I N E S
 *
 * Mostly copies of <sys/param.h>
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#define PHF_BITS(T) (sizeof (T) * CHAR_BIT)
#define PHF_HOWMANY(x, y) (((x) + ((y) - 1)) / (y))
#define PHF_MIN(a, b) (((a) < (b))? (a) : (b))
#define PHF_MAX(a, b) (((a) > (b))? (a) : (b))
#define PHF_ROTL(x, y) (((x) << (y)) | ((x) >> (PHF_BITS(x) - (y))))
#define PHF_COUNTOF(a) (sizeof (a) / sizeof *(a))


/*
 * M O D U L A R  A R I T H M E T I C  R O U T I N E S
 *
 * Two modular reduction schemes are supported: bitwise AND and naive
 * modular division. For bitwise AND we must round up the values r and m to
 * a power of 2.
 *
 * TODO: Implement and test Barrett reduction as alternative to naive
 * modular division.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* round up to nearest power of 2 */
static inline size_t phf_powerup(size_t i) {
#if defined SIZE_MAX
	i--;
	i |= i >> 1;
	i |= i >> 2;
	i |= i >> 4;
	i |= i >> 8;
	i |= i >> 16;
#if SIZE_MAX != 0xffffffffu
	i |= i >> 32;
#endif
	return ++i;
#else
#error No SIZE_MAX defined
#endif
} /* phf_powerup() */

static inline uint64_t phf_a_s_mod_n(uint64_t a, uint64_t s, uint64_t n) {
	uint64_t v;

	assert(n <= UINT32_MAX);

	v = 1;
	a %= n;

	while (s > 0) {
		if (s % 2 == 1)
			v = (v * a) % n;
		a = (a * a) % n;
		s /= 2;
	}

	return v;
} /* phf_a_s_mod_n() */

/*
 * Rabin-Miller primality test adapted from Niels Ferguson and Bruce
 * Schneier, "Practical Cryptography" (Wiley, 2003), 201-204.
 */
static inline bool phf_witness(uint64_t n, uint64_t a, uint64_t s, uint64_t t) {
	uint64_t v, i;

	assert(a > 0 && a < n);
	assert(n <= UINT32_MAX);

	if (1 == (v = phf_a_s_mod_n(a, s, n)))
		return 1;

	for (i = 0; v != n - 1; i++) {
		if (i == t - 1)
			return 0;
		v = (v * v) % n;
	}

	return 1;
} /* phf_witness() */

static inline bool phf_rabinmiller(uint64_t n) {
	/*
	 * Witness 2 is deterministic for all n < 2047. Witnesses 2, 7, 61
	 * are deterministic for all n < 4,759,123,141.
	 */
	static const int witness[] = { 2, 7, 61 };
	uint64_t s, t, i;

	assert(n <= UINT32_MAX);

	if (n < 3 || n % 2 == 0)
		return 0;

	/* derive 2^t * s = n - 1 where s is odd */
	s = n - 1;
	t = 0;
	while (s % 2 == 0) {
		s /= 2;
		t++;
	}

	/* NB: witness a must be 1 <= a < n */
	if (n < 2027)
		return phf_witness(n, 2, s, t);

	for (i = 0; i < PHF_COUNTOF(witness); i++) {
		if (!phf_witness(n, witness[i], s, t))
			return 0;
	}

	return 1;
} /* phf_rabinmiller() */

static inline bool phf_isprime(size_t n) {
	static const char map[] = { 0, 1, 2, 3, 0, 5, 0, 7 };
	size_t i;

	if (n < PHF_COUNTOF(map))
		return map[n];

	for (i = 2; i < PHF_COUNTOF(map); i++) {
		if (map[i] && (n % map[i] == 0))
			return 0;
	}

	return phf_rabinmiller(n);
} /* phf_isprime() */

static inline size_t phf_primeup(size_t n) {
	/* NB: 4294967291 is largest 32-bit prime */
	if (n > 4294967291)
		return 0;

	while (n < SIZE_MAX && !phf_isprime(n))
		n++;

	return n;
} /* phf_primeup() */


/*
 * B I T M A P  R O U T I N E S
 *
 * We use a bitmap to track output hash occupancy when searching for
 * displacement values.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

typedef unsigned long phf_bits_t;

static inline bool phf_isset(phf_bits_t *set, size_t i) {
	return set[i / PHF_BITS(*set)] & ((size_t)1 << (i % PHF_BITS(*set)));
} /* phf_isset() */

static inline void phf_setbit(phf_bits_t *set, size_t i) {
	set[i / PHF_BITS(*set)] |= ((size_t)1 << (i % PHF_BITS(*set)));
} /* phf_setbit() */

static inline void phf_clrbit(phf_bits_t *set, size_t i) {
	set[i / PHF_BITS(*set)] &= ~((size_t)1 << (i % PHF_BITS(*set)));
} /* phf_clrbit() */

static inline void phf_clrall(phf_bits_t *set, size_t n) {
	memset(set, '\0', PHF_HOWMANY(n, PHF_BITS(*set)) * sizeof *set);
} /* phf_clrall() */


/*
 * K E Y  D E D U P L I C A T I O N
 *
 * Auxiliary routine to ensure uniqueness of each key in array.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

namespace PHF {
	namespace Uniq {
		static bool operator!=(const phf_string_t &a, const phf_string_t &b) {
			return a.n != b.n || 0 != memcmp(a.p, b.p, a.n);
		}

		template<typename T>
		static int cmp(const T *a, const T *b) {
			if (*a > *b)
				return -1;
			if (*a < *b)
				return 1;
			return 0;
		} /* cmp() */

		template<>
		int cmp(const phf_string_t *a, const phf_string_t *b) {
			int cmp;
			if ((cmp = memcmp(a->p, b->p, PHF_MIN(a->n, b->n))))
				return cmp;
			if (a->n > b->n)
				return -1;
			if (a->n < b->n)
				return 1;
			return 0;
		} /* cmp<phf_string_t>() */
	} /* Uniq:: */
} /* PHF:: */

template<typename key_t>
PHF_PUBLIC size_t PHF::uniq(key_t k[], const size_t n) {
	using namespace PHF::Uniq;
	size_t i, j;

	qsort(k, n, sizeof *k, reinterpret_cast<int(*)(const void *, const void *)>(&cmp<key_t>));

	for (i = 1, j = 0; i < n; i++) {
		if (k[i] != k[j])
			k[++j] = k[i];
	}

	return (n > 0)? j + 1 : 0;
} /* PHF::uniq() */

template size_t PHF::uniq<uint32_t>(uint32_t[], const size_t);
template size_t PHF::uniq<uint64_t>(uint64_t[], const size_t);
template size_t PHF::uniq<phf_string_t>(phf_string_t[], const size_t);
#if !PHF_NO_LIBCXX
template size_t PHF::uniq<std::string>(std::string[], const size_t);
#endif

PHF_PUBLIC size_t phf_uniq_uint32(uint32_t k[], const size_t n) {
	return PHF::uniq(k, n);
} /* phf_uniq_uint32() */

PHF_PUBLIC size_t phf_uniq_uint64(uint64_t k[], const size_t n) {
	return PHF::uniq(k, n);
} /* phf_uniq_uint64() */

PHF_PUBLIC size_t phf_uniq_string(phf_string_t k[], const size_t n) {
	return PHF::uniq(k, n);
} /* phf_uniq_string() */


/*
 * H A S H  P R I M I T I V E S
 *
 * Universal hash based on MurmurHash3_x86_32. Variants for 32- and 64-bit
 * integer keys, and string keys.
 *
 * We use a random seed to address the non-cryptographic-strength collision
 * resistance of MurmurHash3. A stronger hash like SipHash is just too slow
 * and unnecessary for my particular needs. For some environments a
 * cryptographically stronger hash may be prudent.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static inline uint32_t phf_round32(uint32_t k1, uint32_t h1) {
	k1 *= UINT32_C(0xcc9e2d51);
	k1 = PHF_ROTL(k1, 15);
	k1 *= UINT32_C(0x1b873593);

	h1 ^= k1;
	h1 = PHF_ROTL(h1, 13);
	h1 = h1 * 5 + UINT32_C(0xe6546b64);

	return h1;
} /* phf_round32() */

static inline uint32_t phf_round32(const unsigned char *p, size_t n, uint32_t h1) {
	uint32_t k1;

	while (n >= 4) {
		k1 = (p[0] << 24)
		   | (p[1] << 16)
		   | (p[2] << 8)
		   | (p[3] << 0);

		h1 = phf_round32(k1, h1);

		p += 4;
		n -= 4;
	}

	k1 = 0;

	switch (n & 3) {
	case 3:
		k1 |= p[2] << 8;
	case 2:
		k1 |= p[1] << 16;
	case 1:
		k1 |= p[0] << 24;
		h1 = phf_round32(k1, h1);
	}

	return h1;
} /* phf_round32() */

static inline uint32_t phf_round32(phf_string_t k, uint32_t h1) {
	return phf_round32(reinterpret_cast<const unsigned char *>(k.p), k.n, h1);
} /* phf_round32() */

#if !PHF_NO_LIBCXX
static inline uint32_t phf_round32(std::string k, uint32_t h1) {
	return phf_round32(reinterpret_cast<const unsigned char *>(k.c_str()), k.length(), h1);
} /* phf_round32() */
#endif

static inline uint32_t phf_mix32(uint32_t h1) {
	h1 ^= h1 >> 16;
	h1 *= UINT32_C(0x85ebca6b);
	h1 ^= h1 >> 13;
	h1 *= UINT32_C(0xc2b2ae35);
	h1 ^= h1 >> 16;

	return h1;
} /* phf_mix32() */


/*
 * g(k) & f(d, k)  S P E C I A L I Z A T I O N S
 *
 * For every key we first calculate g(k). Then for every group of collisions
 * from g(k) we search for a displacement value d such that f(d, k) places
 * each key into a unique hash slot.
 *
 * g() and f() are specialized for 32-bit, 64-bit, and string keys.
 *
 * g_mod_r() and f_mod_n() are specialized for the method of modular
 * reduction--modular division or bitwise AND. bitwise AND is substantially
 * faster than modular division, and more than makes up for any space
 * inefficiency, particularly for small hash tables.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/* 32-bit, phf_string_t, and std::string keys */
template<typename T>
static inline uint32_t phf_g(T k, uint32_t seed) {
	uint32_t h1 = seed;

	h1 = phf_round32(k, h1);

	return phf_mix32(h1);
} /* phf_g() */

template<typename T>
static inline uint32_t phf_f(uint32_t d, T k, uint32_t seed) {
	uint32_t h1 = seed;

	h1 = phf_round32(d, h1);
	h1 = phf_round32(k, h1);

	return phf_mix32(h1);
} /* phf_f() */


/* 64-bit keys */
static inline uint32_t phf_g(uint64_t k, uint32_t seed) {
	uint32_t h1 = seed;

	h1 = phf_round32(k, h1);
	h1 = phf_round32(k >> 32, h1);

	return phf_mix32(h1);
} /* phf_g() */

static inline uint32_t phf_f(uint32_t d, uint64_t k, uint32_t seed) {
	uint32_t h1 = seed;

	h1 = phf_round32(d, h1);
	h1 = phf_round32(static_cast<uint32_t>(k), h1);
	h1 = phf_round32(static_cast<uint32_t>(k >> 32), h1);

	return phf_mix32(h1);
} /* phf_f() */


/* g() and f() which parameterize modular reduction */
template<bool nodiv, typename T>
static inline uint32_t phf_g_mod_r(T k, uint32_t seed, size_t r) {
	return (nodiv)? (phf_g(k, seed) & (r - 1)) : (phf_g(k, seed) % r);
} /* phf_g_mod_r() */

template<bool nodiv, typename T>
static inline uint32_t phf_f_mod_m(uint32_t d, T k, uint32_t seed, size_t m) {
	return (nodiv)? (phf_f(d, k, seed) & (m - 1)) : (phf_f(d, k, seed) % m);
} /* phf_f_mod_m() */


/*
 * B U C K E T  S O R T I N G  I N T E R F A C E S
 *
 * For every key [0..n) we calculate g(k) % r, where 0 < r <= n, and
 * associate it with a bucket [0..r). We then sort the buckets in decreasing
 * order according to the number of keys. The sorting is required for both
 * optimal time complexity when calculating f(d, k) (less contention) and
 * optimal space complexity (smaller d).
 *
 * The actual sorting is done in the core routine. The buckets are organized
 * and sorted as a 1-dimensional array to minimize run-time memory (less
 * data structure overhead) and improve data locality (less pointer
 * indirection). The following section merely implements a templated
 * bucket-key structure and the comparison routine passed to qsort(3).
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

static bool operator==(const phf_string_t &a, const phf_string_t &b) {
	return a.n == b.n && 0 == memcmp(a.p, b.p, a.n);
}

template<typename T>
class phf_key {
    public:
	T k;
	phf_hash_t g; /* result of g(k) % r */
	size_t *n;  /* number of keys in bucket g */
    phf_key():
        g(0), n(NULL) {
    }
    
}; /* struct phf_key */

template<typename T>
struct  phf_keycmp {

    bool operator()(const phf_key<T>& a, const phf_key<T>& b) {
        if (a.n < b.n) {
            return true;
        }
        if (a.g < b.g) {
            return true;
        }
        /* duplicate key? */
        if ((a.k == b.k) && (&a != &b)) {
            throw hal_exception("duplicate key in perfect hash");
        }
        return false;
    }
};


/*
 * C O R E  F U N C T I O N  G E N E R A T O R
 *
 * The entire algorithm is contained in PHF:init. Everything else in this
 * source file is either a simple utility routine used by PHF:init, or an
 * interface to PHF:init or the generated function state.
 *
 *  - phf - phf object to initialize
 *  - k - array of keys
 *  - n - number of keys
 *  - l - lambda, number of keys per displacement map bucket (reported as g_load)
 *  - a - alpha, hash table load factor (1%% - 100%%).
 *  - seed - random seed
 *  - nodiv - if true rounds r and m to powers of 2
 *
 * Returns a system error number on failure, or 0 on success. phf is
 * unmodified on failure.
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

template<typename key_t, bool nodiv>
PHF_PUBLIC int PHF::init(struct phf *phf, const key_t k[], const size_t n, const size_t l, const size_t a, const phf_seed_t seed) {
	size_t n1 = PHF_MAX(n, 1); /* for computations that require n > 0 */
	size_t l1 = PHF_MAX(l, 1);
	size_t a1 = PHF_MAX(PHF_MIN(a, 100), 1);
	size_t r; /* number of buckets */
	size_t m; /* size of output array */
        // original code used C arrays, but we change to vector to get initialization.
        std::vector<phf_key<key_t> > B_k; /* linear bucket-slot array */
	std::vector<size_t> B_z;         /* number of slots per bucket */
	phf_key<key_t> *B_p, *B_pe;
	phf_bits_t *T = NULL; /* bitmap to track index occupancy */
	phf_bits_t *T_b;      /* per-bucket working bitmap */
	size_t T_n;
	uint32_t *g = NULL; /* displacement map */
	uint32_t d_max = 0; /* maximum displacement value */
	int error;

	if ((phf->nodiv = nodiv)) {
		/* round to power-of-2 so we can use bit masks instead of modulo division */
		r = phf_powerup(n1 / PHF_MIN(l1, n1));
		m = phf_powerup((n1 * 100) / a1);
	} else {
		r = phf_primeup(PHF_HOWMANY(n1, l1));
		/* XXX: should we bother rounding m to prime number for small n? */
		m = phf_primeup((n1 * 100) / a1);
	}

	if (r == 0 || m == 0)
		return ERANGE;

        B_k.resize(n1);
        B_z.resize(r);

	for (size_t i = 0; i < n; i++) {
		phf_hash_t g = phf_g_mod_r<nodiv>(k[i], seed, r);

		B_k[i].k = k[i];
		B_k[i].g = g;
		B_k[i].n = &B_z[g];
		++*B_k[i].n;
	}

        sort(B_k.begin(), B_k.end(), phf_keycmp<key_t>());

	T_n = PHF_HOWMANY(m, PHF_BITS(*T));
	if (!(T = static_cast<phf_bits_t *>(calloc(T_n * 2, sizeof *T))))
		goto syerr;
	T_b = &T[T_n]; /* share single allocation */

	/*
	 * FIXME: T_b[] is unnecessary. We could clear T[] the same way we
	 * clear T_b[]. In fact, at the end of generation T_b[] is identical
	 * to T[] because we don't clear T_b[] on success.
	 *
	 * We just need to tweak the current reset logic to stop before the
	 * key that failed, and then we can elide the commit to T[] at the
	 * end of the outer loop.
	 */

	if (!(g = static_cast<uint32_t *>(calloc(r, sizeof *g))))
		goto syerr;

	B_p = &B_k[0];
	B_pe = &B_k[n];

	for (; B_p < B_pe && *B_p->n > 0; B_p += *B_p->n) {
		phf_key<key_t> *Bi_p, *Bi_pe;
		size_t d = 0;
		uint32_t f;
retry:
		d++;
		Bi_p = B_p;
		Bi_pe = B_p + *B_p->n;

		for (; Bi_p < Bi_pe; Bi_p++) {
			f = phf_f_mod_m<nodiv>(d, Bi_p->k, seed, m);

			if (phf_isset(T, f) || phf_isset(T_b, f)) {
				/* reset T_b[] */
				for (Bi_p = B_p; Bi_p < Bi_pe; Bi_p++) {
					f = phf_f_mod_m<nodiv>(d, Bi_p->k, seed, m);
					phf_clrbit(T_b, f);
				}

				goto retry;
			} else {
				phf_setbit(T_b, f);
			}
		}

		/* commit to T[] */
		for (Bi_p = B_p; Bi_p < Bi_pe; Bi_p++) {
			f = phf_f_mod_m<nodiv>(d, Bi_p->k, seed, m);
			phf_setbit(T, f);
		}

		/* commit to g[] */
		g[B_p->g] = d;
		d_max = PHF_MAX(d, d_max);
	}

	phf->seed = seed;
	phf->r = r;
	phf->m = m;

	phf->g = g;
	g = NULL;

	phf->d_max = d_max;
	phf->g_op = (nodiv)? phf::PHF_G_UINT32_BAND_R : phf::PHF_G_UINT32_MOD_R;
#if 0
	phf->g_jmp = NULL;
#endif

	error = 0;

	goto clean;
syerr:
	error = errno;
clean:
	free(g);
	free(T);

	return error;
} /* PHF::init() */


/*
 * D I S P L A C E M E N T  M A P  C O M P A C T I O N
 *
 * By default the displacement map is an array of uint32_t. This routine
 * compacts the map by using the smallest primitive type that will fit the
 * largest displacement value.
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

template<typename dst_t, typename src_t>
static inline void phf_memmove(dst_t *dst, src_t *src, size_t n) {
	for (size_t i = 0; i < n; i++) {
		dst_t tmp = src[i];
		dst[i] = tmp;
	}
} /* phf_memmove() */

PHF_PUBLIC void PHF::compact(struct phf *phf) {
	size_t size = 0;
	void *tmp;

	switch (phf->g_op) {
	case phf::PHF_G_UINT32_MOD_R:
	case phf::PHF_G_UINT32_BAND_R:
		break;
	default:
		return; /* already compacted */
	}

	if (phf->d_max <= 255) {
		phf_memmove(reinterpret_cast<uint8_t *>(phf->g), reinterpret_cast<uint32_t *>(phf->g), phf->r);
		phf->g_op = (phf->nodiv)? phf::PHF_G_UINT8_BAND_R : phf::PHF_G_UINT8_MOD_R;
		size = sizeof (uint8_t);
	} else if (phf->d_max <= 65535) {
		phf_memmove(reinterpret_cast<uint16_t *>(phf->g), reinterpret_cast<uint32_t *>(phf->g), phf->r);
		phf->g_op = (phf->nodiv)? phf::PHF_G_UINT16_BAND_R : phf::PHF_G_UINT16_MOD_R;
		size = sizeof (uint16_t);
	} else {
		return; /* nothing to compact */
	}

	/* simply keep old array if realloc fails */
	if ((tmp = realloc(phf->g, phf->r * size)))
		phf->g = static_cast<uint32_t *>(tmp);
} /* PHF::compact() */


/*
 * F U N C T I O N  G E N E R A T O R  &  S T A T E  I N T E R F A C E S
 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

template int PHF::init<uint32_t, true>(struct phf *, const uint32_t[], const size_t, const size_t, const size_t, const phf_seed_t);
template int PHF::init<uint64_t, true>(struct phf *, const uint64_t[], const size_t, const size_t, const size_t, const phf_seed_t);
template int PHF::init<phf_string_t, true>(struct phf *, const phf_string_t[], const size_t, const size_t, const size_t, const phf_seed_t);
#if !PHF_NO_LIBCXX
template int PHF::init<std::string, true>(struct phf *, const std::string[], const size_t, const size_t, const size_t, const phf_seed_t);
#endif

template int PHF::init<uint32_t, false>(struct phf *, const uint32_t[], const size_t, const size_t, const size_t, const phf_seed_t);
template int PHF::init<uint64_t, false>(struct phf *, const uint64_t[], const size_t, const size_t, const size_t, const phf_seed_t);
template int PHF::init<phf_string_t, false>(struct phf *, const phf_string_t[], const size_t, const size_t, const size_t, const phf_seed_t);
#if !PHF_NO_LIBCXX
template int PHF::init<std::string, false>(struct phf *, const std::string[], const size_t, const size_t, const size_t, const phf_seed_t);
#endif

template<bool nodiv, typename map_t, typename key_t>
static inline phf_hash_t phf_hash_(map_t *g, key_t k, uint32_t seed, size_t r, size_t m) {
	if (nodiv) {
		uint32_t d = g[phf_g(k, seed) & (r - 1)];

		return phf_f(d, k, seed) & (m - 1);
	} else {
		uint32_t d = g[phf_g(k, seed) % r];

		return phf_f(d, k, seed) % m;
	}
} /* phf_hash_() */

template<typename T>
PHF_PUBLIC phf_hash_t PHF::hash(const struct phf *phf, T k) {
#if PHF_HAVE_COMPUTED_GOTOS && !PHF_NO_COMPUTED_GOTOS
	static const void *const jmp[] = {
		NULL,
		&&uint8_mod_r, &&uint8_band_r,
		&&uint16_mod_r, &&uint16_band_r,
		&&uint32_mod_r, &&uint32_band_r,
	};

	goto *((phf->g_jmp)? phf->g_jmp : (phf->g_jmp = jmp[phf->g_op]));

	uint8_mod_r:
		return phf_hash_<false>(reinterpret_cast<uint8_t *>(phf->g), k, phf->seed, phf->r, phf->m);
	uint8_band_r:
		return phf_hash_<true>(reinterpret_cast<uint8_t *>(phf->g), k, phf->seed, phf->r, phf->m);
	uint16_mod_r:
		return phf_hash_<false>(reinterpret_cast<uint16_t *>(phf->g), k, phf->seed, phf->r, phf->m);
	uint16_band_r:
		return phf_hash_<true>(reinterpret_cast<uint16_t *>(phf->g), k, phf->seed, phf->r, phf->m);
	uint32_mod_r:
		return phf_hash_<false>(reinterpret_cast<uint32_t *>(phf->g), k, phf->seed, phf->r, phf->m);
	uint32_band_r:
		return phf_hash_<true>(reinterpret_cast<uint32_t *>(phf->g), k, phf->seed, phf->r, phf->m);
#else
	switch (phf->g_op) {
	case phf::PHF_G_UINT8_MOD_R:
		return phf_hash_<false>(reinterpret_cast<uint8_t *>(phf->g), k, phf->seed, phf->r, phf->m);
	case phf::PHF_G_UINT8_BAND_R:
		return phf_hash_<true>(reinterpret_cast<uint8_t *>(phf->g), k, phf->seed, phf->r, phf->m);
	case phf::PHF_G_UINT16_MOD_R:
		return phf_hash_<false>(reinterpret_cast<uint16_t *>(phf->g), k, phf->seed, phf->r, phf->m);
	case phf::PHF_G_UINT16_BAND_R:
		return phf_hash_<true>(reinterpret_cast<uint16_t *>(phf->g), k, phf->seed, phf->r, phf->m);
	case phf::PHF_G_UINT32_MOD_R:
		return phf_hash_<false>(reinterpret_cast<uint32_t *>(phf->g), k, phf->seed, phf->r, phf->m);
	case phf::PHF_G_UINT32_BAND_R:
		return phf_hash_<true>(reinterpret_cast<uint32_t *>(phf->g), k, phf->seed, phf->r, phf->m);
	default:
		abort();
		return 0;
	}
#endif
} /* PHF::hash() */

template phf_hash_t PHF::hash<uint32_t>(const struct phf *, uint32_t);
template phf_hash_t PHF::hash<uint64_t>(const struct phf *, uint64_t);
template phf_hash_t PHF::hash<phf_string_t>(const struct phf *, phf_string_t);
#if !PHF_NO_LIBCXX
template phf_hash_t PHF::hash<std::string>(const struct phf *, std::string);
#endif

/** Deallocates internal tables, but not the struct object itself. */
PHF_PUBLIC void PHF::destroy(struct phf *phf) {
	free(phf->g);
	phf->g = NULL;
} /* PHF::destroy() */



