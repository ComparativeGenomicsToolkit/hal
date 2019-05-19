#ifndef _MMAP_METADATA_H
#define _MMAP_METADATA_H
#include "halMetaData.h"
#include "mmapArray.h"
#include "mmapString.h"

namespace hal {

    // yes. this is a stupid name. but consistency above all!
    class MMapMetaDataData {
      public:
        size_t _keysOffset;
        size_t _valuesOffset;
        char _reserved[256];   // 256 bytes of reserved added in mmap API 1.1
    };

    class MMapMetaData : public MetaData {
      public:
        MMapMetaData(MMapAlignment *alignment)
            : _alignment(alignment), _offset(_alignment->allocateNewArray(sizeof(MMapMetaDataData))),
              _data((MMapMetaDataData *)_alignment->resolveOffset(_offset, sizeof(MMapMetaDataData))), _keys(alignment),
              _values(alignment), _dirty(true) {
            _data->_keysOffset = _keys.getOffset();
            _data->_valuesOffset = _values.getOffset();
        };
        MMapMetaData(MMapAlignment *alignment, size_t offset)
            : _alignment(alignment), _offset(offset),
              _data((MMapMetaDataData *)_alignment->resolveOffset(_offset, sizeof(MMapMetaDataData))),
              _keys(_alignment, _data->_keysOffset), _values(_alignment, _data->_valuesOffset), _dirty(false) {
            read();
        };
        ~MMapMetaData() {
            if (_dirty) {
                write();
            }
        };

        void set(const std::string &key, const std::string &value) {
            _map.insert(make_pair(key, value));
            _dirty = true;
        };
        const std::string &get(const std::string &key) const {
            return _map.find(key)->second;
        }
        bool has(const std::string &key) const {
            return _map.count(key) > 0;
        }
        const std::map<std::string, std::string> &getMap() const {
            return _map;
        }
        size_t getOffset() const {
            return _offset;
        };

      private:
        std::map<std::string, std::string> _map;
        MMapAlignment *_alignment;
        size_t _offset;
        MMapMetaDataData *_data;
        // Arrays of offsets to strings containing keys and values.
        MMapArray<size_t> _keys;
        MMapArray<size_t> _values;
        bool _dirty;
        // Update _map to reflect the file.
        void read() {
            if (_keys.getLength() != _values.getLength()) {
                throw hal_exception("# keys != # values in metadata");
            }
            for (size_t i = 0; i < _keys.getLength(); i++) {
                MMapString key{_alignment, *_keys[i]};
                MMapString value{_alignment, *_values[i]};
                _map.insert(make_pair(key.get(), value.get()));
            }
        }
        // Update the file to reflect _map.
        void write() {
            _data->_keysOffset = _keys.setLength(_map.size());
            _data->_valuesOffset = _values.setLength(_map.size());
            size_t i = 0;
            for (auto kv : _map) {
                MMapString key{_alignment, *_keys[i]};
                MMapString value{_alignment, *_values[i]};
                *_keys[i] = key.set(kv.first);
                *_values[i] = value.set(kv.second);
                i++;
            }
            _dirty = false;
        };
    };
}
#endif
// Local Variables:
// mode: c++
// End:
