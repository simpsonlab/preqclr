/*!
 * @file chain.cpp
 *
 * @brief Chain class source file
 */

#ifdef SPOA_MAIN_

#include <assert.h>

#include "chain.hpp"

namespace SPOA {

std::unique_ptr<Chain> createChain(uint64_t id, const char* name, uint32_t name_length,
    const char* data, uint32_t data_length) {

    assert(name_length);
    assert(data_length);

    return std::unique_ptr<Chain>(new Chain(id, name, name_length, data, data_length));
}

std::unique_ptr<Chain> createChain(uint64_t id, const char* name, uint32_t name_length,
    const char* data, uint32_t data_length, const char* quality, uint32_t quality_length) {

    assert(name_length);
    assert(data_length);
    assert(quality_length && data_length == quality_length);

    return std::unique_ptr<Chain>(new Chain(id, name, name_length, data, data_length,
        quality, quality_length));
}

Chain::Chain(uint64_t id, const char* name, uint32_t name_length, const char* data,
    uint32_t data_length)
        : id_(id), name_(name, name_length), data_(data, data_length), quality_() {
}

Chain::Chain(uint64_t id, const char* name, uint32_t name_length, const char* data,
    uint32_t data_length, const char* quality, uint32_t quality_length)
        : id_(id), name_(name, name_length), data_(data, data_length), quality_(quality, quality_length) {
}

}

#endif
