/*!
 * @file chain.hpp
 *
 * @brief Chain class header file
 */

#ifdef SPOA_MAIN_

#pragma once

#include <stdlib.h>
#include <stdint.h>
#include <memory>
#include <vector>
#include <string>

#include "bioparser/src/bioparser.hpp"

namespace SPOA {

class Chain;

std::unique_ptr<Chain> createChain(uint64_t id, const char* name, uint32_t name_length,
    const char* data, uint32_t data_length);

std::unique_ptr<Chain> createChain(uint64_t id, const char* name, uint32_t name_length,
    const char* data, uint32_t data_length, const char* quality, uint32_t quality_length);

class Chain {
public:

    ~Chain() = default;

    uint64_t id() const {
        return id_;
    }

    const std::string& name() const {
        return name_;
    }

    const std::string& data() const {
        return data_;
    }

    const std::string& quality() const {
        return quality_;
    }

    friend std::unique_ptr<Chain> createChain(uint64_t id, const char* name,
        uint32_t name_length, const char* data, uint32_t data_length);
    friend std::unique_ptr<Chain> createChain(uint64_t id, const char* name,
        uint32_t name_length, const char* data, uint32_t data_length,
        const char* quality, uint32_t quality_length);

    friend BIOPARSER::FastaReader<Chain>;
    friend BIOPARSER::FastqReader<Chain>;

private:

    Chain(uint64_t id, const char* name, uint32_t name_length, const char* data, uint32_t data_length);
    Chain(uint64_t id, const char* name, uint32_t name_length, const char* data, uint32_t data_length,
        const char* quality, uint32_t quality_length);
    Chain(const Chain&) = delete;
    const Chain& operator=(const Chain&) = delete;

    uint64_t id_;
    std::string name_;
    std::string data_;
    std::string quality_;
};

}

#endif
