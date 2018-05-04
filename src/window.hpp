/*!
 * @file window.hpp
 *
 * @brief Window class header file
 */

#pragma once

#include <stdlib.h>
#include <vector>
#include <memory>
#include <string>
#include <utility>
#include <map>

namespace spoa {
    class AlignmentEngine;
}

namespace racon {

enum class WindowType {
    kNGS, // Next Generation Sequencing
    kTGS // Third Generation Sequencing
};

class Window;
std::unique_ptr<Window> createWindow(uint64_t id, uint32_t rank, WindowType type,
    const char* backbone, uint32_t backbone_length, const char* quality,
    uint32_t quality_length);

class Window {
public:
    ~Window();

    uint64_t id() const {
        return id_;
    }
    uint32_t rank() const {
        return rank_;
    }

    const std::map<float,int>& allele_ratio() const {
        return allele_ratio_;
    }

    const std::map<int,std::vector<std::string>>& second_msas() const {
        return second_msas_;
    }

    const std::string& consensus() const {
        return consensus_;
    }
    
    const std::string& msa_consensus() const {
        return msa_consensus_;
    }   
  
    std::map<float,int> allele_ratio_from_msa(std::vector<std::string> &msa, const char * depth, const char * percent_gaps);   

    bool generate_consensus(std::shared_ptr<spoa::AlignmentEngine> alignment_engine, int8_t min_spoa_coverage, int8_t allowed_spoa_gaps_percent);

    void add_layer(const char* sequence, uint32_t sequence_length,
        const char* quality, uint32_t quality_length, uint32_t begin,
        uint32_t end);

    friend std::unique_ptr<Window> createWindow(uint64_t id, uint32_t rank,
        WindowType type, const char* backbone, uint32_t backbone_length,
        const char* quality, uint32_t quality_length);
private:
    Window(uint64_t id, uint32_t rank, WindowType type, const char* backbone,
        uint32_t backbone_length, const char* quality, uint32_t quality_length);
    Window(const Window&) = delete;
    const Window& operator=(const Window&) = delete;

    uint64_t id_;
    uint32_t rank_;
    WindowType type_;
    std::string consensus_;
    std::string msa_consensus_;
    std::map<float,int> allele_ratio_;
    std::map<int,std::vector<std::string>> second_msas_;
    std::vector<std::pair<const char*, uint32_t>> sequences_;
    std::vector<std::pair<const char*, uint32_t>> qualities_;
    std::vector<std::pair<uint32_t, uint32_t>> positions_;
    
};

}
