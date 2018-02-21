/*!
 * @file simd_alignment.hpp
 *
 * @brief SimdAlignment class header file
 */

#pragma once

#include <assert.h>
#include <string>
#include <memory>
#include <vector>

#include "alignment.hpp"

namespace SPOA {

class Graph;
class SimdAlignment;

std::unique_ptr<Alignment> createSimdAlignment(const std::string& sequence,
    std::shared_ptr<Graph> graph, AlignmentParams params);

class SimdAlignment: public Alignment {
public:

    ~SimdAlignment();

    const std::vector<int32_t>& node_ids() const {
        assert(is_aligned_ == true && "No alignment done!");
        return alignment_node_ids_;
    }

    const std::vector<int32_t>& seq_ids() const {
        assert(is_aligned_ == true && "No alignment done!");
        return alignment_seq_ids_;
    }

    void align_sequence_to_graph();
    void adjust_node_ids(const std::vector<int32_t>& mapping);

    friend std::unique_ptr<Alignment> createSimdAlignment(const std::string& sequence,
        std::shared_ptr<Graph> graph, AlignmentParams params);

private:

    SimdAlignment(const std::string& sequence, std::shared_ptr<Graph> graph,
        AlignmentParams params);
    SimdAlignment(const SimdAlignment&) = delete;
    const SimdAlignment& operator=(const SimdAlignment&) = delete;

    std::shared_ptr<Graph> graph_;
    AlignmentParams params_;

    uint32_t matrix_width_;
    uint32_t matrix_height_;

    std::vector<std::vector<int32_t>> sequence_profile_;

    bool is_aligned_;
    std::vector<int32_t> alignment_node_ids_;
    std::vector<int32_t> alignment_seq_ids_;
};

void simd_align_sequence_to_graph(const std::string& sequence, std::shared_ptr<Graph> graph,
    AlignmentParams params);

}
