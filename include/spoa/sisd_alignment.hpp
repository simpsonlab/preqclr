/*!
 * @file sisd_alignment.hpp
 *
 * @brief SisdAlignment class header file
 */

#pragma once

#include <assert.h>
#include <string>
#include <memory>
#include <vector>

#include "alignment.hpp"

namespace SPOA {

class Graph;
class SisdAlignment;

std::unique_ptr<Alignment> createSisdAlignment(const std::string& sequence,
    std::shared_ptr<Graph> graph, AlignmentParams params);

class SisdAlignment: public Alignment {
public:

    ~SisdAlignment();

    int32_t score() const {
        return max_score_;
    }

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

    friend std::unique_ptr<Alignment> createSisdAlignment(const std::string& sequence,
        std::shared_ptr<Graph> graph, AlignmentParams params);

private:

    SisdAlignment(const std::string& sequence, std::shared_ptr<Graph> graph,
        AlignmentParams params);
    SisdAlignment(const SisdAlignment&) = delete;
    const SisdAlignment& operator=(const SisdAlignment&) = delete;

    inline void update_max_score(int32_t* H_row, uint32_t i, uint32_t j);

    void backtrack();
    void print_matrix();

    std::shared_ptr<Graph> graph_;
    AlignmentParams params_;

    uint32_t matrix_width_;
    uint32_t matrix_height_;

    std::vector<std::vector<int32_t>> sequence_profile_;
    std::vector<int32_t> H_;
    std::vector<int32_t> F_;
    std::vector<int32_t> E_;
    std::vector<uint32_t> node_id_to_graph_id_;

    bool is_aligned_;
    int32_t max_i_;
    int32_t max_j_;
    int32_t max_score_;

    std::vector<int32_t> alignment_node_ids_;
    std::vector<int32_t> alignment_seq_ids_;
};

}
