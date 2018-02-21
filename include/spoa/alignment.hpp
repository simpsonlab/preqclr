/*!
 * @file alignment.hpp
 *
 * @brief Alignment, AlignmentParams and AlignmentType classes header file
 */

#pragma once

#include <assert.h>
#include <string>
#include <memory>
#include <vector>

namespace SPOA {

enum class AlignmentType {
    kSW, // Smith Waterman
    kNW, // Needleman Wunsch
    kOV // Overlap
};

class AlignmentParams {
public:

    AlignmentParams(int16_t match, int16_t mismatch, int16_t gap_open,
        int16_t gap_extend, AlignmentType type);
    AlignmentParams(int16_t match, int16_t mismatch, int16_t insertion_open,
        int16_t insertion_extend, int16_t deletion_open, int16_t deletion_extend,
        AlignmentType type);
    ~AlignmentParams();

    int16_t match;
    int16_t mismatch;
    int16_t insertion_open;
    int16_t insertion_extend;
    int16_t deletion_open;
    int16_t deletion_extend;
    AlignmentType type;
};

class Graph;

class Alignment {
public:

    static std::unique_ptr<Alignment> createAlignment(const std::string& sequence,
        std::shared_ptr<Graph> graph, AlignmentParams params);

    virtual ~Alignment() {}

    virtual const std::vector<int32_t>& seq_ids() const = 0;
    virtual const std::vector<int32_t>& node_ids() const = 0;

    virtual void align_sequence_to_graph() = 0;
    virtual void adjust_node_ids(const std::vector<int32_t>& mapping) = 0;
};

}
