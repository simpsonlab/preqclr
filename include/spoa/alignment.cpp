/*!
 * @file alignment.cpp
 *
 * @brief Alignment, AlignmentParams and AlignmentType classes source file
 */

#include <limits>
#include <algorithm>

#include "node.hpp"
#include "edge.hpp"
#include "graph.hpp"
#include "sisd_alignment.hpp"
#include "simd_alignment.hpp"
#include "alignment.hpp"

namespace SPOA {

AlignmentParams::AlignmentParams(int16_t m, int16_t mm, int16_t gap_opn,
    int16_t gap_ext, AlignmentType t) :
        match(m), mismatch(mm), insertion_open(gap_opn), insertion_extend(gap_ext),
        deletion_open(gap_opn), deletion_extend(gap_ext), type(t) {

    assert(type == AlignmentType::kNW || type == AlignmentType::kSW || type == AlignmentType::kOV);
}

AlignmentParams::AlignmentParams(int16_t m, int16_t mm, int16_t ins_opn,
    int16_t ins_ext, int16_t del_opn, int16_t del_ext, AlignmentType t) :
        match(m), mismatch(mm), insertion_open(ins_opn), insertion_extend(ins_ext),
        deletion_open(del_opn), deletion_extend(del_ext), type(t) {

    assert(type == AlignmentType::kNW || type == AlignmentType::kSW || type == AlignmentType::kOV);
}

AlignmentParams::~AlignmentParams() {
}

std::unique_ptr<Alignment> Alignment::createAlignment(const std::string& sequence,
    std::shared_ptr<Graph> graph, AlignmentParams params) {

    auto a = createSimdAlignment(sequence, graph, params);
    if (!a) {
        a = createSisdAlignment(sequence, graph, params);
    }

    return a;
}

}
