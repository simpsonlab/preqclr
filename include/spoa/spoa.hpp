/*!
 * @file poa.hpp
 *
 * @brief Poa header file which encapsulates the implementation
 */

#pragma once

#include <string>
#include <vector>

#include "alignment.hpp"
#include "graph.hpp"

namespace SPOA {

std::shared_ptr<Graph> construct_partial_order_graph(const std::vector<std::string>& sequences,
    AlignmentParams params, bool sorted = false);

std::shared_ptr<Graph> construct_partial_order_graph(const std::vector<std::string>& sequences,
    const std::vector<std::string>& qualities, AlignmentParams params, bool sorted = false);

std::shared_ptr<Graph> construct_partial_order_graph(const std::vector<std::string>& sequences,
    const std::vector<std::string>& qualities, const std::vector<uint32_t>& begin_positions,
    const std::vector<uint32_t>& end_positions, AlignmentParams params);

std::string generate_consensus(const std::vector<std::string>& sequences,
    AlignmentParams params, bool sorted = false);

std::string generate_consensus(const std::vector<std::string>& sequences,
    const std::vector<std::string>& qualities, AlignmentParams params, bool sorted = false);

std::string generate_consensus(const std::vector<std::string>& sequences,
    const std::vector<std::string>& qualities, const std::vector<uint32_t>& begin_positions,
    const std::vector<uint32_t>& end_positions, AlignmentParams params);

void generate_msa(std::vector<std::string>& dst, const std::vector<std::string>& sequences,
    AlignmentParams params, bool sorted = false);

void generate_msa(std::vector<std::string>& dst, const std::vector<std::string>& sequences,
    const std::vector<std::string>& qualities, AlignmentParams params, bool sorted = false);

}
