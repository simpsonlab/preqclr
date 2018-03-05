/*!
 * @file graph.cpp
 *
 * @brief Graph class source file
 */

#include <algorithm>
#include <stack>
#include <stdio.h>
#include <iostream>
#include "node.hpp"
#include "edge.hpp"
#include "alignment.hpp"
#include "graph.hpp"

namespace SPOA {

std::unique_ptr<Graph> createGraph(const std::string& sequence, float weight) {
    std::vector<float> weights(sequence.size(), weight);
    return createGraph(sequence, weights);
}

std::unique_ptr<Graph> createGraph(const std::string& sequence, const std::string& quality) {
    std::vector<float> weights;
    for (const auto& q: quality) {
        weights.emplace_back((float) (q - 33)); // PHRED quality
    }
    return createGraph(sequence, weights);
}

std::unique_ptr<Graph> createGraph(const std::string& sequence, const std::vector<float>& weights) {
    assert(sequence.size() != 0);
    assert(sequence.size() == weights.size());
    return std::unique_ptr<Graph>(new Graph(sequence, weights));
}

Graph::Graph() :
        num_sequences_(), num_nodes_(), nodes_(), alphabet_(), is_sorted_(false),
        sorted_nodes_ids_(), sequences_start_nodes_ids_(), consensus_() {
}
Graph::Graph(const std::string& sequence, const std::vector<float>& weights) :
        num_sequences_(), num_nodes_(), nodes_(), alphabet_(), is_sorted_(false),
        sorted_nodes_ids_(), sequences_start_nodes_ids_(), consensus_() {

    for (const auto& c: sequence) {
        alphabet_.insert(c);
    }

    int32_t start_node_id = this->add_sequence(sequence, weights, 0, sequence.size());

    sequences_start_nodes_ids_.emplace_back(start_node_id);
    ++num_sequences_;
}

Graph::~Graph() {
}

uint32_t Graph::add_node(char letter) {
    nodes_.emplace_back(createNode(num_nodes_, letter));
    return num_nodes_++;
}

void Graph::add_edge(uint32_t begin_node_id, uint32_t end_node_id, float weight) {
    assert(begin_node_id < num_nodes_ && end_node_id < num_nodes_);

    for (const auto& edge: nodes_[begin_node_id]->out_edges()) {
        if (edge->end_node_id() == end_node_id) {
            edge->add_sequence(num_sequences_, weight);
            return;
        }
    }

    std::shared_ptr<Edge> edge = createEdge(begin_node_id, end_node_id, num_sequences_, weight);
    nodes_[begin_node_id]->add_out_edge(edge);
    nodes_[end_node_id]->add_in_edge(edge);
}

void Graph::topological_sort(bool rigorous) {

    if (is_sorted_) {
        return;
    }
    sorted_nodes_ids_.clear();

    // 0 - unmarked, 1 - temporarily marked, 2 - permanently marked
    std::vector<uint8_t> marks(num_nodes_, 0);
    // false - do not check aligned nodes, true - do it
    std::vector<bool> check(num_nodes_, true);
    std::stack<uint32_t> nodes_to_visit;

    for (uint32_t i = 0; i < num_nodes_; ++i) {
        if (marks[i] != 0) {
            continue;
        }

        nodes_to_visit.push(i);
        while (nodes_to_visit.size() != 0) {
            uint32_t node_id = nodes_to_visit.top();
            const auto& node = nodes_[node_id];
            bool valid = true;

            if (marks[node_id] != 2) {
                for (const auto& edge: node->in_edges()) {
                    if (marks[edge->begin_node_id()] != 2) {
                        nodes_to_visit.push(edge->begin_node_id());
                        valid = false;
                    }
                }

                if (rigorous && check[node_id]) {
                    for (const auto& aid: node->aligned_nodes_ids()) {
                        if (marks[aid] != 2) {
                            nodes_to_visit.push(aid);
                            check[aid] = false;
                            valid = false;
                        }
                    }
                }

                assert((valid || marks[node_id] != 1) && "Graph is not a DAG!");
                if (valid) {
                    marks[node_id] = 2;
                    if (!rigorous) {
                        sorted_nodes_ids_.push_back(node_id);
                    } else if (check[node_id]) {
                        sorted_nodes_ids_.push_back(node_id);
                        for (const auto& aid: node->aligned_nodes_ids()) {
                            sorted_nodes_ids_.emplace_back(aid);
                        }
                    }
                } else {
                    marks[node_id] = 1;
                }
            }

            if (valid) {
                nodes_to_visit.pop();
            }
        }
    }

    assert(this->is_topologically_sorted() == true);
    is_sorted_ = true;
}

bool Graph::is_topologically_sorted() const {
    assert(nodes_.size() == sorted_nodes_ids_.size());

    std::vector<bool> visited_nodes(num_nodes_, false);
    for (uint32_t node_id: sorted_nodes_ids_) {
        for (const auto& edge: nodes_[node_id]->in_edges()) {
            if (visited_nodes[edge->begin_node_id()] == false) {
                return false;
            }
        }
        visited_nodes[node_id] = true;
    }

    return true;
}

// backtracing from right to left!
void Graph::extract_subgraph_nodes(std::vector<bool>& dst, uint32_t start_node_id,
    uint32_t end_node_id) const {

    dst.resize(num_nodes_, false);

    std::stack<uint32_t> nodes_to_visit;
    nodes_to_visit.push(start_node_id);

    while (nodes_to_visit.size() != 0) {
        uint32_t node_id = nodes_to_visit.top();
        nodes_to_visit.pop();

        if (dst[node_id] == false && node_id >= end_node_id) {
            for (const auto& edge: nodes_[node_id]->in_edges()) {
                nodes_to_visit.push(edge->begin_node_id());
            }
            for (const auto& aid: nodes_[node_id]->aligned_nodes_ids()) {
                nodes_to_visit.push(aid);
            }

            dst[node_id] = true;
        }
    }
}

std::unique_ptr<Graph> Graph::subgraph(uint32_t begin_node_id, uint32_t end_node_id,
    std::vector<int32_t>& subgraph_to_graph_mapping) {

    std::vector<bool> is_subgraph_node;
    extract_subgraph_nodes(is_subgraph_node, end_node_id, begin_node_id);

    // init subgraph
    auto subgraph = std::unique_ptr<Graph>(new Graph());
    subgraph->alphabet_ = std::unordered_set<uint8_t>(this->alphabet_);

    // create mapping from subgraph to graph and vice versa
    // and add nodes to subgraph
    subgraph_to_graph_mapping.resize(this->num_nodes_, -1);
    std::vector<int32_t> graph_to_subgraph_mapping(this->num_nodes_, -1);

    for (uint32_t i = 0; i < is_subgraph_node.size(); ++i) {
        if (is_subgraph_node[i] == false) {
            continue;
        }

        uint32_t subgraph_id = subgraph->add_node(nodes_[i]->letter());
        graph_to_subgraph_mapping[i] = subgraph_id;
        subgraph_to_graph_mapping[subgraph_id] = i;
    }

    // add edges and aligned nodes
    for (uint32_t i = 0; i < is_subgraph_node.size(); ++i) {
        if (is_subgraph_node[i] == false) {
            continue;
        }
        const auto& node = nodes_[i];
        uint32_t subgraph_id = graph_to_subgraph_mapping[i];
        const auto& subgraph_node = subgraph->node(subgraph_id);

        for (const auto& edge: node->in_edges()) {
            if (graph_to_subgraph_mapping[edge->begin_node_id()] == -1) continue;
            subgraph->add_edge(graph_to_subgraph_mapping[edge->begin_node_id()], subgraph_id, edge->total_weight());
        }
        for (const auto& aid: node->aligned_nodes_ids()) {
            if (graph_to_subgraph_mapping[aid] == -1) continue;
            subgraph_node->add_aligned_node_id(graph_to_subgraph_mapping[aid]);
        }
    }

    return subgraph;
}

void Graph::add_alignment(std::shared_ptr<Alignment> alignment, const std::string& sequence,
    float weight) {

    std::vector<float> weights(sequence.size(), weight);
    this->add_alignment(alignment, sequence, weights);
}

void Graph::add_alignment(std::shared_ptr<Alignment> alignment, const std::string& sequence,
    const std::string& quality) {

    std::vector<float> weights;
    for (const auto& q: quality) {
        weights.emplace_back((float) (q - 33)); // PHRED quality
    }
    this->add_alignment(alignment, sequence, weights);
}

void Graph::add_alignment(std::shared_ptr<Alignment> alignment, const std::string& sequence,
    const std::vector<float>& weights) {

    assert(sequence.size() != 0);
    assert(sequence.size() == weights.size());

    const auto& node_ids = alignment->node_ids();
    const auto& seq_ids = alignment->seq_ids();

    assert(node_ids.size() == seq_ids.size());

    for (const auto& c: sequence) {
        alphabet_.insert(c);
    }

    if (seq_ids.size() == 0) { // no local alignment!
        int32_t start_node_id = this->add_sequence(sequence, weights, 0, sequence.size());
        ++num_sequences_;
        sequences_start_nodes_ids_.emplace_back(start_node_id);

        is_sorted_ = false;
        return;
    }

    std::vector<uint32_t> valid_seq_ids;
    for (const auto& seq_id: seq_ids) {
        if (seq_id != -1) {
            valid_seq_ids.emplace_back(seq_id);
        }
    }

    uint32_t tmp = num_nodes_;
    int32_t start_node_id = this->add_sequence(sequence, weights, 0, valid_seq_ids.front());
    int32_t head_node_id = tmp == num_nodes_ ? -1 : num_nodes_ - 1;

    int32_t tail_node_id = this->add_sequence(sequence, weights, valid_seq_ids.back() + 1, sequence.size());

    int32_t new_node_id = -1;
    float prev_weight = head_node_id == -1 ? 0 : weights[valid_seq_ids.front() - 1];

    for (uint32_t i = 0; i < seq_ids.size(); ++i) {
        if (seq_ids[i] == -1) {
            continue;
        }

        char letter = sequence[seq_ids[i]];
        if (node_ids[i] == -1) {
            new_node_id = this->add_node(letter);

        } else {
            auto node = nodes_[node_ids[i]];
            if (node->letter() == letter) {
                new_node_id = node_ids[i];

            } else {
                int32_t aligned_to_node_id = -1;
                for (const auto& aid: node->aligned_nodes_ids()) {
                    if (nodes_[aid]->letter() == letter) {
                        aligned_to_node_id = aid;
                        break;
                    }
                }

                if (aligned_to_node_id == -1) {
                    new_node_id = this->add_node(letter);

                    for (const auto& aid: node->aligned_nodes_ids()) {
                        nodes_[new_node_id]->add_aligned_node_id(aid);
                        nodes_[aid]->add_aligned_node_id(new_node_id);
                    }

                    nodes_[new_node_id]->add_aligned_node_id(node_ids[i]);
                    node->add_aligned_node_id(new_node_id);

                } else {
                    new_node_id = aligned_to_node_id;
                }
            }
        }

        if (start_node_id == -1) {
            start_node_id = new_node_id;
        }

        if (head_node_id != -1) {
            this->add_edge(head_node_id, new_node_id, prev_weight + weights[seq_ids[i]]); // both nodes contribute to edge weight
        }

        head_node_id = new_node_id;
        prev_weight = weights[seq_ids[i]];
    }

    if (tail_node_id != -1) {
        this->add_edge(head_node_id, tail_node_id, prev_weight + weights[valid_seq_ids.back() + 1]); // both nodes contribute to edge weight
    }

    ++num_sequences_;
    sequences_start_nodes_ids_.emplace_back(start_node_id);

    is_sorted_ = false;
}

int32_t Graph::add_sequence(const std::string& sequence, const std::vector<float>& weights,
    uint32_t begin, uint32_t end) {

    if (begin == end) {
        return -1;
    }

    assert(begin < sequence.size() && end <= sequence.size());

    int32_t first_node_id = this->add_node(sequence[begin]);

    uint32_t node_id;
    for (uint32_t i = begin + 1; i < end; ++i) {
        node_id = this->add_node(sequence[i]);
        this->add_edge(node_id - 1, node_id, weights[i - 1] + weights[i]); // both nodes contribute to edge weight
    }

    return first_node_id;
}

void Graph::generate_msa(std::vector<std::string>& dst, bool include_consensus) {

    // force rigorous topological sort
    is_sorted_ = false;
    this->topological_sort(true);

    // assign msa id to each node
    std::vector<int32_t> msa_node_ids(num_nodes_, -1);
    int32_t base_counter = 0;
    for (uint32_t i = 0; i < num_nodes_; ++i) {
        uint32_t node_id = sorted_nodes_ids_[i];

        msa_node_ids[node_id] = base_counter;
        for (uint32_t j = 0; j < nodes_[node_id]->aligned_nodes_ids().size(); ++j) {
            msa_node_ids[sorted_nodes_ids_[++i]] = base_counter;
        }
        ++base_counter;
    }

    // extract sequences from graph and create msa strings (add indels(-) where necessary)
    for (uint32_t i = 0; i < num_sequences_; ++i) {
        std::string alignment_str(base_counter, '-');
        uint32_t curr_node_id = sequences_start_nodes_ids_[i];

        while (true) {
            alignment_str[msa_node_ids[curr_node_id]] = nodes_[curr_node_id]->letter();

            uint32_t prev_node_id = curr_node_id;
            for (const auto& edge: nodes_[prev_node_id]->out_edges()) {
                for (const auto& label: edge->sequence_labels()) {
                    if (label == i) {
                        curr_node_id = edge->end_node_id();
                        break;
                    }
                }
                if (prev_node_id != curr_node_id) {
                    break;
                }
            }

            if (prev_node_id == curr_node_id) {
                break;
            }
        }

        dst.emplace_back(alignment_str);
    }

    if (include_consensus) {
        // do the same for consensus sequence
        this->traverse_heaviest_bundle();

        std::string alignment_str(base_counter, '-');
        for (const auto& node_id: consensus_) {
            alignment_str[msa_node_ids[node_id]] = nodes_[node_id]->letter();
        }
        dst.emplace_back(alignment_str);
    }
}

void Graph::check_msa(const std::vector<std::string>& msa, const std::vector<std::string>& sequences,
    const std::vector<uint32_t>& indices) const {

    for (uint32_t i = 0; i < sequences.size(); ++i) {
        std::string temp = "";
        for (const auto& c: msa[i]) {
            if (c != '-') temp += c;
        }
        assert(temp.size() == sequences[indices[i]].size() && "different lenghts");
        assert(temp.compare(sequences[indices[i]]) == 0 && "different sequence");
    }
}

std::string Graph::generate_consensus() {

    this->traverse_heaviest_bundle();
    std::string consensus_str = "";
    for (const auto& node_id: consensus_) {
        consensus_str += nodes_[node_id]->letter();
    }

    return consensus_str;
}

std::string Graph::generate_consensus(std::vector<uint32_t>& dst) {

    auto consensus_str = this->generate_consensus();

    auto calculate_coverage = [&](uint32_t node_id) -> uint32_t {
        std::unordered_set<uint32_t> label_set;
        const auto& node = nodes_[node_id];
        for (const auto& edge: node->in_edges()) {
            for (const auto& label: edge->sequence_labels()) {
                label_set.insert(label);
            }
        }
        for (const auto& edge: node->out_edges()) {
            for (const auto& label: edge->sequence_labels()) {
                label_set.insert(label);
            }
        }
        return label_set.size();
    };

    for (const auto& node_id: consensus_) {
        uint32_t total_coverage = calculate_coverage(node_id);
        for (const auto& aid: nodes_[node_id]->aligned_nodes_ids()) {
            total_coverage += calculate_coverage(aid);
        }
        dst.emplace_back(total_coverage);
    }

    return consensus_str;
}

void Graph::traverse_heaviest_bundle() {

    this->topological_sort();

    std::vector<int32_t> predecessors(num_nodes_, -1);
    std::vector<float> scores(num_nodes_, 0);

    uint32_t max_score_id = 0;
    for (const auto& node_id: sorted_nodes_ids_) {
        for (const auto& edge: nodes_[node_id]->in_edges()) {
            if (scores[node_id] < edge->total_weight() ||
                (scores[node_id] == edge->total_weight() &&
                scores[predecessors[node_id]] <= scores[edge->begin_node_id()])) {

                scores[node_id] = edge->total_weight();
                predecessors[node_id] = edge->begin_node_id();
            }
        }

        if (predecessors[node_id] != -1) {
            scores[node_id] += scores[predecessors[node_id]];
        }

        if (scores[max_score_id] < scores[node_id]) {
            max_score_id = node_id;
        }
    }

    if (nodes_[max_score_id]->out_edges().size() != 0) {

        std::vector<uint32_t> node_id_to_rank(num_nodes_, 0);
        for (uint32_t i = 0; i < num_nodes_; ++i) {
            node_id_to_rank[sorted_nodes_ids_[i]] = i;
        }

        while (nodes_[max_score_id]->out_edges().size() != 0) {
            max_score_id = this->branch_completion(scores, predecessors,
                node_id_to_rank[max_score_id]);
        }
    }

    // traceback
    consensus_.clear();
    while (predecessors[max_score_id] != -1) {
        consensus_.emplace_back(max_score_id);
        max_score_id = predecessors[max_score_id];
    }
    consensus_.emplace_back(max_score_id);

    std::reverse(consensus_.begin(), consensus_.end());
}

uint32_t Graph::branch_completion(std::vector<float>& scores,
    std::vector<int32_t>& predecessors, uint32_t rank) {

    uint32_t node_id = sorted_nodes_ids_[rank];
    for (const auto& edge: nodes_[node_id]->out_edges()) {
        for (const auto& o_edge: nodes_[edge->end_node_id()]->in_edges()) {
            if (o_edge->begin_node_id() != node_id) {
                scores[o_edge->begin_node_id()] = -1;
            }
        }
    }

    float max_score = 0;
    uint32_t max_score_id = 0;
    for (uint32_t i = rank + 1; i < sorted_nodes_ids_.size(); ++i) {

        uint32_t node_id = sorted_nodes_ids_[i];
        scores[node_id] = -1;
        predecessors[node_id] = -1;

        for (const auto& edge: nodes_[node_id]->in_edges()) {
            if (scores[edge->begin_node_id()] == -1) {
                continue;
            }

            if (scores[node_id] < edge->total_weight() ||
                (scores[node_id] == edge->total_weight() &&
                scores[predecessors[node_id]] <= scores[edge->begin_node_id()])) {

                scores[node_id] = edge->total_weight();
                predecessors[node_id] = edge->begin_node_id();
            }
        }

        if (predecessors[node_id] != -1) {
            scores[node_id] += scores[predecessors[node_id]];
        }

        if (max_score < scores[node_id]) {
            max_score = scores[node_id];
            max_score_id = node_id;
        }
    }

    return max_score_id;
}

void Graph::print() const {

    std::vector<int32_t> in_consensus(num_nodes_, -1);
    int32_t rank = 0;
    for (const auto& id: consensus_) {
        in_consensus[id] = rank++;
    }

    printf("digraph %d {\n", num_sequences_);
    printf("    graph [rankdir=LR]\n");
    for (uint32_t i = 0; i < num_nodes_; ++i) {
        printf("    %d [label = \"%d|%c\"", i, i, nodes_[i]->letter());
        if (in_consensus[i] != -1) {
            printf(", style=filled, fillcolor=goldenrod1");
        }
        printf("]\n");

        for (const auto& edge: nodes_[i]->out_edges()) {
            printf("    %d -> %d [label = \"%.3f\"", i, edge->end_node_id(),
                edge->total_weight());
            if (in_consensus[i] + 1 == in_consensus[edge->end_node_id()]) {
                printf(", color=goldenrod1");
            }
            printf("]\n");
        }
        for (const auto& aid: nodes_[i]->aligned_nodes_ids()) {
            if (aid > i) {
                printf("    %d -> %d [style = dotted, arrowhead = none]\n", i, aid);
            }
        }
    }
    printf("}\n");
}

void Graph::printtofile() const {

    std::cout <<"Inside printtofile" <<"\n";
    FILE * pFile;
    pFile = fopen ("spoa.dot","w");
    std::vector<int32_t> in_consensus(num_nodes_, -1);
    int32_t rank = 0;
    for (const auto& id: consensus_) {
        in_consensus[id] = rank++;
    }

    fprintf(pFile, "digraph %d {\n", num_sequences_);
    fprintf(pFile, "    graph [rankdir=LR]\n");
    for (uint32_t i = 0; i < num_nodes_; ++i) {
        fprintf(pFile, "    %d [label = \"%d|%c\"", i, i, nodes_[i]->letter());
        if (in_consensus[i] != -1) {
            fprintf(pFile, ", style=filled, fillcolor=goldenrod1");
        }
        fprintf(pFile,"]\n");

        for (const auto& edge: nodes_[i]->out_edges()) {
            fprintf(pFile,"    %d -> %d [label = \"%.3f\"", i, edge->end_node_id(),
                edge->total_weight());
            if (in_consensus[i] + 1 == in_consensus[edge->end_node_id()]) {
                fprintf(pFile,", color=goldenrod1");
            }
            fprintf(pFile,"]\n");
        }
        for (const auto& aid: nodes_[i]->aligned_nodes_ids()) {
            if (aid > i) {
                fprintf(pFile,"    %d -> %d [style = dotted, arrowhead = none]\n", i, aid);
            }
        }
    }
    fprintf(pFile,"}\n");
}

}
