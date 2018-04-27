/*!
 * @file window.cpp
 *
 * @brief Window class source file
 */

#include <algorithm>

#include <map>
#include <stdio.h>
#include <iostream>
#include "window.hpp"
#include <vector>
#include <string>
#include <assert.h>
#include "spoa/spoa.hpp"
#include <string.h>
/*
namespace htzy
{
    //Map structure to store allele ratios
    std::map<float,int> allele_ratio;
}
*/

namespace racon {

//std::map<float,int> allele_ratio;
  
std::unique_ptr<Window> createWindow(uint64_t id, uint32_t rank, WindowType type,
    const char* backbone, uint32_t backbone_length, const char* quality,
    uint32_t quality_length) {

    if (backbone_length == 0 || backbone_length != quality_length) {
        fprintf(stderr, "[racon::createWindow] error: "
            "empty backbone sequence/unequal quality length!\n");
        exit(1);
    }

    return std::unique_ptr<Window>(new Window(id, rank, type, backbone,
        backbone_length, quality, quality_length));
}

Window::Window(uint64_t id, uint32_t rank, WindowType type, const char* backbone,
    uint32_t backbone_length, const char* quality, uint32_t quality_length)
        : id_(id), rank_(rank), type_(type), consensus_(), msa_consensus_(), allele_ratio_(), sequences_(),
        qualities_(), positions_() {

    sequences_.emplace_back(backbone, backbone_length);
    qualities_.emplace_back(quality, quality_length);
    positions_.emplace_back(0, 0);
}

Window::~Window() {
}

void Window::add_layer(const char* sequence, uint32_t sequence_length,
    const char* quality, uint32_t quality_length, uint32_t begin, uint32_t end) {

    if (quality != nullptr && sequence_length != quality_length) {
        fprintf(stderr, "[racon::Window::add_layer] error: "
            "unequal quality size!\n");
        exit(1);
    }
    if (begin >= end || begin > sequences_.front().second || end > sequences_.front().second) {
        fprintf(stderr, "[racon::Window::add_layer] error: "
            "layer begin and end positions are invalid!\n");
        exit(1);
    }

    sequences_.emplace_back(sequence, sequence_length);
    qualities_.emplace_back(quality, quality_length);
    positions_.emplace_back(begin, end);
}



std::map<float,int> Window::allele_ratio_from_msa(std::vector<std::string> &msa, const char * depth_threshold, const char * percent_gaps){
    /* 
    ===========================================================================
    Calculate the allele ratio per column of an MSA
    ----------------------------------------------------------------------------
    Input:      MSA (vector of strings of equal length), depth threshold, percent 
                allowed gap threshold in a column
    Output:     Prints allele ratios and writes it to allele_ratio map structure
    =============================================================================
    */        
    //std::cout << "Inside allele ratio from msa"  << std::endl;
    int MSA_len = (msa.front()).length();
    std::map<float,int> allele_ratio;
    std::vector<std::map<char, int>> allele_count;
    //cout << "MSA_len" << MSA_len<<endl;
    //cout << "Thread ID inside allele_ratio_from_msa = " << omp_get_thread_num() << endl;    

    for (std::vector<std::string>::const_iterator v = msa.begin(); v != msa.end(); ++v){
         //std::cout <<"SEQ = "<<*v << "\n";
         assert ((*v).length()==MSA_len);
        
         //Store first and last occurence of ATGC to know true gaps
         int first_found = (*v).find_first_of("ATGC");
         int last_found = (*v).find_last_of("ATGC");         
         //cout << "first_found = " << first_found << endl;
         //cout << "last_found = " << last_found << endl;         

         if(v==msa.begin()){
             for (unsigned int q=0; q<MSA_len; q++){
                 //cout << "q=" << q << endl;
                 std::map<char, int> t;
                 if(((*v).at(q))=='-' && q>=first_found && q<=last_found){
                     //cout << "\nTrue gap "  << (*v).at(q) << endl;
                     t.insert(std::pair<char,int>((*v).at(q),1));
                     allele_count.push_back(t);
                 }
                 if( ( (((*v).at(q))=='-') && (q<first_found) ) || ( (((*v).at(q))=='-') && (q>last_found) )) {
                     //cout << "False gap "  << (*v).at(q) << endl;
                     t.insert(std::pair<char,int>((*v).at(q),0));                   
                     allele_count.push_back(t);
                 }
                 if((((*v).at(q))!='-')) {
                     //cout << "Char = " << (*v).at(q) << endl;
                     t.insert(std::pair<char,int>((*v).at(q),1));
                     allele_count.push_back(t);
                 }
             } 

         }
         else {
             for (unsigned int q=first_found; q<=last_found; q++){

             std::map<char, int> t;          
             auto f = allele_count[q].find((*v).at(q));
                 if ( f == allele_count[q].end() ) {
                     allele_count[q][(*v).at(q)]=1;
                 }
                 else {
                     allele_count[q][(*v).at(q)] +=1;  
                 }            
             }   
         }

    }
    //cout << "\n\nallele_count.size() = " << allele_count.size() << endl;
    
    /*
    Print the allele count vector
    */     
    int count = 0;
    for (std::vector<std::map<char, int>>::const_iterator r = allele_count.begin(); r != allele_count.end(); ++r){ 
         //cout << "POS = " << count << endl; count ++; 
         std::pair<char,int> max_allele ('X', 0), second_max_allele('Y', 0);           
         for(std::map<char, int>::const_iterator s = (*r).begin(); s!= (*r).end(); ++s){
             //std::cout << s->first << " " << s->second << " ";
             if((s->second > max_allele.second) && (s->second > second_max_allele.second) && (s->first!='-')){
                 second_max_allele.first = max_allele.first; second_max_allele.second= max_allele.second;
                 max_allele.first = s->first; max_allele.second= s->second;               
                }
             if((s->second <= max_allele.second) && (s->second > second_max_allele.second) && (s->first!='-') && (s->first!=max_allele.first)){
                 second_max_allele.first = s->first; second_max_allele.second= s->second;
                }              
         }
         //cout << "\nmax_allele = " << max_allele.first <<"="<< max_allele.second << endl;
         //cout << "second_max_allele = " << second_max_allele.first <<"="<< second_max_allele.second << endl;
         
         
         auto g = (*r).find('-');
         if ( g == (*r).end()){
             if ((msa.size()>=(atoi(depth_threshold))) &&  (msa.size()<=80)  &&  (second_max_allele.first!='X' && second_max_allele.first!='Y' 
                  && max_allele.first!='X' && max_allele.first!='Y')){
                  auto ntgar = roundf((float(max_allele.second)/float(second_max_allele.second+max_allele.second))*100)/100;
                  //std::cout << "No true gaps, Allele Ratio=" << ntgar << std::endl;
                  auto f = allele_ratio.find(ntgar);
                  if ( f == allele_ratio.end() ) {
                     allele_ratio[ntgar]=1;
                   }
                  else {
                     allele_ratio[ntgar]+=1;
                  }
             
              }          

         }
         else {
              //cout << "Gaps ratio=" << roundf((float((*r).at('-'))/float(msa.size()))*100)/100  << endl;
              if((float((*r).at('-'))/(float(msa.size())))<=((atof(percent_gaps))/100.00) && (msa.size()>=(atoi(depth_threshold))) && (msa.size()<=80) && 
                  (second_max_allele.first!='X' && second_max_allele.first!='Y' && max_allele.first!='X' && max_allele.first!='Y')){
                  auto ar = roundf((float(max_allele.second)/float(second_max_allele.second+max_allele.second))*100)/100;
                  //std::cout << "Allele Ratio=" << ar  << std::endl; 
                  auto f = allele_ratio.find(ar);
                  if ( f == allele_ratio.end() ) {
                     allele_ratio[ar]=1;
                   }
                  else {
                     allele_ratio[ar]+=1;
                  }

              } 
         } 
        // cout << endl;
    }    
  return allele_ratio;
} 



bool Window::generate_consensus(std::shared_ptr<spoa::AlignmentEngine> alignment_engine, int8_t min_spoa_coverage, int8_t allowed_spoa_gaps_percent) {

    if (sequences_.size() < 3) {
        consensus_ = std::string(sequences_.front().first, sequences_.front().second);
        return false;
    }
    //std::cerr << "Inside Window::generate_consensus" << "\n";
    //std::cerr << "sequences_.size() = " << sequences_.size() << "\n";
    //std::cerr << "sequences_.front().first = " << sequences_.front().first << "\n";
    //std::cerr << "sequences_.front().second = " << sequences_.front().second << "\n";
    //std::cerr << "qualities_.front().first = " << qualities_.front().first << "\n";
    //std::cerr << "qualities_.front().second = " << qualities_.front().second << "\n";
     
   /* 
    for (uint32_t i = 0; i < sequences_.size(); ++i) {     
        std::cerr << "sequences_["<<i<<"].second = " << sequences_[i].second << "  ";
        std::cerr << "positions_["<<i<<"].first = " <<  positions_[i].first << "  ";
        std::cerr << "positions_["<<i<<"].second = " <<  positions_[i].second << "\n";    
    }*/

    auto graph = spoa::createGraph();
    graph->add_alignment(spoa::Alignment(), sequences_.front().first,
        sequences_.front().second, qualities_.front().first,
        qualities_.front().second);

    std::vector<uint32_t> rank;
    rank.reserve(sequences_.size());
    for (uint32_t i = 0; i < sequences_.size(); ++i) {
        rank.emplace_back(i);
    }

    std::sort(rank.begin() + 1, rank.end(), [&](uint32_t lhs, uint32_t rhs) {
        return positions_[lhs].first < positions_[rhs].first; });

    uint32_t offset = 0.01 * sequences_.front().second;
    for (uint32_t j = 1; j < sequences_.size(); ++j) {
        uint32_t i = rank[j];

        spoa::Alignment alignment;
        if (positions_[i].first < offset && positions_[i].second >
            sequences_.front().second - offset) {
            alignment = alignment_engine->align_sequence_with_graph(
                sequences_[i].first, sequences_[i].second, graph);
        } else {
            std::vector<int32_t> mapping;
            auto subgraph = graph->subgraph(positions_[i].first,
                positions_[i].second, mapping);
            alignment = alignment_engine->align_sequence_with_graph(
                sequences_[i].first, sequences_[i].second, subgraph);
            subgraph->update_alignment(alignment, mapping);
        }


        if (qualities_[i].first == nullptr) {
            graph->add_alignment(alignment, sequences_[i].first,
                sequences_[i].second);
        } else {
            graph->add_alignment(alignment, sequences_[i].first,
                sequences_[i].second, qualities_[i].first,
                qualities_[i].second);
        }
    }

    std::vector<uint32_t> coverages;
    consensus_ = graph->generate_consensus(coverages);
        
    std::vector<std::string> msa;
    graph->generate_multiple_sequence_alignment(msa, true);
    //msa_consensus_ = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT" + msa.at(0);   
    msa_consensus_ = msa.at(0);
    msa_consensus_.erase(std::remove(msa_consensus_.begin(), msa_consensus_.end(), '-'), msa_consensus_.end());    

    //Second graph construction. TODO: Make a function to avoid redundant code
    auto second_graph = spoa::createGraph();
    sequences_.erase(sequences_.begin());
    qualities_.erase(qualities_.begin());
    positions_.erase(positions_.begin());
    
    sequences_.insert(sequences_.begin(), std::make_pair(msa_consensus_.c_str(), msa_consensus_.length()));
    qualities_.insert(qualities_.begin(), std::make_pair(std::string((msa_consensus_.length()), '~').c_str(), (msa_consensus_.length())));
    //positions_.insert(positions_.begin(), std::make_pair(0, (msa_consensus_.length())));
    positions_.insert(positions_.begin(), std::make_pair(0, 0));


    //std::cerr << "SECOND Inside Window::generate_consensus" << "\n";
    //std::cerr << "SECOND sequences_.size() = " << sequences_.size() << "\n";
    //std::cerr << "SECOND sequences_.front().first = " << sequences_.front().first << "\n";
    //std::cerr << "SECOND sequences_.front().second = " << sequences_.front().second << "\n";
    //for (uint32_t i = 0; i < sequences_.size(); ++i) {
    //    std::cerr << "SECOND sequences_["<<i<<"].second = " << sequences_[i].second << "  ";
    //    std::cerr << "SECOND positions_["<<i<<"].first = " <<  positions_[i].first << "  ";
    //    std::cerr << "SECOND positions_["<<i<<"].second = " <<  positions_[i].second << "\n";
    //}


    second_graph->add_alignment(spoa::Alignment(), sequences_.front().first,
        sequences_.front().second, qualities_.front().first,
        qualities_.front().second);
    
    rank.clear();
    rank.reserve(sequences_.size());
    for (uint32_t i = 0; i < sequences_.size(); ++i) {
        rank.emplace_back(i);
    }
    std::sort(rank.begin() + 1, rank.end(), [&](uint32_t lhs, uint32_t rhs) {
        return positions_[lhs].first < positions_[rhs].first; });

    offset = 0.01 * sequences_.front().second;
    for (uint32_t j = 1; j < sequences_.size(); ++j) {
        uint32_t i = rank[j];

        spoa::Alignment second_alignment;
        if (positions_[i].first < offset && positions_[i].second >
            sequences_.front().second - offset) {
            second_alignment = alignment_engine->align_sequence_with_graph(
                sequences_[i].first, sequences_[i].second, second_graph);
        } else {
            std::vector<int32_t> mapping;
            auto subgraph = second_graph->subgraph(positions_[i].first,
                positions_[i].second, mapping);
            second_alignment = alignment_engine->align_sequence_with_graph(
                sequences_[i].first, sequences_[i].second, subgraph);
            subgraph->update_alignment(second_alignment, mapping);
        }

       // std::cerr << "qualities_["<<i<<"].first = " <<  qualities_[i].first << "\n";

        if (qualities_[i].first == nullptr) {
            second_graph->add_alignment(second_alignment, sequences_[i].first,
                   sequences_[i].second);
   
        } else {
            second_graph->add_alignment(second_alignment, sequences_[i].first,
                sequences_[i].second, qualities_[i].first,
                qualities_[i].second);
        }
    }
     
    std::vector<uint32_t> second_coverages;
    //consensus_ = second_graph->generate_consensus(coverages);

    std::vector<std::string> second_msa;
    second_graph->generate_multiple_sequence_alignment(second_msa);

 
    //allele_ratio_ = allele_ratio_from_msa(msa, const_cast<char*>("20"), const_cast<char*>("30")); //msa, cov, allowed percent gap
    allele_ratio_ = allele_ratio_from_msa(second_msa, const_cast<char*>(std::to_string(min_spoa_coverage).c_str()), const_cast<char*>(std::to_string(allowed_spoa_gaps_percent).c_str()));
    fprintf(stdout, "Multiple sequence alignment for %d seqs\n", sequences_.size());
    for (const auto& it: msa) {
        fprintf(stdout, "%s\n", it.c_str());
    }  
  
    fprintf(stdout, "Second Multiple sequence alignment for %d seqs\n", sequences_.size());
    for (const auto& it: second_msa) {
        fprintf(stdout, "%s\n", it.c_str());
    }

    if (type_ == WindowType::kTGS) {
        uint32_t average_coverage = (sequences_.size() - 1) / 2;

        int32_t begin = 0, end = consensus_.size() - 1;
        for (; begin < static_cast<int32_t>(consensus_.size()); ++begin) {
            if (coverages[begin] >= average_coverage) {
                break;
            }
        }
        for (; end >= 0; --end) {
            if (coverages[end] >= average_coverage) {
                break;
            }
        }

        if (begin >= end) {
            fprintf(stderr, "[racon::Window::generate_consensus] warning: "
                "contig %lu might be chimeric in window %u!\n", id_, rank_);
        } else {
            consensus_ = consensus_.substr(begin, end - begin + 1);
        }

    } else if (type_ == WindowType::kNGS) {
        uint32_t i = 0;
        for (uint32_t j = 0; i < consensus_.size(); ++i) {
            if (consensus_[i] != 'N') {
                continue;
            }

            j = std::max(j, i);
            while (j < consensus_.size() && consensus_[j] == 'N') {
                ++j;
            }

            if (j >= consensus_.size()) {
                break;
            } else if (i != j) {
                consensus_[i] = consensus_[j];
                consensus_[j] = 'N';
            }
        }
        if (i < consensus_.size()) {
            consensus_.resize(i);
        }
    }

    return true;
}

}
