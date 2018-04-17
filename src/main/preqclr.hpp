//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Joanna Pineda (joanna.pineda@oicr.on.ca) and Hamza Khan (hamza.khan@oicr.on.ca)
//---------------------------------------------------------
//
// preqc-lr.hpp 
// calculates basic QC statistics

#include <map>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "sequence.hpp"

#include "readpaf/paf.h"

#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/filewritestream.h"
#include "spoa/spoa.hpp"


//using namespace htzy;
using namespace std;
using namespace rapidjson;

typedef PrettyWriter<StringBuffer> JSONWriter;

double calculate_est_cov_and_est_genome_size( map<string, sequence> paf, JSONWriter* writer );
void calculate_read_length( vector <pair <double, int>> fq, JSONWriter* writer);
void calculate_GC_content( vector <pair <double, int>> fq, JSONWriter* writer);
void calculate_tot_bases( map<string, sequence> paf, JSONWriter* writer);
void calculate_ngx( vector<double> contig_lengths, double genome_size_est, JSONWriter* writer);
void calculate_total_num_bases_vs_min_cov( map<double, long long int, greater<double>> per_cov_total_num_bases, JSONWriter* writer);

int getopt( int argc, char* const* argv[], const char *optstring);
enum { OPT_VERSION, OPT_KEEP_LOW_COV, OPT_KEEP_HIGH_COV, OPT_REMOVE_DUPS };
void calculate_total_num_bases_vs_min_cov( map<double, long long int, greater<double>> per_cov_total_num_bases, JSONWriter* writer);
void parse_args( int argc, char *argv[]);
map<string, sequence> parse_paf();
vector<double> parse_gfa();
vector<pair<double,int>> parse_fq( string readsFile );

vector<string> random_reads(int i, map<string, vector<sequence>> &temp_map);
void allele_ratio_from_msa(vector<string> &msa, const char * depth, const char * percent_gaps);

void run_racon(const std::string& sequences_path,const std::string& overlaps_path, const std::string& target_path,
int32_t window_length, double quality_threshold, double error_threshold,int8_t match,int8_t mismatch,int8_t gap, int num_threads);

void estimate_heterozygosity();
vector<string> calculate_heterozygosity(vector<string> &sequences, const char * a, const char * b, const char * c, const char * d, const char * e);
void allele_ratio_to_json(JSONWriter* writer);

