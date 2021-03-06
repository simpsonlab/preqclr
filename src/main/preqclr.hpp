//---------------------------------------------------------
// Copyright 2018 Ontario Institute for Cancer Research
// Written by Joanna Pineda (joanna.pineda@oicr.on.ca)
//---------------------------------------------------------
//
// preqclr.hpp 
// calculates basic QC statistics

#include <map>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "sequence.hpp"
#include "contig.hpp"

#include "readpaf/paf.h"

#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/filewritestream.h"

using namespace std;
using namespace rapidjson;

typedef PrettyWriter<StringBuffer> JSONWriter;

double calculate_est_cov_and_est_genome_size(map<string, sequence> paf, JSONWriter* writer);
void write_read_length(vector <pair <double, int>> fq, JSONWriter* writer);
void calculate_GC_content(vector <pair <double, int>> fq, JSONWriter* writer);
void calculate_tot_bases(map<string, sequence> paf, JSONWriter* writer);
void calculate_ngx(map<string, contig> contigs, double genome_size_est, JSONWriter* writer);
void calculate_total_num_bases_vs_min_cov(map<double, long long int, greater<double>> per_cov_total_num_bases, JSONWriter* writer);
void calculate_repetitivity(map<string, contig> ctg, double g, int n, JSONWriter* writer);
double calculateDustScore(const string& seq);
map<string, contig> calculate_ctgs();

int getopt( int argc, char* const* argv[], const char *optstring);
enum { OPT_VERSION, OPT_KEEP_LOW_COV, OPT_KEEP_HIGH_COV, OPT_KEEP_DUPS, OPT_REMOVE_INT_MATCHES, OPT_MAX_OVERHANG, OPT_MAX_OVERHANG_RATIO, OPT_REMOVE_CONTAINED, OPT_PRINT_READ_COV, OPT_KEEP_SELF_OVERLAPS, OPT_PRINT_GSE_STAT, OPT_PRINT_NEW_PAF };
void parse_args(int argc, char *argv[]);
map<string, sequence> parse_paf(JSONWriter* writer);
void parse_gfa(map<string, contig> ctgs);
vector<pair<double,int>> parse_fq(string readsFile, JSONWriter* writer);
