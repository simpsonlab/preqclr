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

double calculate_est_cov_and_est_genome_size( map<string, sequence> paf, JSONWriter* writer );
void calculate_read_length( vector <pair <double, int>> fq, JSONWriter* writer);
void calculate_GC_content( vector <pair <double, int>> fq, JSONWriter* writer);
void calculate_tot_bases( map<string, sequence> paf, JSONWriter* writer);
void calculate_ngx( vector<double> contig_lengths, double genome_size_est, JSONWriter* writer);
void calculate_total_num_bases_vs_min_cov( map<double, long long int, greater<double>> per_cov_total_num_bases, JSONWriter* writer);

int getopt( int argc, char* const* argv[], const char *optstring);
enum { OPT_VERSION, OPT_KEEP_LOW_COV, OPT_KEEP_HIGH_COV, OPT_REMOVE_DUPS };
void parse_args( int argc, char *argv[]);
map<string, sequence> parse_paf();
vector<double> parse_gfa();
vector<pair<double,int>> parse_fq( string readsFile );
