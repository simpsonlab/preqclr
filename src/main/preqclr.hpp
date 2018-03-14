//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Joanna Pineda (joanna.pineda@oicr.on.ca)
//---------------------------------------------------------
//
// preqc-lr.hpp 
// calculates basic QC statistics

#include <map>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "contig.hpp"
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
void calculate_ngx( map<string, contig> contigs, double genome_size_est, JSONWriter* writer);
void calculate_repetitivity(map<string, contig> ctg, double g, int n, JSONWriter* writer );

int getopt( int argc, char* const* argv[], const char *optstring);
enum { OPT_VERSION, OPT_DVCUTOFF, OPT_RLENCUTOFF };
void parse_args( int argc, char *argv[]);
map<string, sequence> parse_paf();
map<string, contig> parse_gfa();
vector<pair<double,int>> parse_fq( string readsFile );
