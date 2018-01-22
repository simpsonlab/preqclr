#include <map>
#include <string>
#include <stdio.h>
#include <stdlib.h>
#include "read_seq.hpp"

#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/filewritestream.h"

using namespace std;
using namespace rapidjson;

typedef PrettyWriter<StringBuffer> JSONWriter;

void calculate_est_cov_and_est_genome_size( map<string, read_seq> paf, JSONWriter* writer );
void calculate_read_length( map<string, read_seq> paf, JSONWriter* writer);
void calculate_GC_content( string readsFile, JSONWriter* writer);
void calculate_tot_bases( map<string, read_seq> paf, JSONWriter* writer);

int getopt( int argc, char* const* argv[], const char *optstring);
enum { OPT_VERSION };
void parse_args( int argc, char *argv[]);
map<string, read_seq> parse_paf();

