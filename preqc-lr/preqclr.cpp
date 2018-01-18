#include <vector>
#include <iostream>
#include <list>
#include <fstream>
#include <functional>
#include <algorithm>
#include <map>
#include <sstream>
#include <string>
#include "include/rapidjson/document.h"
#include "include/rapidjson/writer.h"
#include "include/rapidjson/stringbuffer.h"
#include "include/rapidjson/prettywriter.h"
#include "include/rapidjson/filereadstream.h"
#include "include/rapidjson/filewritestream.h"
#include <math.h>
#include <zlib.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include "include/readfq/kseq.h"
KSEQ_INIT(gzFile, gzread)

#define VERSION "2.0"
#define SUBPROGRAM "calculate"
#define read seq_segs

using namespace std;
using namespace rapidjson;

typedef PrettyWriter<StringBuffer> JSONWriter;

namespace opt
{
    static unsigned int verbose;
    static string reads_file;
    static string paf_file;
    static string gfa_file;
    static string type;
    static string sample_name;
}

class read{
  public:
    string read_id;
    int read_len;
    int total_len_overlaps;
    int total_num_overlaps;
  void set(string i, int l, int lol)
  {
    read_id = i;
    read_len = l;
    total_len_overlaps = lol;
    total_num_overlaps = 1; 
  }
  void updateOverlap(int new_ol_len)
  {
    total_len_overlaps+=new_ol_len;
    total_num_overlaps+=1;
  }
};

void calculate_est_cov_and_est_genome_size( map<string, read> paf, JSONWriter* writer );
void calculate_read_length( map<string, read> paf, JSONWriter* writer);
void calculate_GC_content( string readsFile, JSONWriter* writer);
void calculate_tot_bases( map<string, read> paf, JSONWriter* writer);

int getopt( int argc, char *const arv[], const char *optstring);
enum { OPT_VERSION };

int main( int argc, char *argv[]){
  ifstream inFile;
  int i;

  // getopt
  extern char *optarg;
  extern int optind, opterr, optopt;
  int rflag=0, tflag=0, nflag=0, pflag=0, gflag=0, verboseflag=0, versionflag=0;
  int c;
  const char* const short_opts = "hvr:t:n:p:g:";
  // longopts requires the long options table below...
  // the first element is a const char *name: this is the name of the option without any leading dashes
  // the sec. element specifies  whether the long option has an argument and if so what kind of argument
  // // no_argument = the option does not take an argument
  // // required_argument = the option requires an argument
  // // optional_argument = the option's argument is optional
  // the third element is a pointer that returns the value in the val field of the structure if set to NULL
  // the fourth element is the value returned if the long option is seen 
  const option long_opts[] = {
    {"verbose", 	no_argument, 		NULL,	'v'},
    {"version", 	no_argument, 		NULL, 	OPT_VERSION},
    {"reads",		required_argument, 	NULL,	'r'},
    {"type", 		required_argument, 	NULL,	't'},
    {"sample_name",	required_argument, 	NULL,	'n'},
    {"paf",			required_argument, 	NULL, 	'p'},
    {"gfa", 		optional_argument, 	NULL, 	'g'},
    {"help", 		no_argument, 		NULL, 	'h'},
    { NULL, 		0, 					NULL,	 0}
  };

  static const char* PREQCLR_CALCULATE_VERSION_MESSAGE = 
  "preqclr-" SUBPROGRAM " version " VERSION "\n"
  "Written by Joanna Pineda.\n"
  "\n"
  "Copyright 2017 Ontario Institute for Cancer Research\n"; 

  static const char* PREQCLR_CALCULATE_USAGE_MESSAGE =
  "Usage: preqclr version " VERSION " " SUBPROGRAM " [OPTIONS] --reads reads.fa --type {ont|pb} --paf overlaps.paf --gfa layout.gfa \n"
  "Calculate information for preqclr report\n"
  "\n"
  "-v, --verbose				display verbose output\n"
  "    --version				display version\n"
  "-r, --reads				Fasta, fastq, fasta.gz, or fastq.gz files containing reads\n"
  "-t, --type				Type of long read sequencer. Either pacbio (pb) or oxford nanopore technology data (ont)\n"
  "-n, --sample_name			Sample name; you can use the name of species for example. This will be used as output prefix\n"
  "-p, --paf				Minimap2 pairwise alignment file (PAF). This is produced using \'minimap2 -x ava-ont sample.fastq sample.fastq\'\n"
  "-g, --gfa				Miniasm graph gragment assembly (GFA) file. This is produced using \'miniasm -f reads.fasta overlaps.paf\'\n"
  "\n";

  while ( (c = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1 ) {
    // getopt will loop through arguments, returns -1 when end of options, and store current arg in optarg
    // if optarg is not null, keep as optarg else ""
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch(c) {
      case 'v':
        opt::verbose = 1; // set verbose flag
        break;
      case OPT_VERSION:
        cout << PREQCLR_CALCULATE_VERSION_MESSAGE << endl;
        exit(0);
      case 'r':
        rflag = 1;
        arg >> opt::reads_file;
        break;
      case 't':
        tflag = 1;
        if (( arg.str().compare("ont")  != 0 ) && ( arg.str().compare("pb") != 0 )) {
          fprintf(stderr, "preqclr %s: option -t,--type is missing a valid argument {ont,pb}. \n\n", SUBPROGRAM);
          fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE, argv[0]);
          exit(1);
        }
        arg >> opt::type;
        break;
      case 'n':
        nflag = 1;
        arg >> opt::sample_name;
        break;
      case 'p':
        pflag = 1;
        arg >> opt::paf_file;
        break;
      case 'g':
        gflag = 1;
        arg >> opt::gfa_file;
        break;
      case 'h':
        cout << PREQCLR_CALCULATE_USAGE_MESSAGE << endl;
        exit(0);
      case ':':
        fprintf(stderr, "preqclr %s: option `-%c' is missing a required argument\n", SUBPROGRAM, optopt);
        fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE, argv[0]);
        exit(1);
      case '?':
        // invalid option: getopt_long already printed an error message
        fprintf(stderr, "preqclr %s: option `-%c' is invalid: ignored\n", SUBPROGRAM, optopt);
        fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE, argv[0]);
        break;
    }
  }

  if( argc < 4 ){
    cerr << PREQCLR_CALCULATE_USAGE_MESSAGE << endl;
    return -1;
  }
  
  // print any remaining command line argumentes
  if (optind < argc) { 
    for (; optind < argc; optind++)
      cerr << "preqclr " << SUBPROGRAM << ": too many arguments: argv[optind]" << endl;
  }

  // check mandatory variables and assign defaults
  if ( rflag == 0 ) { 
    fprintf(stderr, "preqclr %s: missing -r,--reads option\n\n", SUBPROGRAM);
    fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE, argv[0]);
    exit(1);
  } 
  if ( nflag == 0 ) {
    fprintf(stderr, "preqclr %s: missing -n,--sample_name option\n\n", SUBPROGRAM);
    fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE);
    exit(1);
  }
  if ( pflag == 0 ) {
    fprintf(stderr, "preqclr %s: missing -p,--paf option\n\n", SUBPROGRAM);
    fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE);
    exit(1);
  }
  if ( tflag == 0 ) {
    fprintf(stderr, "preqclr %s: missing -t,--type option\n\n", SUBPROGRAM);
    fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE);
    exit(1);
  }

  // check if an option is used more than once
  // check if files exist and are readable ...

  // read each line in paf file passed on by user
  i=0;
  string line;
  ifstream infile(opt::paf_file);
  map<string, read> paf_records;
  float sum_left_clip = 0;
  float sum_right_clip = 0;
  float sum_ol = 0;
  while(getline(infile, line)) {
    // read each line/overlap and save each column into variable
    string qname;
    int qlen=0, qstart=0, qend=0;
    char strand;
    string tname;
    int tlen=0, tstart=0, tend=0;
    int qprefix_len=0, qsuffix_len=0, tprefix_len=0, tsuffix_len=0, overlap_len=0;

    stringstream ss(line);
    ss >> qname >> qlen >> qstart >> qend >> strand >> tname >> tlen >> tstart >> tend;
    //cout<<qname<<"\t"<<qlen<<"\t"<<qstart<<"\t"<<qend<<"\t"<<strand<<"\t"<<tname<<"\t"<<tlen<<"\t"<<tstart<<"\t"<<tend<<"\n";

    if ( qname != tname ) {
      qprefix_len = qstart;
      qsuffix_len = qlen - qend;
      tprefix_len = tstart;
      tsuffix_len = tlen - tend;
      overlap_len = qend - qstart;
      //cout << overlap_len << endl;
      // calculate overlap length, we need to take into account minimap2's softclipping
      int left_clip = 0, right_clip = 0;
      if ( ( qstart != 0 ) && ( tstart != 0 ) ){
        if ( strand == '+' ) {
          left_clip = min(qprefix_len, tprefix_len);
        } else {
          left_clip = min(qprefix_len, tsuffix_len);
        }
      }
      if ( ( qend != 0 ) && ( tend != 0 ) ){
        if ( strand == '+' ) {
          right_clip = min(qsuffix_len, tsuffix_len);
        } else {
          right_clip = min(qsuffix_len, tprefix_len);
        }
      }
      overlap_len+=left_clip;
      overlap_len+=right_clip;
      //cout << overlap_len << endl;
      sum_left_clip+=left_clip;
      sum_right_clip+=right_clip;
      // add this information to paf_records dictionary
      auto i = paf_records.find(qname);
      if ( i == paf_records.end() ) {
            // if read not found initialize in paf_records
            read qr;
            qr.set(qname, qlen, overlap_len);
            paf_records.insert(pair<string,read>(qname, qr));
        } else {
            // if read found, update the overlap info
            i->second.updateOverlap(overlap_len);
        }
        auto j = paf_records.find(tname);
        if ( j == paf_records.end() ) {
            // if target read not found initialize in paf_records
            read tr;
            tr.set(tname, tlen, overlap_len);
            paf_records.insert(pair<string,read>(tname, tr));
        } else {
            // if target read found, update the overlap info
            j->second.updateOverlap(overlap_len);
        }
  }
  }
  
  // let's check the records...
  int sum_tot_len_ol = 0; 
  auto it = paf_records.begin();
  while (it != paf_records.end())
  {
    // Accessing KEY from element pointed by it
    string read_id = it->first;
    read temp = it->second;
//    cout << read_id << ": " << temp.read_len << "," << temp.total_num_overlaps << "," << temp.total_len_overlaps<< "\n";
    sum_tot_len_ol+=temp.total_len_overlaps;  

    // Increment the Iterator to point to next entry
    it++;
  }
//  cout <<paf_records.size() <<endl;
//  printf( "\n%f", sum_left_clip );
//  printf( "\n%f", sum_right_clip );
  printf( "\nsum_tot_len_ol: %i", sum_tot_len_ol );
  // start the JSON object
  StringBuffer s;
  JSONWriter writer(s);

  writer.StartObject();

  // add input arguments
  writer.String("sample_name");
  writer.String(opt::sample_name.c_str());

  // start calculations
//  calculate_read_length( paf_records, &writer);
//  calculate_est_cov_and_est_genome_size( paf_records, &writer);
//  calculate_GC_content( opt::reads_file, &writer);
//  calculate_tot_bases( paf_records, &writer);

  // convert JSON document to string and print
  writer.EndObject();
  //cout << s.GetString() << endl;

}

void calculate_tot_bases( map<string, read> paf, JSONWriter* writer)
{
  // bin the reads by read length
  map < int, int, greater<int>> read_lengths;
  for( auto it = paf.begin(); it != paf.end(); it++)
  {
    string id = it->first;
    int r_len = it->second.read_len;

    // add as new read length if read length not in map yet
    auto j = read_lengths.find(r_len);
    if ( j == read_lengths.end() ) {
      // if read length not found
      read_lengths.insert( pair<int,int>(r_len, 1) );
    } else {
      // if read length found
      j->second+=1;
    }
  }

  // sort read lengths
  // Display items in sorted order of keys
  int curr_longest;
  int tot_num_bases = 0;
  writer->Key("total_num_bases_vs_min_read_length");
  writer->StartObject();
  for (const auto& p : read_lengths){
    curr_longest = p.first;
    tot_num_bases += ( p.second * p.first);
    string key = to_string( curr_longest );
    writer->Key(key.c_str());
    writer->Int(tot_num_bases);
    // cout << curr_longest << ',' << tot_num_bases << "\n";
  }
  writer->EndObject();
}

void calculate_GC_content( string file, JSONWriter* writer )
{
  // add GC content info
  // read fasta/fastq file
  // count the number of Cs and Gs in the sequence
  // store in vector
  vector < double > GC_content;
  gzFile fp;
  kseq_t *seq;
  const char *c =  file.c_str();
  fp = gzopen(c, "r");
  seq = kseq_init(fp);
  writer->Key("read_counts_per_GC_content");
  writer->StartArray();
  while (kseq_read(seq) >= 0) {
     string id = seq->name.s;
     string sequence = seq->seq.s;
     //cout << "name: " <<  seq->name.s << endl;
     //cout << "seq: " << seq->seq.s << endl;
     size_t C_count = count(sequence.begin(), sequence.end(), 'C');
     size_t G_count = count(sequence.begin(), sequence.end(), 'G');
     double r_len = sequence.length();
     double gc_cont = (double( C_count + G_count ) / r_len) *100.0;
     writer->Double(gc_cont);
     //array.PushBack(gc_cont, allocator);
  }
  writer->EndArray();
  kseq_destroy(seq);
  gzclose(fp); 
}

void calculate_est_cov_and_est_genome_size( map<string, read> paf, JSONWriter* writer )
{

  vector< pair < float, int > > covs;
  double r_len, r_tot_len_ol, r_tot_num_ol, r_cov;

  // make an object that will hold pair of coverage and read length
  writer->Key("per_read_est_cov_and_read_length");
  writer->StartObject();
  float sum_tot_len_ol = 0;
  float sum_tot_ol = 0;
  for( auto it = paf.begin(); it != paf.end(); ++it)
  {
    float r_len, r_tot_len_ol, r_tot_num_ol, r_cov;
    string id = it->first;
    read* r = &it->second;
    r_len = float(r->read_len);
    sum_tot_len_ol += float(r->total_len_overlaps);
    r_tot_len_ol = float(r->total_len_overlaps);
    r_tot_num_ol = float(r->total_num_overlaps);
    sum_tot_ol += r_tot_num_ol;
    r_cov = r_tot_len_ol / r_len;
    //cout<<id<<'\t'<<r_len<<'\t'<<r_tot_num_ol<<'\n';
    string key = to_string(r_cov);
    writer->Key(key.c_str());
    writer->Int(r_len);
    covs.push_back(make_pair(r_cov,r_len));
  }
  printf( "%f\n", sum_tot_ol );
  printf( "%f\n", sum_tot_len_ol);
  cout<<"sum_tot_ol: "<<sum_tot_ol<<endl;
  cout<<"sum_tot_len_ol: "<<sum_tot_len_ol<<endl;
  writer->EndObject();
  // filter the coverage: remove if outside 1.5*interquartile_range
  // calculate IQR
  // sort the estimated coverages
  sort(covs.begin(), covs.end());

  // get the index of the 25th and 75th percentile item
  int i25 = ceil(covs.size() * 0.25);
  int i75 = ceil(covs.size() * 0.75);
  //cout<<covs[i25].first<<endl;

  // find IQR
  int IQR = covs[i75].first - covs[i25].first;
  double bd = double(IQR)*1.5;
  double lowerbound = double(covs[i25].first) - bd;
  double upperbound = double(covs[i75].first) + bd;
  cout<<lowerbound<<","<<upperbound<<endl;
  // create a new list with filters
  // stores info of set of reads after filtering:
  float sum_len = 0;
  float sum_cov = 0;
  float tot_reads = 0;
  lowerbound = 39.5;
  upperbound = 99.5;
  for( auto it = covs.begin(); it != covs.end(); ++it) {
    if (( it->first > lowerbound ) && ( it->first < upperbound ))
      sum_len+=it->second;
      sum_cov+=it->first;
      tot_reads+=1;
  }

  // calculate mean read length
  float mean_read_len;
  mean_read_len = sum_len / tot_reads;

  float mean_cov;
  mean_cov = sum_cov / tot_reads;

  // calculate the estimated genome size  
  float n, l, c, est_genome_size;
  n = tot_reads;
  l = mean_read_len;
  c = mean_cov;
  est_genome_size = ( n * l ) / c;
  cout<<"sum_cov"<<sum_cov<<"sum_len"<<sum_len<<","<<tot_reads<<"n,"<<l<<"l,"<<c<<"c,"<<est_genome_size<<std::endl;

  // now store in JSON object
  writer->Key("est_cov_post_filter_info");
  writer->StartArray();
  writer->Double(lowerbound);
  writer->Double(upperbound);
  writer->Int(n);
  writer->Double(IQR);
  writer->EndArray();

  writer->Key("est_genome_size");
  writer->Double(est_genome_size); 
 
}

void calculate_read_length( map<string, read> paf, JSONWriter* writer)
{
  writer->Key("per_read_read_length");
  writer->StartArray();
  // loop through the map of reads and get read lengths
  for(auto it = paf.begin(); it != paf.end(); it++) {
    read r;
    r = it->second;
    int r_len = r.read_len;
    writer->Int(r_len);
    // cout<<r_len<<'\n';
  }
  writer->EndArray();
}
