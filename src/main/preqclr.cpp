//---------------------------------------------------------
// Copyright 2018 Ontario Institute for Cancer Research
// Written by Joanna Pineda (joanna.pineda@oicr.on.ca)
// Hamza Khan (Hamza.Khan@oicr.on.ca)
//---------------------------------------------------------
//
// preqclr.cpp -- main program
// calculates basic QC statistics



#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include "Sequence.hpp"
#include "overlap.hpp"
#include "polisher.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "edlib.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <functional>
#include <iostream>
#include <list>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <limits.h>
#include <math.h>
#include "preqclr.hpp"

#include <zlib.h>
#include <stdio.h>
#include <getopt.h>

#include "kseq.h"
#include "readpaf/paf.h"
#include "readpaf/sdict.h"

#include "zstr.hpp"
#include "strict_fstream.hpp"

#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/filewritestream.h"
#include <typeinfo>

KSEQ_INIT(gzFile, gzread)

#define VERSION "2.0"
#define SUBPROGRAM "calculate"

using namespace std;
using namespace rapidjson;

typedef PrettyWriter<StringBuffer> JSONWriter;
typedef std::chrono::duration<float> fsec;

namespace htzy
{
    //Vector to store n random reads 
    unordered_set<string> rreads;
    //Map structure to store each PAF record
    unordered_map<string, vector<sequence>> full_paf_records;     
    //Map structure to store fastq records
    map<string, string> parsed_fq;
    //Map structure to store reads and target seqs to construct SPOA
    map<string,vector<string>> seqs_for_spoa; 
    //Map structure to store SPOA MSA for a read and its targets
    map<string,vector<string>> read_msa;
    //Map structure to store allele ratios
    map<float,int> allele_ratio;
    std::vector<std::unique_ptr<racon::Overlap>> dst;
    std::vector<std::unique_ptr<racon::Sequence>> fqdst;
}

namespace opt
{
    unsigned nThrd=1;
    unsigned num_random_reads=1000;
    static unsigned int verbose;
    static unsigned int htzgy;
    static string reads_file;
    static string paf_file;
    static string gfa_file = "";
    static string sample_name;
    static int rlen_cutoff = 0;
    static bool remove_dups = false;
    static bool filter_high_cov = true;
    static bool filter_low_cov = true;
    uint32_t window_length = 500;
    double quality_threshold = 10.0;
    double error_threshold = 0.3;
    int8_t match = 5;
    int8_t mismatch = -4;
    int8_t gap = -8;
    int8_t allowed_spoa_gaps_percent = 20;
    int8_t min_spoa_coverage = 20;
}


bool endFile = false;

inline void rightStrip(const char* src, uint32_t& src_length) {
    while (src_length > 0 && isspace(src[src_length - 1])) {
        --src_length;
        }
    }

inline void rightStripHard(const char* src, uint32_t& src_length) {
    for (uint32_t i = 0; i < src_length; ++i) {
        if (isspace(src[i])) {
            src_length = i;
            break;
        }
    }
}



bool preqc_to_racon_parsepaf(std::vector<std::unique_ptr<racon::Overlap>>& d){

   d = std::move(htzy::dst);
   return true;

}

bool preqc_to_racon_parsefq(std::vector<std::unique_ptr<racon::Sequence>>& dq){

   dq = std::move(htzy::fqdst);
   return true;

}

void allele_ratio_to_json(JSONWriter* writer){
    writer->Key("allele_ratio");
    writer->StartObject(); 
    for( auto it = htzy::allele_ratio.begin(); it != htzy::allele_ratio.end(); it++){
        float ratio = it->first;
        int count = it->second;
        writer->Key((to_string(ratio)).c_str()); 
        writer->Int(count);  
    } 
    writer->EndObject();
}


unordered_set<string> random_reads(unordered_map<string,string> &temp_map){
    /* 
    ========================================================
    Returns opt::num_random_reads random reads from the given set of reads
    --------------------------------------------------------
    Input:      set of read IDs
    Output:     set of opt::num_random_reads read IDs
    ========================================================
    */
    
    vector <string> temp_vec;
    for(auto it = temp_map.begin(); it != temp_map.end(); ++it) {
        temp_vec.push_back((it->first));
    }
    std::srand(std::time(0));
    std::random_shuffle (temp_vec.begin(), temp_vec.end());
    std::vector<string> ran_reads(temp_vec.begin(), temp_vec.begin() + opt::num_random_reads);
 
    //Debugging with just two reads
    //ran_reads.clear();
    //ran_reads.push_back("S1HapC_528");
    //ran_reads.push_back("S1HapA_3");

    unordered_set<string> ran_reads_set;
    for (auto i: ran_reads) {
        ran_reads_set.insert(i);
        ran_reads_set.insert(temp_map[i]);
    } 

    cerr << "Size of ran_reads = " << ran_reads.size() << "\n";
    cerr << "Size of ran_reads_set = " << ran_reads_set.size() << "\n";
    //Print all reads id in set
   // for (auto i = ran_reads_set.begin(); i!=ran_reads_set.end();++i) {
   //    std::cerr << *i << "\n";
   // } 
   
    return ran_reads_set;
    
    
    /*Testing OMP
    std::vector<string> test_ran_reads = {"S1HapC_999", "S1HapA_1", "S1HapA_1002"};
    return test_ran_reads;
    */
}

/*
void allele_ratio_from_msa(vector<string> &msa, const char * depth_threshold, const char * percent_gaps){
     
    ===========================================================================
    Calculate the allele ratio per column of an MSA
    ----------------------------------------------------------------------------
    Input:      MSA (vector of strings of equal length), depth threshold, percent 
                allowed gap threshold in a column
    Output:     Prints allele ratios and writes it to allele_ratio map structure
    =============================================================================
            
    int MSA_len = (msa.front()).length();
    vector<map<char, int>> allele_count;
    //cout << "MSA_len" << MSA_len<<endl;
    //cout << "Thread ID inside allele_ratio_from_msa = " << omp_get_thread_num() << endl;    

    for (vector<string>::const_iterator v = msa.begin(); v != msa.end(); ++v){
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
                 map<char, int> t;
                 if(((*v).at(q))=='-' && q>=first_found && q<=last_found){
                     //cout << "\nTrue gap "  << (*v).at(q) << endl;
                     t.insert(pair<char,int>((*v).at(q),1));
                     allele_count.push_back(t);
                 }
                 if( ( (((*v).at(q))=='-') && (q<first_found) ) || ( (((*v).at(q))=='-') && (q>last_found) )) {
                     //cout << "False gap "  << (*v).at(q) << endl;
                     t.insert(pair<char,int>((*v).at(q),0));                   
                     allele_count.push_back(t);
                 }
                 if((((*v).at(q))!='-')) {
                     //cout << "Char = " << (*v).at(q) << endl;
                     t.insert(pair<char,int>((*v).at(q),1));
                     allele_count.push_back(t);
                 }
             } 

         }
         else {
             for (unsigned int q=first_found; q<=last_found; q++){
             map<char, int> t;          
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
    
    
    //Print the allele count vector
         
    int count = 0;
    for (vector<map<char, int>>::const_iterator r = allele_count.begin(); r != allele_count.end(); ++r){ 
         //cout << "POS = " << count << endl; count ++; 
         pair<char,int> max_allele ('X', 0), second_max_allele('Y', 0);           
         for(map<char, int>::const_iterator s = (*r).begin(); s!= (*r).end(); ++s){
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
             if ((msa.size()>=(atoi(depth_threshold))) &&  (second_max_allele.first!='X' && second_max_allele.first!='Y' 
                  && max_allele.first!='X' && max_allele.first!='Y')){
                  auto ntgar = roundf((float(max_allele.second)/float(second_max_allele.second+max_allele.second))*100)/100;
                  //cout << "No true gaps, Allele Ratio=" << ntgar << endl;
                  auto f = htzy::allele_ratio.find(ntgar);
                  if ( f == htzy::allele_ratio.end() ) {
                     htzy::allele_ratio[ntgar]=1;
    inline void rightStrip(const char* src, uint32_t& src_length) {
        while (src_length > 0 && isspace(src[src_length - 1])) {
            --src_length;
            }
        }

    inline void rightStripHard(const char* src, uint32_t& src_length) {
        for (uint32_t i = 0; i < src_length; ++i) {
            if (isspace(src[i])) {
                src_length = i;
                break;
            }
        }
    }


                   }
                  else {
                     htzy::allele_ratio[ntgar]+=1;
                  }
             
              }          

         }
         else {
              //cout << "Gaps ratio=" << roundf((float((*r).at('-'))/float(msa.size()))*100)/100  << endl;
              if((float((*r).at('-'))/(float(msa.size())))<=((atof(percent_gaps))/100.00) && (msa.size()>=(atoi(depth_threshold))) && 
                  (second_max_allele.first!='X' && second_max_allele.first!='Y' && max_allele.first!='X' && max_allele.first!='Y')){
                  auto ar = roundf((float(max_allele.second)/float(second_max_allele.second+max_allele.second))*100)/100;
                  //cout << "Allele Ratio=" << ar  << endl; 
                  auto f = htzy::allele_ratio.find(ar);
                  if ( f == htzy::allele_ratio.end() ) {
                     htzy::allele_ratio[ar]=1;
                   }
                  else {
                     htzy::allele_ratio[ar]+=1;
                  }

              } 
         } 
        // cout << endl;
    }    
}     
     
*/

void out( string o )
{
    // Let's handle the verbose option
    // SO: https://stackoverflow.com/questions/10150468/how-to-redirect-cin-and-cout-to-files
    string logfile = opt::sample_name + "_preqclr.log";
    ofstream out;
    out.open( logfile, ios::app );
    out << o << "\n";
    if ( opt::verbose == 1 ) {
        cout << o << "\n";
    }    
    if (endFile) {
        out.close();
    }
}

int main( int argc, char *argv[]) 
{
    // parse the input arguments, if successful it will save all the arguments
    // in the global struct opts
    parse_args(argc, argv);

    #ifdef _OPENMP
    omp_set_num_threads(opt::nThrd);
    #endif
    
    // clear any previous log files with same name
    ofstream ofs;
    ofs.open( opt::sample_name + "_preqclr.log", ofstream::out | ios::trunc );
    ofs.close();

    out("========================================================");
    out("RUNNING PREQC-LR CALCULATE");
    out("========================================================");
    auto tot_start = chrono::system_clock::now();
    auto tot_start_cpu = clock();
 
    // parse the input PAF file and return a map with key = read id, 
    // and value = read object with only needed overlap info
    // SO: https://stackoverflow.com/questions/11062804/measuring-the-runtime-of-a-c-code
    // SO1 to get cast to milliseconds: https://stackoverflow.com/questions/30131181/calculate-time-to-execute-a-function 
    out("[ Parse PAF file ] ");
    auto swc = chrono::system_clock::now();    
    auto scpu = clock();
    map<string, sequence> paf_records = parse_paf();
    auto ewc = chrono::system_clock::now();
    auto ecpu = clock();
    fsec elapsedwc = ewc - swc;
    double elapsedcpu = (ecpu - scpu)/(double)CLOCKS_PER_SEC;;
    out("[+] Time elapsed: " + to_string(elapsedwc.count()) + "s, CPU time: "  + to_string(elapsedcpu) + "s");

    // start json object
    StringBuffer s;
    JSONWriter writer(s);
    writer.StartObject();

    // add input arguments
    writer.String("sample_name");
    writer.String(opt::sample_name.c_str());

    // start calculations
    // SO: Calculating CPU time. (https://stackoverflow.com/questions/17432502/how-can-i-measure-cpu-time-and-wall-clock-time-on-both-linux-windows)
    out("\n[ Parse reads file ]");
    swc = chrono::system_clock::now();
    scpu = clock();
    auto fq_records = parse_fq ( opt::reads_file );
    ewc = chrono::system_clock::now();
    ecpu = clock();
    elapsedwc = ewc - swc;
    elapsedcpu = (ecpu - scpu)/(double)CLOCKS_PER_SEC;;
    out("[+] Time elapsed: " + to_string(elapsedwc.count()) + "s, CPU time: "  + to_string(elapsedcpu) +"s");

    swc = chrono::system_clock::now();
    scpu = clock();
    out("\n[ Calculating read length distribution ]");
    calculate_read_length( fq_records, &writer);
    ewc = chrono::system_clock::now();
    ecpu = clock();
    elapsedwc = ewc - swc;
    elapsedcpu = (ecpu - scpu)/(double)CLOCKS_PER_SEC;;
    out("[+] Time elapsed: " + to_string(elapsedwc.count()) + "s, CPU time: "  + to_string(elapsedcpu) + "s");

    out("\n[ Calculating est cov per read and est genome size ]");
    swc = chrono::system_clock::now();
    scpu = clock();
    int genome_size_est = calculate_est_cov_and_est_genome_size( paf_records, &writer);
    ewc = chrono::system_clock::now();
    ecpu = clock();
    elapsedwc = ewc - swc;
    elapsedcpu = (ecpu - scpu)/(double)CLOCKS_PER_SEC;
    out("[+] Time elapsed: " + to_string(elapsedwc.count()) + "s, CPU time: "  + to_string(elapsedcpu) +"s");

    out("\n[ Calculating GC-content per read ]");
    swc = chrono::system_clock::now();
    scpu = clock();
    calculate_GC_content( fq_records, &writer);
    ewc = chrono::system_clock::now();
    ecpu = clock();
    elapsedwc = ewc - swc;
    elapsedcpu = (ecpu - scpu)/(double)CLOCKS_PER_SEC;;
    out("[+] Time elapsed: " + to_string(elapsedwc.count()) + "s, CPU time: "  + to_string(elapsedcpu) +"s");

    out("\n[ Calculating total number of bases as a function of min read length ]");
    swc = chrono::system_clock::now();
    scpu = clock();
    calculate_tot_bases( paf_records, &writer);
    ewc = chrono::system_clock::now();
    ecpu = clock();
    elapsedwc = ewc - swc;
    elapsedcpu = (ecpu - scpu)/(double)CLOCKS_PER_SEC;;
    out("[+] Time elapsed: " + to_string(elapsedwc.count()) + "s, CPU time: "  + to_string(elapsedcpu) +"s");

    if ( !opt::gfa_file.empty() ) {
        out("\n[ Parse GFA file ] ");
        swc = chrono::system_clock::now();
        scpu = clock();
        auto contigs = parse_gfa();
        ewc = chrono::system_clock::now();
        ecpu = clock();
        elapsedwc = ewc - swc;
        elapsedcpu = (ecpu - scpu)/(double)CLOCKS_PER_SEC;;
        out("[+] Time elapsed: " + to_string(elapsedwc.count()) + "s, CPU time: "  + to_string(elapsedcpu) +"s");

        out("\n[ Calculating NGX ]");
        swc = chrono::system_clock::now();
        scpu = clock();
        calculate_ngx( contigs, genome_size_est, &writer );
        ewc = chrono::system_clock::now();
        ecpu = clock();
        elapsedwc = ewc - swc;
        elapsedcpu = (ecpu - scpu)/(double)CLOCKS_PER_SEC;;
        out("[+] Time elapsed: " + to_string(elapsedwc.count()) + "s, CPU time: "  + to_string(elapsedcpu) +"s");
    }

    //Calculate heterozygosity
    run_racon(opt::reads_file.c_str(), opt::paf_file.c_str(),  opt::reads_file.c_str(), opt::window_length,
                  opt::quality_threshold, opt::error_threshold, opt::match, opt::mismatch, opt::gap, opt::nThrd);


    if(opt::htzgy == 1){
        out("[ Calculating Heterozygosity ]");
        swc = chrono::system_clock::now();
        scpu = clock();
        //estimate_heterozygosity();   
        //allele_ratio_to_json(&writer);
        ewc = chrono::system_clock::now();
        ecpu = clock();
        elapsedwc = ewc - swc;
        elapsedcpu = (ecpu - scpu)/(double)CLOCKS_PER_SEC;;
        out("[+] Time elapsed: " + to_string(elapsedwc.count()) + "s, CPU time: "  + to_string(elapsedcpu) +"s");
    }

    // convert JSON document to string and print
    string filename = opt::sample_name + ".preqclr";

    // wrap it up
    out("\n[ Done ]");
    out("[+] Resulting preqclr file: " + filename );
    auto tot_end = chrono::system_clock::now();
    auto tot_end_cpu = clock();
    fsec tot_elapsed = tot_end - tot_start;
    double tot_elapsed_cpu = (tot_end_cpu - tot_start_cpu)/(double)CLOCKS_PER_SEC;

    writer.Key("tot_cpu_time");
    writer.Double(tot_elapsed_cpu);

    writer.EndObject();
    ofstream preqclrFILE;
    preqclrFILE.open( filename );
    preqclrFILE << s.GetString() << endl;
    preqclrFILE.close();

    endFile = true;
    out("[+] Total time: " + to_string(tot_elapsed.count()) + "s, CPU time: " + to_string(tot_elapsed_cpu) + "s");
}

vector<double> parse_gfa()
{
    // parse gfa to get the contig lengths in MB
    string line;
    ifstream infile(opt::gfa_file);
    if (!infile.is_open()) {
        fprintf(stderr, "ERROR: GFA failed to open. Check to see if it exists, is readable, and is non-empty.\n\n");
        exit(1);
    }
    vector <double> contig_lengths;
    while( getline(infile, line) ) {
        char spec;
        stringstream ss(line);
        ss >> spec;

        // get only lines that have summary information on the contig length
        if ( spec == 'x' ) {
            string ctgName;
            int ctgLen, nreads;
            ss >> spec >> ctgName >> ctgLen >> nreads;
            double len = double(ctgLen)/1000000;
            contig_lengths.push_back(len);
        }
    }
    return contig_lengths;
}

void parse_args ( int argc, char *argv[])
{
    // getopt
    extern char *optarg;
    extern int optind, opterr, optopt;
    const char* const short_opts = "t:g:c:hvzr:n:p:d:w:q:e:m:x:l:s:i:b:j:";
    const option long_opts[] = {
        {"verbose",         no_argument,        NULL,   'v'},
        {"threads",     required_argument,      NULL,   't' },
        {"htzgy",           no_argument,        NULL,   'z'},
        {"version",         no_argument,        NULL,   OPT_VERSION},
        {"reads",           required_argument,  NULL,   'r'},
        {"sample_name", required_argument,  NULL,   'n'},
        {"paf",         required_argument,  NULL,   'p'},
        {"gfa",         required_argument,  NULL,   'g'},
        {"window-length", required_argument, 0, 'w'},
        {"quality-threshold", required_argument, 0, 'q'},
        {"error-threshold", required_argument, 0, 'e'},
        {"match", required_argument, 0, 'm'},
        {"mismatch", required_argument, 0, 'x'},
        {"gap", required_argument, 0, 's'},
        {"num-random-reads", required_argument, 0, 'i'},
        {"spoa-min-cov", required_argument, 0, 'j'},
        {"spoa-gaps-percent", required_argument, 0, 'b'},
        {"help",            no_argument,    NULL,   'h'},
        {"min_rlen",     required_argument,    NULL,   'l'},
        {"keep_low_cov",     no_argument,    NULL,   OPT_KEEP_LOW_COV},	
        {"keep_high_cov",     no_argument,    NULL,   OPT_KEEP_HIGH_COV},
        {"remove_dups",     no_argument,    NULL,   OPT_REMOVE_DUPS},
        { NULL, 0, NULL, 0 }
    };

    static const char* PREQCLR_CALCULATE_VERSION_MESSAGE =
    "preqclr " SUBPROGRAM " version " VERSION "\n"
    "Written by Joanna Pineda and Hamza Khan.\n"
    "\n"
    "Copyright 2018 Ontario Institute for Cancer Research\n";

    static const char* PREQCLR_CALCULATE_USAGE_MESSAGE =
    "Usage: ./preqclr [OPTIONS] --sample_name ecoli --reads reads.fa --paf overlaps.paf --gfa layout.gfa \n"
    "Calculate information for preqclr report\n"
    "\n"
    "-v, --verbose                      Display verbose output\n"
    "    --version                      Display version\n"
    "-r, --reads                        Fasta, fastq, fasta.gz, or fastq.gz files containing reads\n"
    "-n, --sample_name                  Sample name; you can use the name of species for example (This will be used as an output prefix)\n"
    "-p, --paf                          Minimap2 Pairwise mapping Format (PAF) file \n"
    "                                   (This is produced using \'minimap2 -x ava-ont sample.fastq sample.fasta)\n"
    "-g, --gfa                          Miniasm Graph Fragment Assembly (GFA) file\n"
    "                                   This file is produced using \'miniasm -f reads.fasta overlaps.paf\'\n"
    "-z, --htzgy                        Calculate heterozygosity\n"
    "-t, --threads=N                    Use N parallel threads [1]\n"
    "-l, --min_rlen=INT                 Use overlaps with read lengths >= INT\n"   
    "--keep_low_cov          Keep reads with low coverage (<= Q25 - IQR*1.25) for genome size est. calculations \n"
    "--keep_high_cov         Keep reads with high coverage (>= Q75 + IQR*1.25) for genome size est. calculations \n\n"      
    "--remove_dups           Remove duplicate overlaps. Between duplicate overlaps choose longest alignment. \n"
    "-t, --threads=N         Use N parallel threads [1]\n"

                               
    "\nHeterozygosity calculation options (Should be used along with the -z or --htzgy option):\n\n"
    "-i, --num-random-reads=INT         Number of random reads to be used for estimating heterozygosity[1000]\n"
    "-w, --window-length=INT            Size of window on which POA is performed[500]\n"
    "-q, --quality-threshold=float      Threshold for average base quality of windows used in poa[10.0]\n"
    "-e, --error-threshold=float        Maximum allowed error rate used for filtering overlaps[0.3]\n"
    "-m, --match=INT                    Score for matching bases[5]\n"
    "-x, --mismatch=INT                 Score for mismatching bases[4]\n"
    "-s, --gap=INT                      Gap penalty (must be negative)[-8]\n"
    "-j, --spoa-min-cov=INT             Minimum coverage a column of an spoa to have for calculating allele ratio [20]\n"
    "-b, --spoa-gaps-percent=INT        Maximum percentage of gaps allowed in an spoa column while calculating allele ratios [20]\n"
    
    "\n";

    int rflag=0, nflag=0, pflag=0, gflag=0, verboseflag=0, versionflag=0;
    int c;
    uint32_t window_length = 500;
    double quality_threshold = 10.0;
    double error_threshold = 0.3;

    int8_t match = 5;
    int8_t mismatch = -4;
    int8_t gap = -8;

    while ( (c = getopt_long(argc, argv, short_opts, long_opts, NULL)) != -1 ) {
    // getopt will loop through arguments, returns -1 when end of options, and store current arg in optarg
    // if optarg is not null, keep as optarg else ""
    std::istringstream arg(optarg != NULL ? optarg : "");
    switch(c) {
        case 'v':
            opt::verbose = 1; // set verbose flag
            break;
        case 'z':
            opt::htzgy = 1; // set verbose flag
            break;
        case 't':
            arg >> opt::nThrd;
            break;
        case OPT_VERSION:
            cout << PREQCLR_CALCULATE_VERSION_MESSAGE << endl;
            exit(0);
        case 'r':
            if ( rflag == 1 ) {
                fprintf(stderr, "./preqclr: multiple instances of option -r,--reads. \n\n");
                fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE, argv[0]); 
                exit(1);
            }
            rflag = 1;
            arg >> opt::reads_file;
            break;
        case 'n':
            if ( nflag == 1 ) {
                fprintf(stderr, "./preqclr: multiple instances of option -n,--sample_name. \n\n");
                fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE, argv[0]);
                exit(1);
            }
            nflag = 1;
            arg >> opt::sample_name;
            break;
        case 'p':
            if ( pflag == 1 ) {
                fprintf(stderr, "./preqclr: multiple instances of option -p,--paf. \n\n");
                fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE, argv[0]);
                exit(1);
            }
            pflag = 1;
            arg >> opt::paf_file;
            break;
        case 'g':
            if ( gflag == 1 ) {
                fprintf(stderr, "./preqclr: multiple instances of option -g,--gfa. \n\n");
                fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE, argv[0]);
                exit(1);
            }
            gflag = 1;
            arg >> opt::gfa_file;
            break;
        case 'w':
            window_length = atoi(optarg);
            break;
        case 'q':
            quality_threshold = atof(optarg);
            break;
        case 'e':
            error_threshold = atof(optarg);
            break;
        case 'm':
            match = atoi(optarg);
            break;
        case 'x':
            mismatch = atoi(optarg);
            break;
        case 's':
            gap = atoi(optarg);
            break;
        case 'j':
            opt::min_spoa_coverage = atoi(optarg);
            break;
        case 'b':
            opt::allowed_spoa_gaps_percent = atoi(optarg);
            break;
        case 'i':
            opt::num_random_reads = atoi(optarg);
            break;
        case 'h':
            cout << PREQCLR_CALCULATE_USAGE_MESSAGE << endl;
            exit(0);
        case 'l':
            arg >> opt::rlen_cutoff;
            break;
        case OPT_KEEP_LOW_COV:		
            opt::filter_low_cov = false;		
            break;		
        case OPT_KEEP_HIGH_COV:		
            opt::filter_high_cov = false;		
            break;
        case OPT_REMOVE_DUPS:		
            opt::remove_dups = true;		
            break;
        case ':':
            if (optopt == 'c') {
                fprintf(stderr, "./preqclr: option `-%c' is missing a required argument\n", optopt);
            } else if(isprint(optopt)) {
                fprintf(stderr, "./preqclr: option `-%c' is missing a required argument\n",optopt);
            } else {
                fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
            }
                fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE, argv[0]);
            exit(1);
        case '?':
            // invalid option: getopt_long already printed an error message
            if(isprint(optopt)) {
                fprintf(stderr, "./preqclr: option `-%c' is invalid thus ignored\n",optopt);
            } else {
                fprintf(stderr, "Unknown option character `\\x%x'.\n", optopt);
            }
            break;
        }
    }
    if( argc < 4 ) {
        cerr << PREQCLR_CALCULATE_USAGE_MESSAGE << endl;
        exit(1);
    }

    // check mandatory variables and assign defaults
    if ( rflag == 0 ) {
        fprintf(stderr, "./preqclr: missing -r,--reads option\n\n");
        fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE, argv[0]);
        exit(1);
    }
    if ( nflag == 0 ) {
        fprintf(stderr, "./preqclr: missing -n,--sample_name option\n\n");
        fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE);
        exit(1);
    }
    if ( pflag == 0 ) {
        fprintf(stderr, "./preqclr: missing -p,--paf option\n\n");
        fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE);
        exit(1);
    }

};

map<string, sequence> parse_paf()
{
        /*
    ========================================================
    Parse PAF (in 2 passes)
    --------------------------------------------------------
    In first pass, we note which lines/overlaps in the PAF
    file we don't want to include.
    In second pass, we read all "good" lines/overlaps and
    perform necessary calculations (like cov, read length)
    Input:    PAF file
    Output:   Dictionary: (each entry is a read)
              key = read id
              value = read (cov, length)
    ========================================================
    */  

    // PASS 1: note overlaps we do not want
    string line1;
    const char *c1 = opt::paf_file.c_str();
    paf_file_t *fp1;
    paf_rec_t r1;
    sdict_t *d1;
    fp1 = paf_open(c1);
    if (!fp1) {
        fprintf(stderr, "ERROR: PAF file failed to open. Check to see if it exists, is readable, and is non-empty.\
n\n");
        exit(1);
    }
    d1 = sd_init();
    
    // we need to filter self overlaps and duplicate overlaps
    // self overlap: minimap2 reports an overlap between the same reads (i.e. when query read id == target read id in the PAF file)
    // duplicate overlap: the same two reads are reported to overlap (i.e. when the same query read id and target read id pair are 
    //                    seen in multiple lines in the PAF file)
    vector<int> badlines; // stores all the lines we do not want
    map<size_t, pair<int, int>> h; // stores all the hashed query read name + target read name pairs with line number and alignment length
    int ln = 0; // current line number
    unordered_map<string,string> check_read_headers;
    while (paf_read(fp1, &r1) >= 0) {
        string qname = r1.qn;
        string tname = r1.tn;          
        // check for self overlap
        if ( qname.compare(tname) == 0 ) { 
            // self-overlap: the same read
            badlines.push_back(ln);
        } 

        if ( opt::remove_dups ) {
            // create a hashkey with lexicographically smallest combination of read names
            size_t hashkey = min(hash<string>{}(qname + tname), hash<string>{}(tname + qname));
            // check if we've seen this overlap between these two reads before
            auto it = h.find(hashkey);
            if (it != h.end()) {
                // YES, duplicate detected
                // let's compare the length of overlaps. We want the longer overlap.
                int aln_len = int(r1.bl);
                if ( it->second.first < aln_len ) {
                    // if the prev. overlap between these 2 reads is shorter, we use the current line instead
                    it->second.first = aln_len;
                    // prev. overlap's line number is recorded as "bad". it will be skipped in second pass.
                    badlines.push_back(it->second.second);
                    it->second.second = ln;
                } else {
                    badlines.push_back(ln);
                }
            } else {
                // First time we've seen this pair
                int aln_len = int(r1.bl);
                h.insert(make_pair(hashkey, make_pair(aln_len, ln)));
            }
        }
        ln+=1; // read next line 
        if (check_read_headers.find(qname) == check_read_headers.end() && check_read_headers.size()<=(3*opt::num_random_reads)){
            check_read_headers[qname]=tname;}           
    }
    //for (auto v : read_headers_set)
    //    std::cout << v << "\n";    
   
   
    htzy::rreads = random_reads(check_read_headers);   
    std::cerr << "Random reads size = "<< htzy::rreads.size() << std::endl;
    
    //for (auto v : htzy::rreads)
    //    std::cerr << v << "\n"; 
    //std::cerr << "END random reads" << endl; 

    h.clear(); // free up memory
    check_read_headers.clear();
    
    sort(badlines.begin(), badlines.end());
    
    //for (auto v: badlines)
      // std::cout << v << "\n"; 

    // PASS 2: read each line in PAF file that has NOT been noted in PASS 1 as a "bad" line
    //std::vector<std::unique_ptr<racon::Overlap>> dst;
    string line;
    const char *c = opt::paf_file.c_str();
    paf_file_t *fp;
    paf_rec_t r;
    sdict_t *d;
    fp = paf_open(c);
    if (!fp) {
        fprintf(stderr, "ERROR: PAF file failed to open. Check to see if it exists, is readable, and is non-empty.\n\n");
        exit(1);
    }
    d = sd_init();
    map<string, sequence> paf_records;
    int ln1 = 0;
    int iv = 0; // index in vector		
    // the vector is sorted numerically		
    // we can loop through the vector once by storing which is the next line to avoid		
    // once we have reached this line, we can move on to the next bad line and
    
    // look out for that one while going through the next lines		
    int bad = int(badlines.at(iv)); // first bad line to watch out for		
    // read good lines in PAF  

    int count = 0;
    while (paf_read(fp, &r) >= 0) { 
       // if (ln1<bad) {
            // read each line/overlap and save each column into variable
            string qname = r.qn;
            string tname = r.tn;
            unsigned int qlen = r.ql;
            unsigned int qstart = r.qs;
            unsigned int qend = r.qe;
            unsigned int strand = r.rev;
            unsigned int tlen = r.tl;
            unsigned int tstart = r.ts;
            unsigned int tend = r.te;
            


            const char* q_name = r.qn;
            uint32_t q_name_length = strlen(q_name);
            uint32_t q_length = r.ql;
            uint32_t q_begin = r.qs; 
            uint32_t q_end = r.qe;  
            char orientation = r.rev==0 ? '+' : '-';
            const char* t_name = r.tn;
            uint32_t t_name_length = strlen(t_name);
            uint32_t t_length = r.tl; 
            uint32_t t_begin = r.ts;
            uint32_t t_end = r.te; 
            uint32_t matching_bases = r.ml; 
            uint32_t overlap_length = r.bl;
            uint32_t mapping_quality = r.mq;
            
            constexpr uint32_t kSSS = 1024;
            q_name_length = std::min(q_name_length, kSSS);
            t_name_length = std::min(t_name_length, kSSS);
            rightStripHard(q_name, q_name_length);
            rightStripHard(t_name, t_name_length);
	    //cout << "q_name = " << q_name << ", q_name_length = " << q_name_length<<", q_length = " << q_length << ",q_begin = " << q_begin <<",q_end=" <<q_end <<",orientation=" << orientation<<",t_name=" << t_name<<",t_name_length=" <<t_name_length <<",t_length=" <<t_length <<",t_begin=" <<t_begin <<",t_end=" <<t_end <<",matching_bases=" <<matching_bases <<",overlap_length=" <<overlap_length <<",mapping_quality="<<mapping_quality << "\n";
            //cerr << qname<< "," << tname << "," << qlen << ","  << qstart << "," << qend << "\n";
             // filter reads by read length
            if (( qlen >= opt::rlen_cutoff ) && ( tlen >= opt::rlen_cutoff )) {
                unsigned int qprefix_len = qstart;
                unsigned int qsuffix_len = qlen - qend - 1;
                unsigned int tprefix_len = tstart;
                unsigned int tsuffix_len = tlen - tend - 1;

      
                // calculate overlap length, we need to take into account minimap2's softclipping
                int left_clip = 0, right_clip = 0;
                if ( ( qstart != 0 ) && ( tstart != 0 ) ){
                    if ( strand == 0 ) {
                        left_clip += min(qprefix_len, tprefix_len);
                    } else {
                        left_clip += min(qprefix_len, tsuffix_len);
                    }
                }
                if ( ( qend != 0 ) && ( tend != 0 ) ){
                    if ( strand == 0 ) {
                        right_clip += min(qsuffix_len, tsuffix_len);
                    } else {
                        right_clip += min(qsuffix_len, tprefix_len);
                    }  
                }
                            
                              

            //Check if the query or target read IDs are present in rreads set,
            // if yes, then add the target and query read IDs to the rreads set
            // and store their records in full_paf_records 
            unordered_set<string>::const_iterator qfinder = htzy::rreads.find(qname);
            unordered_set<string>::const_iterator tfinder = htzy::rreads.find(tname);
            if (qfinder != htzy::rreads.end() || tfinder != htzy::rreads.end()){ 
                //Declare two sequence objects, one for query and one for target record    
                count++; 
                sequence qrec;                                                               
                sequence trec;                                                               
                qrec.set_paf(qname,tname,qlen,qstart,qend,strand,tlen,tstart,tend);          
                trec.set_paf(tname,qname,tlen,tstart,tend,strand,qlen,qstart,qend);          
                //cout << "CHECK qrec.qname = " << qrec.qname << endl;                             
                //cout << "prec.strand = " << prec.strand << endl;                           
                //cout << "prec.tname = " << prec.tname << endl;                             
                //cout << "prec.qlen = " << prec.qlen << endl;                               
                //cout << "prec.qstart = " << prec.qstart << endl;                           
                //cout << "prec.qend = " << prec.qend << endl;                 
                htzy::full_paf_records[qname].push_back(qrec);
                htzy::full_paf_records[tname].push_back(trec);
                //std::cout << "Type of orientation ="<< typeid(orientation).name() << '\n'; 
                //cout << "q_name = " << q_name << ", q_name_length = " << q_name_length<<", q_length = " << q_length << ",q_begin = " << q_begin <<",q_end=" <<q_end <<",orientation="<< orientation<<",t_name=" << t_name<<",t_name_length=" <<t_name_length <<",t_length=" <<t_length <<",t_begin=" <<t_begin <<",t_end=" <<t_end <<",matching_bases=" <<matching_bases <<",overlap_length=" <<overlap_length <<",mapping_quality="<<mapping_quality << "\n";
 
                htzy::dst.emplace_back(new racon::Overlap(q_name,  q_name_length,  q_length, q_begin,  
                q_end, orientation, t_name,
                t_name_length,  t_length,  t_begin,
                t_end,  matching_bases,  overlap_length, mapping_quality));
       
                
                htzy::dst.emplace_back(new racon::Overlap(t_name,  t_name_length,  t_length, t_begin,
                t_end, orientation, q_name,
                q_name_length,  q_length,  q_begin,
                q_end,  matching_bases,  overlap_length, mapping_quality));
                

            }
            
            // calculate coverage per read               
            // add this information to paf_records dictionary
            auto i = paf_records.find(qname);
            unsigned int overlap_len = abs(qend - qstart) + left_clip + right_clip;
            double cov = double(overlap_len) / double(qlen);
            if ( i == paf_records.end() ) {
                // if read not found initialize in paf_records
                sequence qr;
                qr.set(qlen, cov);
                paf_records.insert(pair<string,sequence>(qname, qr));
            } else {
                // if read found, update the overlap info
                i->second.updateCov(cov);
            }

            auto j = paf_records.find(tname);
            overlap_len = abs(tend - tstart) + left_clip + right_clip;
            cov = double(overlap_len) / double(tlen); 
            if ( j == paf_records.end() ) {
                // if target read not found initialize in paf_records
                sequence tr;
                tr.set(tlen, cov);
                paf_records.insert(pair<string,sequence>(tname, tr));
            } else {
                // if target read found, update the overlap info
                j->second.updateCov(cov);
            }
           }
   
   // } else if ( iv+1 < badlines.size() ) {
   //     iv+=1;
   //     bad = int(badlines.at(iv));

   // }
    ln1+=1;
    } 
    std::cerr << "Size of full_paf_records = " << htzy::full_paf_records.size() << ", Count variable = " << count << ", Size of htzy::dst" << (htzy::dst).size() << "\n";

 

 /*   //Print full_paf_records
   for(auto it = htzy::full_paf_records.begin(); it != htzy::full_paf_records.end(); ++it) {
       std::cout << "Key = " << (*it).first ;
       for(auto it2 = 0; it2 < (*it).second.size(); ++it2)
           cout << ",Value No." << it2 << " = " <<(*it).second[it2].tname << " ";
       cout <<"\n";
   }
 */  


  

   //htzy::rreads = random_reads(1000, htzy::full_paf_records);

   // XXXXXXXXXXXXXXXXXXX
   // DEBUGGING ZONE
   // XXXXXXXXXXXXXXXXXXX
   //cout << "TOTAL NUMBER OF READS: "<< paf_records.size() << endl;
   //for ( auto const& r : paf_records ) {
   //    sequence temp = r.second;
   //    cout << r.first << "\t" << temp.cov << "\t" << temp.read_len << endl;
   //} 
   // XXXXXXXXXXXXXXXXXXX

   return paf_records;
}

void calculate_ngx( vector<double> contig_lengths, double genome_size_est, JSONWriter* writer ){
    /*
    ========================================================
    Calculating NGX
    --------------------------------------------------------
    Uses GFA information to evaluate the assembly quality
    Input:      All the contig lengths in MB ****
    Output:     NGX values in a dictionary:
                key   = X
                value = contig length where summing contigs 
                with length greater than or equal to this 
                length is Xth percentile
                of the genome size estimate....
    ========================================================
    */
    // we will need to add contig lengths downstream
    // so first we need to check if addition of contig lens would cause overflow
    genome_size_est = double(genome_size_est)/1000000;
    int x = 0;
    int nx;
    // this is going to hold key = x percent of genome size estimate, value = x
    map<double, int> gx;
    while ( x <= 100 ) {
        gx.insert( make_pair((double(x) * genome_size_est)/100, x) );
        x += 1;
    }
    
    // sort in descending order the contig_lengths
    sort(contig_lengths.rbegin(), contig_lengths.rend());
   
    // this is going to hold key = x, value = ngx
    map<int, double> ngx;
    double start = 0, end = 0;
    for ( auto const& c : contig_lengths ) {
       end += double(c);
       // for all values that are less then the curr sum
       for ( auto& p : gx ) {
           if ( ( p.first >= start ) && ( p.first <= end ) ) {
               ngx.insert( make_pair(p.second, c) );
           }           
       }
       start += c;
    }

    writer->Key("ngx_values");
    writer->StartObject();
    for ( auto& p : ngx ) {
        string key = to_string(p.first);
        // x value:
        writer->Key(key.c_str());
        // ngx value:
        writer->Double(p.second);
    }
    writer->EndObject();   
}

void calculate_tot_bases( map<string, sequence> paf, JSONWriter* writer)
{
    /*
    ========================================================
    Calculating total number of bases as a function of 
    min read length
    --------------------------------------------------------
    Shows the total number of bases with varying minimum 
    read length cut offs.
    Input:      Dictionary of reads with read length info in value
    Output:     Dictionary:
                key   = read length cut off
                value = total number of bases
    ========================================================
    */

    // bin the reads by read length in BASES, and sort in decreasing order
    map < unsigned int, int, greater<unsigned int>> read_lengths;
    for( auto it = paf.begin(); it != paf.end(); it++) {
        string id = it->first;
        unsigned int r_len = it->second.read_len;

        // add as new read length if read length not in map yet
        auto j = read_lengths.find(r_len);
        if ( j == read_lengths.end() ) {
            // if read length not found
            read_lengths.insert( pair<unsigned int,int>(r_len, 1) );
        } else {
            // if read length found
            j->second += 1;
        }
    }

    writer->Key("total_num_bases_vs_min_read_length");
    writer->StartObject();
    double curr_longest;   // current longest readlength in BASES
    unsigned int nr;       // number of reads with read length
    double nb;             // total number of KILOBASES of reads with current longest read length
    double tot_num_bases = 0;
    for (const auto& p : read_lengths) {
        curr_longest = double(p.first);
        nr = p.second;
        nb = (curr_longest/double(1000.0)) * nr;

        // detect for potential overflow issues:
        // SO: https://stackoverflow.com/questions/199333/how-to-detect-integer-overflow
        // curr_longest * nr may have encountered an overflow issue
        // leading to a negative number. We do not include negative nb. 
            if (!(nb > 0) || !(tot_num_bases > INT_MAX - nb)) {
                // would not overflow 
                tot_num_bases += nb;
                string key = to_string( curr_longest );
                writer->Key(key.c_str());
                writer->Int(tot_num_bases);
            }
    }
    writer->EndObject();
}



vector <pair< double, int >> parse_fq( string file )
{
    gzFile fp;
    kseq_t *seq;
    const char *c = file.c_str();
    fp = gzopen(c, "r");
    if (fp == 0) {
        fprintf(stderr, "ERROR: reads file failed to open. Check to see if it exists, is readable, and is non-empty.\n\n");
        exit(1);
    }
    seq = kseq_init(fp);
    vector <pair< double, int >> fq_records;
    while (kseq_read(seq) >= 0) {
         string id = seq->name.s;
         string sequence = seq->seq.s;
         uint32_t name_length = strlen(id.c_str());
         uint32_t sequence_length = strlen(sequence.c_str());           


         //Added checks for fastq and fasta based on the value stored in qual.s
         if (seq->qual.s != nullptr && htzy::full_paf_records.find(id)!=htzy::full_paf_records.end()){
             htzy::parsed_fq[id] = sequence;
             //cerr << "id = " << id << ", (const char*)id.c_str() " << (const char*)id.c_str() << endl;
             const char* quality = seq->qual.s;
             uint32_t quality_length = strlen(quality);

             rightStripHard(id.c_str(), name_length);
             rightStrip(sequence.c_str(), sequence_length);
             rightStrip(quality, quality_length);
              
             //cerr << "id=" << id << ",name_length = "<<name_length << ",Seq = " <<sequence.c_str() << ",sequence_length = " << sequence_length << endl;
             htzy::fqdst.emplace_back(std::unique_ptr<racon::Sequence>(new racon::Sequence(
                    id.c_str(), name_length ,
                    (const char*)sequence.c_str(), sequence_length,
                    (const char*)quality, quality_length)));
           
         }

         if (seq->qual.s == nullptr && htzy::full_paf_records.find(id)!=htzy::full_paf_records.end()){
             htzy::parsed_fq[id] = sequence;   
             //cerr << "FASTA id = " << id << ", (const char*)id.c_str() " << (const char*)id.c_str() << endl;
             //cerr << "id=" << id << ",name_length = "<<name_length << ",Seq = " <<sequence.c_str() << ",sequence_length = " << sequence_length << endl; 
             rightStripHard(id.c_str(), name_length);
             rightStrip(sequence.c_str(), sequence_length);
             htzy::fqdst.emplace_back(std::unique_ptr<racon::Sequence>(new racon::Sequence(
                    id.c_str(), name_length ,
                    (const char*)sequence.c_str(), sequence_length
                    )));
         }

         int r_len = sequence.length();
         int gc = 0;
         // only read 40% of sequences
         if (((rand() % 10) + 1) < 4) {
             for ( int i=0; i<r_len; i++) {
                 gc += sequence[i] == 'G' || sequence[i] == 'C' ? 1 : 0;
             }
             double gc_cont = (double(gc) / double(r_len)) *100.0;
             fq_records.push_back(make_pair(gc_cont, r_len));
         } else {
             fq_records.push_back(make_pair(0, r_len));
         }
    }
    kseq_destroy(seq);
    gzclose(fp);        
    cerr << "htzy::parsed_fq.size() = " << htzy::parsed_fq.size() << endl; 
    //for(auto it = htzy::parsed_fq.begin(); it != htzy::parsed_fq.end(); ++it) {
    //   std::cout << (*it).first << "\n";
    //}    
    return fq_records;
}

void calculate_GC_content( vector <pair< double, int >> fq, JSONWriter* writer )
{
    /*
    ========================================================
    Calculating GC-content per read
    --------------------------------------------------------
    Parses through read file and counts the Cs and Gs,
    then divides by the read length.
    Input:     Path to reads FILE
    Output:    List of GC content of all reads
    ========================================================
    */
    writer->Key("read_counts_per_GC_content");
    writer->StartArray();
    map<double, int> freq; // store counts for each GC content 
    int max = 0;
    double mode = 0.0;
    for (auto &r : fq) {
         if ( r.first != 0 ) {
             writer->Double(r.first);
             auto i = freq.find(round(r.first * 10.0)/10.0); // round to nearest 100th decimal place
             if ( i == freq.end() ){ 
                 freq.insert(make_pair(round(r.first*10.0)/10.0, 1));
             } else {
                 i->second+=1;
                 if ( i->second > max ){
                     max = i->second;
                     mode = i->first;
                 }
             }
         }
    }
    writer->EndArray();
    writer->Key("peak_GC_content");
    writer->Double(mode);
}


void calculate_total_num_bases_vs_min_cov( map<double, long long int, greater<double>> cov_info, JSONWriter* writer ) {
    /*
    ========================================================
    Calculate total number of bases per min cov
    --------------------------------------------------------
    Calculate the total number of bases at varying min
    cov cut-offs
    Input:     sorted in desc by cov, tot bases per cov dict
    Output:    tot bases per min cov
    ========================================================
    */
    writer->Key("total_num_bases_vs_min_cov");
    writer->StartObject();
    long long int tot_bases = 0;
    for (auto &c : cov_info) {
         tot_bases += c.second;
         string key = to_string( c.first );
         writer->Key(key.c_str());
         writer->Int(tot_bases);
    }
    writer->EndObject();
}

void calculate_median_cov_vs_min_read_length( vector <pair<double, int>> covs, JSONWriter* writer ) {
    /*
    ========================================================
    Calculate total number of bases per min cov
    --------------------------------------------------------
    Calculate the total number of bases at varying min
    cov cut-offs
    Input:     sorted in desc by cov, tot bases per cov dict
    Output:    tot bases per min cov
    ========================================================
    */
    // sort by descending read length
    sort(covs.begin(), covs.end(), [](const pair<double,int> &left, const pair<double,int> &right) { return left.second > right.second; });

    writer->Key("median_cov_vs_min_read_length");
    writer->StartObject();
    double median_cov = 0;
    int i = 0;   // index of current read length
    int i50 = 0; // 50% 
    for (auto &c : covs) {
         string key = to_string( c.second );
         writer->Key(key.c_str());
         writer->Double(covs[i50].first);
         i+=1;
         i50 = floor(double(i)/2);
    }
    writer->EndObject();
}


double calculate_est_cov_and_est_genome_size( map<string, sequence> paf, JSONWriter* writer )
{
    /*
    ========================================================
    Calculating est cov per read and est genome size
    --------------------------------------------------------
    For each read uses length and sum of lengths of all 
    overlaps.
    Input:    PAF records dictionary
    Output:   Dictionary: (each entry is a read)
              key = est coverage
              value = read length 
    ========================================================
    */

    vector <pair<double, int>> covs;

    // make an object that will hold pair of coverage and read length
    writer->Key("per_read_est_cov_and_read_length");
    writer->StartObject();
    long long sum_len = 0;
    long double sum_cov = 0;
    int tot_reads = 0;
    map<double,long long int, greater<double>> per_cov_total_num_bases;
    for( auto it = paf.begin(); it != paf.end(); it++)
    {
        string id = it->first;
        sequence r = it->second;
        int r_len = r.read_len;
        double r_cov = r.cov;
        string key = to_string(r_cov);
        writer->Key(key.c_str());
        writer->Int(r_len);
        covs.push_back(make_pair(r_cov,r_len));        
        sum_cov += r_cov;
        tot_reads += 1;
        sum_len += r_len;

        // save total bases for each coverage level
        auto j = per_cov_total_num_bases.find(round(r_cov));
        if ( j == per_cov_total_num_bases.end() ){
            per_cov_total_num_bases.insert(pair<double, long long int>(round(r_cov), r_len));
        } else {
            j->second += r_len;
        }
    }
    writer->EndObject();

    // write total number of bases vs min cov to json
    calculate_total_num_bases_vs_min_cov(per_cov_total_num_bases, writer);

    // calculate median coverage vs min read length
    calculate_median_cov_vs_min_read_length(covs, writer);
    // calculate IQR to use as limits in plotting script
    // sort the estimated coverages
    sort(covs.begin(), covs.end());
   
   /**Print Covs 
   cout <<"PRINTING Covs" << endl;
   for (std::vector<double>::const_iterator i = covs.begin(); i != covs.end(); ++i)
        std::cout << *i << ' ';
   */ 

    // get the index of the 25th and 75th percentile item
    int i25 = ceil(covs.size() * 0.25);
    int i75 = ceil(covs.size() * 0.75); 
    double IQR = covs[i75].first - covs[i25].first;
    double bd = IQR*1.5;
    double upperbound = round(double(covs[i75].first) + bd);
    double lowerbound = round(double(covs[i25].first) - bd);
    if ( !opt::filter_low_cov ) {
        lowerbound = -1;
    }
    if ( !opt::filter_high_cov ) {
        upperbound = 100000000000;
    }


    // filter outliers by cov: include reads with coverage [Q25-IQR*1.5, Q75+IQR*1.5]
    long long sum_len_f = 0;
    long double sum_cov_f = 0;
    int tot_reads_f = 0;
    vector<double> covs_f;
    // the following are used to get the mode of distribution
    // we bin the reads by coverage incrementing by 0.25x
    // we start binning from the lowest value of cov calculated
    double l = covs[0].first;
    double u = covs[0].first + 0.25;
    double curr_largest = -1000.0;
    double mode_cov = 0;
    int count = 0;
    int i = 0;
    while ( i < covs.size() ){ // iterate through reads
        // look at reads that fall within current bin
        while ( covs[i].first >= l && covs[i].first < u ) { 
            // filter outliers: [Q25-IQR*1.5, Q75+IQR*1.5]
            if ( (covs[i].first >= lowerbound) && (covs[i].first <= upperbound) ) {
                // count how many reads have coverage in current bin
                count += 1;
                tot_reads_f += 1;
                sum_len_f += covs[i].second;
                sum_cov_f += covs[i].first;
                covs_f.push_back(covs[i].first);
            }
            i += 1;
        }
    //    cout << u << ": " << count << "\n";
        // if this bin has the most amount of reads, the coverage is the mode
        if ( count > curr_largest ) {
            curr_largest = count;
            mode_cov = u;
        }
        // next bin 
        u += 0.25;
        l += 0.25;
        count = 0;
    }
 
    // get the mean read length
    double mean_read_len = sum_len_f/double(tot_reads_f);

    // get the median coverage
    int i50 = ceil(covs_f.size() * 0.50);
    double median_cov = double(covs_f[i50]);

    // calculate estimated genome size
    double est_genome_size = ( tot_reads_f * mean_read_len ) / double(mode_cov);
    double est_genome_size1 = ( tot_reads_f * mean_read_len ) / double(median_cov);
    out("mode_cov: " + to_string(mode_cov));
    out("median_cov: " + to_string(median_cov));
    out("mean_read_length: " + to_string(mean_read_len));
    out("est_genome_size_with_mode_cov: " + to_string(est_genome_size));
    out("est_genome_size_with_median_cov: " + to_string(est_genome_size1));
    out("tot_reads: " + to_string(tot_reads_f) );
    cout <<  opt::sample_name << "," << mode_cov << "," <<  median_cov << "," << mean_read_len << "," <<  tot_reads_f << ","<< est_genome_size << "," <<  est_genome_size1 << "\n";

    
    // now store in JSON object
    writer->Key("est_cov_post_filter_info");
    writer->StartArray();
    writer->Double(lowerbound);
    writer->Double(upperbound);
    writer->Int(tot_reads);
    writer->Double(IQR);
    writer->EndArray();

    writer->Key("est_genome_size");
    writer->Double(est_genome_size);    

    writer->Key("mean_read_len");
    writer->Double(mean_read_len);

    writer->Key("median_cov");
    writer->Double(median_cov);

    writer->Key("mode_cov");
    writer->Double(mode_cov);

    writer->Key("peak_cov");
    writer->Double(round(mode_cov*100.0)/100.0);

    writer->Key("tot_reads");
    writer->Int(tot_reads_f);

    return est_genome_size;
}



void calculate_read_length( vector<pair<double, int>> fq, JSONWriter* writer)
{
    /*
    ========================================================
    Calculating read lengths
    --------------------------------------------------------
    For each read add read lengths to JSON object
    Input:        PAF records dictionary
    Output:     List of all read lengths
    ========================================================
    */

    writer->Key("read_lengths");
    writer->StartArray();

    // loop through the map of reads and get read lengths
    for(auto &r : fq) {
        writer->Int(r.second);
    }
    writer->EndArray();
}

/*
void estimate_heterozygosity(){
        
    ========================================================
    Estimating heterozygosity
    --------------------------------------------------------
    For each read in the PAF file generate SPOA MSA and calculate 
    allele ratio while considering cutoffs
    Input:      - 
    Output:     File with allele ratio for each column
    ========================================================
    
   
    //cout << "Random reads total= " <<htzy::rreads.size() << endl;
    #pragma omp parallel for schedule(guided) shared (htzy::full_paf_records, htzy::parsed_fq, htzy::read_msa, htzy::allele_ratio)       
    for(auto it = 0; it < htzy::rreads.size(); ++it){
        //cout <<"\nRead = "<< htzy::rreads[it] << " \n";
        //cout << "Thread ID  = " << omp_get_thread_num() << endl;
        vector <int> qstarts, qends;
        vector <string> subreads;
        string qn;
        for (auto it2 = 0; it2 < htzy::full_paf_records[htzy::rreads[it]].size(); ++it2){
            qn = htzy::full_paf_records[htzy::rreads[it]][it2].qname;
            string tn = htzy::full_paf_records[htzy::rreads[it]][it2].tname;          
            unsigned int qs = htzy::full_paf_records[htzy::rreads[it]][it2].qstart;
            unsigned int qe = htzy::full_paf_records[htzy::rreads[it]][it2].qend;
            unsigned int ts = htzy::full_paf_records[htzy::rreads[it]][it2].tstart;
            unsigned int te = htzy::full_paf_records[htzy::rreads[it]][it2].tend;
            unsigned int ql = htzy::full_paf_records[htzy::rreads[it]][it2].qlen;
            unsigned int tl = htzy::full_paf_records[htzy::rreads[it]][it2].tlen;
            unsigned int sd = htzy::full_paf_records[htzy::rreads[it]][it2].strand;

            //cout << qn <<" " << ql << " "<< qs <<" "<< qe <<" "<<tn <<" "<< tl<<" "<< ts <<" "<< te <<" "<<sd << endl;
            //string test = "CATAAAAGAACG";
            //string rctest = reverseComplement(test);
            //cout <<"String = " << test << " Reverse Comp = "<< rctest <<endl;  
          
            qstarts.push_back(qs);
            qends.push_back(qe);
            //cout <<"Target read is = " << htzy::parsed_fq[tn] << "\nSize of read is = " << htzy::parsed_fq[tn].size()<< endl;
            //cout << "Substr start = "<< ts << ", Substr to =" << (te-ts+1)<< endl;
            //cout << "Substr string is = " <<htzy::parsed_fq[tn].substr(ts,(te-ts+1)) << endl;
            if(sd==0)
                subreads.push_back(htzy::parsed_fq[tn].substr(ts,(te-ts+1))); 
            else
                subreads.push_back(reverseComplement(htzy::parsed_fq[tn].substr(ts,(te-ts+1))));
                
        }
                 
        auto min_qstart = *min_element(std::begin(qstarts), std::end(qstarts)); 
        auto max_qend = *max_element(std::begin(qends), std::end(qends));
        //cout <<"min_qstart = " << min_qstart << " , max_qend =" << max_qend<<endl;
        subreads.push_back((htzy::parsed_fq[qn]).substr(min_qstart,(max_qend-min_qstart)));  
        vector<string> msa = calculate_heterozygosity(subreads, const_cast<char*>("5"), const_cast<char*>("-4"),\
        const_cast<char*>("-8"), const_cast<char*>("-6"), const_cast<char*>("0"));
        htzy::read_msa[qn] = msa;
        allele_ratio_from_msa(msa, const_cast<char*>("20"), const_cast<char*>("10"));
    }
    
}



vector<std::string> calculate_heterozygosity( vector<string> &sequences, const char * a, const char * b, const char * c, const char * d, const char * e) {
    
    ========================================================
    Calculate SPOA MSA from a set of reads
    --------------------------------------------------------
    Input:      Set of reads/strings, score for matching bases
                score for mis-matching bases, gap-open penalty, gap-extend penalty,
                alignment mode: 0 - local, 1 - global, 2 - semi-global
    Output:     Prints MSA to stdout
                Returns MSA in the form of vector of strings
    ========================================================
    */

     /*std::vector<std::string> sequences = {
            "CATAAAAGAACGTAGGTCGCCCGTCCGTAACCTGTCGGATCACCGGAAAGGACCCGTAAAGTGATAATGAT",
            "ATAAAGGCAGTCGCTCTGTAAGCTGTCGATTCACCGGAAAGATGGCGTTACCACGTAAAGTGATAATGATTAT",
            "ATCAAAGAACGTGTAGCCTGTCCGTAATCTAGCGCATTTCACACGAGACCCGCGTAATGGG",
            "CGTAAATAGGTAATGATTATCATTACATATCACAACTAGGGCCGTATTAATCATGATATCATCA",
            "GTCGCTAGAGGCATCGTGAGTCGCTTCCGTACCGCAAGGATGACGAGTCACTTAAAGTGATAAT",
            "CCGTAACCTTCATCGGATCACCGGAAAGGACCCGTAAATAGACCTGATTATCATCTACAT"
     };*/

     /*
     //cout << "Thread ID inside calulate_heterozygosity = " << omp_get_thread_num() << endl;
     //auto params = SPOA::AlignmentParams(atoi(a), atoi(b), atoi(c),
       //  atoi(d), (SPOA::AlignmentType) atoi(e));

     //std::string consensus = SPOA::generate_consensus(sequences, params, true);

     //fprintf(stderr, "Consensus (%zu)\n", consensus.size());
     //fprintf(stderr, "%s\n", consensus.c_str());

     std::vector<std::string> msa;
     //SPOA::generate_msa(msa, sequences, params, true);
     
     //
    // fprintf(stderr, "Multiple sequence alignment\n");
     //for (const auto& it: msa) {
     //   fprintf(stderr, "%s\n", it.c_str());
     //}
     //
     return msa;
        
}

*/


void run_racon(const std::string& sequences_path,const std::string& overlaps_path, const std::string& target_path,
int32_t window_length, double quality_threshold, double error_threshold,int8_t match,int8_t mismatch,int8_t gap, int num_threads){

    cout << "Inside run_racon" << endl;
    bool drop_unpolished_sequences = true;
    auto polisher = racon::createPolisher(sequences_path, overlaps_path,
        target_path, racon::PolisherType::kF, window_length, quality_threshold,
        error_threshold, match, mismatch, gap, num_threads);

    polisher->initialize();

    std::vector<std::unique_ptr<racon::Sequence>> polished_sequences;
    std::vector<std::map<float,int>> ar;
    std::vector<std::map<int, std::vector<std::string>>> s_msa;
    polisher->polish(polished_sequences, ar, s_msa, drop_unpolished_sequences, opt::min_spoa_coverage, opt::allowed_spoa_gaps_percent);
    /*
    for (const auto& it: polished_sequences) {
        fprintf(stdout, ">%s\n%s\n", it->name().c_str(), it->data().c_str());
    }
    */
    std::map<float,int> final_allele_ratio;
    for(auto it = ar.begin(); it != ar.end(); ++it) {
        for( auto si = (*it).begin(); si != (*it).end(); si++) {
           auto j = (final_allele_ratio).find(si->first);
           if ( j == final_allele_ratio.end()){
               final_allele_ratio[si->first] = si->second;
           }
           else{
               final_allele_ratio[si->first]+= si->second;    
           }
         //fprintf(stdout, ">%s,%s\n", si->first().c_str(), si->second().c_str());
       }

    }

    for(auto ir = s_msa.begin(); ir != s_msa.end(); ++ir){
        for( auto is = (*ir).begin(); is != (*ir).end(); is++) {
        fprintf(stdout, "Multiple sequence alignment for %d seqs\n", is->first);
           for (const auto& it: is->second) {
               fprintf(stdout, "%s\n", it.c_str());
           }
        }
    }

    for( auto si =  final_allele_ratio.begin(); si !=  final_allele_ratio.end(); si++) {
           cout << si->first <<","<< si->second << endl;
    }


//    for (const auto& it: ar) {
//         for(const auto& si: (*it)){
//            fprintf(stdout, ">%s,%s\n", si->first().c_str(), si->second().c_str());
//         }
 //   }

}
