#include <algorithm>
#include <fstream>
#include <functional>
#include <iostream>
#include <list>
#include <map>
#include <math.h>
#include <sstream>
#include <string>
#include <vector>

#include "preqclr.hpp"
#include <chrono>

#include <zlib.h>
#include <stdio.h>
#include <getopt.h>

#include "readfq/kseq.h"

#include "rapidjson/document.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include "rapidjson/prettywriter.h"
#include "rapidjson/filereadstream.h"
#include "rapidjson/filewritestream.h"

KSEQ_INIT(gzFile, gzread)

#define VERSION "2.0"
#define SUBPROGRAM "calculate"
#define read read_seq

using namespace std;
using namespace rapidjson;

typedef PrettyWriter<StringBuffer> JSONWriter;

namespace opt
{
    static unsigned int verbose;
    static string reads_file;
    static string paf_file;
    static string gfa_file = "";
    static string type;
    static string sample_name;
}

int main( int argc, char *argv[]) 
{
    // parse the input arguments, if successful it will save all the arguments
    // in the global struct opts
    parse_args(argc, argv);

    // Let's handle the verbose option
    SO: https://stackoverflow.com/questions/10150468/how-to-redirect-cin-and-cout-to-files
    if ( opt::verbose != 1 ) {
        // if verbose option not found, then redirect all cout to preqclr.log file
        ofstream out( "preqclr.log" );
        streambuf *coutbuf = cout.rdbuf();
        cout.rdbuf(out.rdbuf());
    }
    cout << "========================================================" << endl;
    cout << "RUNNING PREQC-LR CALCULATE" << endl;
    cout << "========================================================" << endl;
    auto tot_start = chrono::system_clock::now();
 
    // parse the input PAF file and return a map with key = read id, 
    // and value = read object with all needed read info
    // SO: https://stackoverflow.com/questions/11062804/measuring-the-runtime-of-a-c-code
    // SO1 to get cast to milliseconds: https://stackoverflow.com/questions/30131181/calculate-time-to-execute-a-function 
    cout << "[ Parse PAF file ] " << endl;
    auto start = chrono::system_clock::now();    
    map<string, read> paf_records = parse_paf();
    auto end = chrono::system_clock::now();
    auto elapsed = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    cout << "[+] Time elapsed: " << elapsed << " milliseconds" << endl;

    // start json object
    StringBuffer s;
    JSONWriter writer(s);
    writer.StartObject();

    // add input arguments
    writer.String("sample_name");
    writer.String(opt::sample_name.c_str());

    // start calculations
    cout << "[ Calculating read length distribution ]" << endl;
    start = chrono::system_clock::now();
    calculate_read_length( paf_records, &writer);
    end = chrono::system_clock::now();
    elapsed = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    cout << "[+] Time elapsed: " << elapsed << " milliseconds" << endl;

    cout << "[ Calculating est cov per read and est genome size ]" << endl;
    start = chrono::system_clock::now();
    float genome_size_est = calculate_est_cov_and_est_genome_size( paf_records, &writer);
    end = chrono::system_clock::now();
    elapsed = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    cout << "[+] Time elapsed: " << elapsed << " milliseconds" << endl;

    cout << "[ Calculating GC-content per read ]" << endl;
    start = chrono::system_clock::now();
    calculate_GC_content( opt::reads_file, &writer);
    end = chrono::system_clock::now();
    elapsed = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    cout << "[+] Time elapsed: " << elapsed << " milliseconds" << endl;

    cout << "[ Calculating total number of bases as a function of min read length ]" << endl;
    start = chrono::system_clock::now();
    calculate_tot_bases( paf_records, &writer);
    end = chrono::system_clock::now();
    elapsed = chrono::duration_cast<chrono::milliseconds>(end - start).count();
    cout << "[+] Time elapsed: " << elapsed << " milliseconds" << endl;

    if ( opt::gfa_file != "" ) {
        cout << "[ Parse GFA file ] " << endl;
        start = chrono::system_clock::now();
        vector<int> contigs = parse_gfa();
        end = chrono::system_clock::now();
        elapsed = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        cout << "[+] Time elapsed: " << elapsed << " milliseconds" << endl;

        cout << "[ Calculating NGX ]" << endl;
        start = chrono::system_clock::now();
        calculate_ngx( contigs, genome_size_est, &writer );
        end = chrono::system_clock::now();
        elapsed = chrono::duration_cast<chrono::milliseconds>(end - start).count();
        cout << "[+] Time elapsed: " << elapsed << " milliseconds" << endl;
    }
    
    // convert JSON document to string and print
    writer.EndObject();
    ofstream preqclrFILE;
    string filename = opt::sample_name + ".preqclr"; 
    preqclrFILE.open( filename );
    preqclrFILE << s.GetString() << endl;
    preqclrFILE.close();
    cout << "[+] Resulting preqclr file: " << filename << endl;
 
    // wrap it up
    auto tot_end = chrono::system_clock::now();
    auto tot_elapsed = chrono::duration_cast<chrono::milliseconds>(tot_end - tot_start).count();
    cout << "[ Done ]" << endl;
    cout << "[+] Total time: " << tot_elapsed << " milliseconds" << endl;
}

vector<int> parse_gfa()
{
    // parse gfa to get the contig lengths
    string line;
    ifstream infile(opt::gfa_file);
    vector <int> contig_lengths;
    while( getline(infile, line) ) {
        char spec;
        stringstream ss(line);
        ss >> spec;

        // get only lines that have information on the contig length
        if ( spec == 'S' ) {
            string contig_id, contig_seq, contig_len;
            ss >> spec >> contig_id >> contig_seq >> contig_len;
            const string toErase = "LN:i:";
            size_t pos = contig_len.find(toErase);

            // Search for the substring in string in a loop untill nothing is found
            if (pos != string::npos)
            {
                // If found then erase it from string
                contig_len.erase(pos, toErase.length());
            }
            int len = stoi(contig_len);
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
    const char* const short_opts = "hvr:t:n:p:g:";
    const option long_opts[] = {
        {"verbose",         no_argument,        NULL,   'v'},
        {"version",         no_argument,        NULL,   OPT_VERSION},
        {"reads",           required_argument,  NULL,   'r'},
        {"type",            required_argument,  NULL,   't'},
        {"sample_name", required_argument,  NULL,   'n'},
        {"paf",         required_argument,  NULL,   'p'},
        {"gfa",         optional_argument,  NULL,   'g'},
        {"help",            no_argument,    NULL,   'h'},
        { NULL,         0,  NULL,   0}
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
    "-v, --verbose              display verbose output\n"
    "        --version              display version\n"
    "-r, --reads                Fasta, fastq, fasta.gz, or fastq.gz files containing reads\n"
    "-t, --type             Type of long read sequencer. Either pacbio (pb) or oxford nanopore technology data (ont)\n"
    "-n, --sample_name          Sample name; you can use the name of species for example. This will be used as output prefix\n"
    "-p, --paf              Minimap2 pairwise alignment file (PAF). This is produced using \'minimap2 -x ava-ont sample.fastq sample.fasta"
    "\n"
    "-g, --gfa              Miniasm graph gragment assembly (GFA) file. This is produced using \'miniasm -f reads.fasta overlaps.paf\'\n"
    "\n";

    int rflag=0, tflag=0, nflag=0, pflag=0, gflag=0, verboseflag=0, versionflag=0;
    int c;
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
            if (( arg.str().compare("ont")    != 0 ) && ( arg.str().compare("pb") != 0 )) {
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
    if( argc < 4 ) {
        cerr << PREQCLR_CALCULATE_USAGE_MESSAGE << endl;
        exit(1);
    }

    // print any remaining command line argumentes
    if (optind < argc) {
        for (; optind < argc; optind++)
            cerr << "preqclr " << SUBPROGRAM << ": too many arguments:"
                 << " argv[optind]" << endl;
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

};

map<string, read> parse_paf()
{
    // read each line in paf file passed on by user
    string line;
    ifstream infile(opt::paf_file);
    map<string, read> paf_records;
    while( getline(infile, line) ) {

        // read each line/overlap and save each column into variable
        string qname;
        unsigned int qlen, qstart, qend;
        char strand;
        string tname;
        unsigned int tlen, tstart, tend;
        stringstream ss(line);
        ss >> qname >> qlen >> qstart >> qend >> strand >> tname >> tlen >> tstart >> tend;

        if ( qname.compare(tname) != 0 ) {
            unsigned int qprefix_len = qstart;
            unsigned int qsuffix_len = qlen - qend;
            unsigned int tprefix_len = tstart;
            unsigned int tsuffix_len = tlen - tend;

            // calculate overlap length, we need to take into account minimap2's softclipping
            int left_clip = 0, right_clip = 0;
            if ( ( qstart != 0 ) && ( tstart != 0 ) ){
                if ( strand == '+' ) {
                    left_clip += min(qprefix_len, tprefix_len);
                } else {
                    left_clip += min(qprefix_len, tsuffix_len);
                }
            }
            if ( ( qend != 0 ) && ( tend != 0 ) ){
                if ( strand == '+' ) {
                    right_clip += min(qsuffix_len, tsuffix_len);
                } else {
                    right_clip += min(qsuffix_len, tprefix_len);
                }
            }

            unsigned long overlap_len = (qend - qstart) + left_clip + right_clip;

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
    for ( auto const& r : paf_records ) {
        // Accessing KEY from element pointed by it
        string read_id = r.first;
        read temp = r.second;
    }

   return paf_records;
}

void calculate_ngx( vector<int> contig_lengths, int genome_size_est, JSONWriter* writer ){
    /*
    ========================================================
    Calculating NGX
    --------------------------------------------------------
    Uses GFA information to evaluate the assembly quality
    Input:      All the contig lengths
    Output:     NGX values in a dictionary:
                key   = X
                value = contig length where summing contigs 
                with length greater than or equal to this 
                length is Xth percentile
                of the genome size estimate....
    ========================================================
    */
    int x = 0;
    int nx;
    map<int,int> ngx;
    while ( x < 100 ) {
        ngx.insert( make_pair((float(x) * genome_size_est)/100, 0) );
        x += 1;
    }
    
    // sort in descending order the contig_lengths
    sort(contig_lengths.rbegin(), contig_lengths.rend());
    
    int start = 0, end = 0;
    for ( auto const& c : contig_lengths ) {
       end += c;
       // for all percentile values that are less then the curr sum
       for ( auto& p : ngx ) {
           if ( ( p.first >= float(start) ) && ( p.first <= float(end) ) ) {
               p.second = c;
           }           
       }
    }
    writer->Key("ngx_values");
    writer->StartObject();
    for ( auto& p : ngx ) {
        string key = to_string(p.first);
        writer->Key(key.c_str());
        writer->Int(p.second);
    }
    writer->EndObject();   
}

void calculate_tot_bases( map<string, read> paf, JSONWriter* writer)
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

    // bin the reads by read length
    map < unsigned long, unsigned long, greater<unsigned long>> read_lengths;
    for( auto it = paf.begin(); it != paf.end(); it++) {
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
    // display items in sorted order of keys
    writer->Key("total_num_bases_vs_min_read_length");
    writer->StartObject();
    unsigned long curr_longest;
    unsigned long tot_num_bases = 0;
    for (const auto& p : read_lengths) {
        curr_longest = p.first;
        tot_num_bases += ( p.second * p.first);
        string key = to_string( curr_longest );
        writer->Key(key.c_str());
        writer->Int(tot_num_bases);
    }
    writer->EndObject();
}

void calculate_GC_content( string file, JSONWriter* writer )
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

    vector < double > GC_content;
    gzFile fp;
    kseq_t *seq;
    const char *c = file.c_str();
    fp = gzopen(c, "r");
    seq = kseq_init(fp);
    writer->Key("read_counts_per_GC_content");
    writer->StartArray();
    while (kseq_read(seq) >= 0) {
         string id = seq->name.s;
         string sequence = seq->seq.s;
         size_t C_count = count(sequence.begin(), sequence.end(), 'C');
         size_t G_count = count(sequence.begin(), sequence.end(), 'G');
         double r_len = sequence.length();
         double gc_cont = (double( C_count + G_count ) / r_len) *100.0;
         writer->Double(gc_cont);
    }
    writer->EndArray();
    kseq_destroy(seq);
    gzclose(fp); 
}

float calculate_est_cov_and_est_genome_size( map<string, read> paf, JSONWriter* writer )
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

    vector< pair < int, int > > covs;
    double r_len, r_tot_len_ol, r_tot_num_ol, r_cov;

    // make an object that will hold pair of coverage and read length
    writer->Key("per_read_est_cov_and_read_length");
    writer->StartObject();
    for( auto it = paf.begin(); it != paf.end(); it++)
    {
        string id = it->first;
        read r = it->second;
        float r_len = float(r.read_len);
        unsigned long r_tot_len_ol = r.total_len_overlaps;
        float r_cov = r_tot_len_ol / r_len;
        string key = to_string(r_cov);
        writer->Key(key.c_str());
        writer->Int(r_len);
        covs.push_back(make_pair(r_cov,r_len));
    }
    writer->EndObject();

    // filter the coverage: remove if outside 1.5*interquartile_range
    // calculate IQR
    // sort the estimated coverages
    sort(covs.begin(), covs.end());

    // get the index of the 25th and 75th percentile item
    int i25 = ceil(covs.size() * 0.25);
    int i75 = ceil(covs.size() * 0.75);

    // find IQR
    int IQR = covs[i75].first - covs[i25].first;
    double bd = IQR*1.5;
    int lowerbound = double(covs[i25].first) - bd;
    int upperbound = double(covs[i75].first) + bd;

    // create a new set after applying this filter
    // stores info of set of reads after filtering:
    unsigned long sum_len = 0;
    unsigned long sum_cov = 0;
    float tot_reads = 0;
    for( auto it = covs.begin(); it != covs.end(); ++it) {
        if (( it->first > lowerbound ) && ( it->first < upperbound )) {
            unsigned long co = it->first;
            unsigned long le = it->second;
            sum_len += le;
            sum_cov += co;
            tot_reads += 1;
        }
    }

    // calculate estimated genome size
    float mean_read_len = sum_len / tot_reads;
    float mean_cov = sum_cov / tot_reads;
    float est_genome_size = ( tot_reads * mean_read_len ) / mean_cov;

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
    return est_genome_size;
}

void calculate_read_length( map<string, read> paf, JSONWriter* writer)
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

    writer->Key("per_read_read_length");
    writer->StartArray();

    // loop through the map of reads and get read lengths
    for(auto it = paf.begin(); it != paf.end(); it++) {
        read r = it->second;
        int r_len = r.read_len;
        writer->Int(r_len);
    }
    writer->EndArray();
}
