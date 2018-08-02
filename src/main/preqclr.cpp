//---------------------------------------------------------
// Copyright 2018 Ontario Institute for Cancer Research
// Written by Joanna Pineda (joanna.pineda@oicr.on.ca)
//---------------------------------------------------------
//
// preqclr.cpp -- main program
// calculates basic QC statistics

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

KSEQ_INIT(gzFile, gzread)

#define VERSION "2.0"
#define SUBPROGRAM "calculate"
using namespace std;
using namespace rapidjson;

typedef PrettyWriter<StringBuffer> JSONWriter;
typedef std::chrono::duration<float> fsec;

namespace opt
{
    static unsigned int verbose;
    static string reads_file;
    static string paf_file;
    static string gfa_file = "";
    static string sample_name;
    static unsigned int rlen_cutoff = 0;
    static unsigned int olen_cutoff = 0;
    //static double olen_ratio_cutoff = 0;
    static bool keep_dups = false;
    static bool filter_high_cov = true;
    static bool filter_low_cov = true;
    static bool remove_internal_matches = false;
    static double max_overhang = 1000.0;
    static double max_overhang_ratio = 0.80;
    static bool remove_contained = false;
    static bool print_read_cov = false;
    static bool print_gse_stat = false;
    static bool keep_self_overlaps = false;
    static bool print_new_paf = false;
    static double min_iden = 0;
    static unsigned int min_match = 100;
}


bool endFile = false;
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

    // clear any previous log files with same name
    ofstream ofs;
    ofs.open( opt::sample_name + ".preqclr.log", ofstream::out | ios::trunc );
    ofs.close();

    out("========================================================");
    out("Run preqclr");
    out("========================================================");
    auto tot_start = chrono::system_clock::now();
    auto tot_start_cpu = clock();

    out("\n[ Parse reads file ]");
    auto swc = chrono::system_clock::now();
    auto scpu = clock();
    auto fq_records = parse_fq ( opt::reads_file );
    auto ewc = chrono::system_clock::now();
    auto ecpu = clock();
    fsec elapsedwc = ewc - swc;
    double elapsedcpu = (ecpu - scpu)/(double)CLOCKS_PER_SEC;;
    out("[+] Time elapsed: " + to_string(elapsedwc.count()) + "s, CPU time: "  + to_string(elapsedcpu) +"s");
 
    // parse the input PAF file and return a map with key = read id, 
    // and value = read object with only needed overlap info
    // SO: https://stackoverflow.com/questions/11062804/measuring-the-runtime-of-a-c-code
    // SO1 to get cast to milliseconds: https://stackoverflow.com/questions/30131181/calculate-time-to-execute-a-function 
    out("[ Parse PAF file ] ");
    swc = chrono::system_clock::now();    
    scpu = clock();
    map<string, sequence> paf_records = parse_paf();
    ewc = chrono::system_clock::now();
    ecpu = clock();
    elapsedwc = ewc - swc;
    elapsedcpu = (ecpu - scpu)/(double)CLOCKS_PER_SEC;;
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
    extern int optind, optopt;
    const char* const short_opts = ":g:c:hvr:n:p:";
    const option long_opts[] = {
        {"verbose",             no_argument,        NULL,   'v'},
        {"version",             no_argument,        NULL,   OPT_VERSION},
        {"reads",               required_argument,  NULL,   'r'},
        {"sample_name",         required_argument,  NULL,   'n'},
        {"paf",                 required_argument,  NULL,   'p'},
        {"gfa",                 required_argument,  NULL,   'g'},
        {"help",                no_argument,        NULL,   'h'},
        {"min-rlen",            required_argument,  NULL,   'l'},
        {"min-olen",            required_argument,  NULL,   'm'},
        {"keep-low-cov",        no_argument,        NULL,   OPT_KEEP_LOW_COV},
        {"keep-high-cov",       no_argument,        NULL,   OPT_KEEP_HIGH_COV},
        {"keep-dups",         no_argument,        NULL,   OPT_KEEP_DUPS},
        {"remove-int-matches",  no_argument,        NULL,   OPT_REMOVE_INT_MATCHES},
        {"remove-contained",    no_argument,        NULL,   OPT_REMOVE_CONTAINED},
        {"max-overhang",        required_argument,  NULL,   OPT_MAX_OVERHANG},
        {"max-overhang-ratio",  required_argument,  NULL,   OPT_MAX_OVERHANG_RATIO},
        {"print-read-cov",      no_argument,        NULL,   OPT_PRINT_READ_COV},
        {"print-gse-stat",      no_argument,        NULL,   OPT_PRINT_GSE_STAT},
        {"keep-self-overlaps",  no_argument,        NULL,   OPT_KEEP_SELF_OVERLAPS},
        {"print-new-paf",		no_argument,		NULL,	OPT_PRINT_NEW_PAF},
        { NULL, 0, NULL, 0 }
    };

    static const char* PREQCLR_CALCULATE_VERSION_MESSAGE =
    "preqclr " SUBPROGRAM " version " VERSION "\n"
    "Written by Joanna Pineda.\n"
    "\n"
    "Copyright 2018 Ontario Institute for Cancer Research\n";

    static const char* PREQCLR_CALCULATE_USAGE_MESSAGE =
    "Usage: preqclr [OPTIONS] --sample_name ecoli --reads reads.fa --paf overlaps.paf --gfa layout.gfa \n"
    "Calculate information for preqclr report\n"
    "\n"
    "    -v, --verbose		Display verbose output\n"
    "        --version		Display version\n"
    "    -r, --reads			Fasta, fastq, fasta.gz, or fastq.gz files containing reads\n"
    "    -n, --sample_name		Sample name; we recommend using the name of species for example\n" 
    "				This will be used as output prefix\n"
    "    -p, --paf			Minimap2 Pairwise mApping Format (PAF) file \n"
    "				This is produced using \'minimap2 -x ava-ont sample.fasta sample.fasta\'\n"
    "    -g, --gfa			Miniasm Graph Fragment Assembly (GFA) file\n"
    "				This file is produced using \'miniasm -f reads.fasta overlaps.paf\'\n"
    "    -l, --min-rlen=INT		Use overlaps with read lengths >= INT\n"
    "    -m, --min-olen=INT 		Use overlaps longer than >=INT\n"
    "        --keep-low-cov		Keep reads with low coverage (<= Q25 - IQR*1.25) for genome size est. calculations \n" 
    "        --keep-high-cov		Keep reads with high coverage (>= Q75 + IQR*1.25) for genome size est. calculations \n"
    "        --keep-self-overlaps	Keep overlaps where the query read and target read are the same \n"
    "        --keep-dups		Keep duplicate overlaps \n"
    "        --remove-contained	Remove contained overlaps \n"
    "        --remove-int-matches	Remove internal matches (overlaps where it is a short match in the middle of both reads) \n"
    "        --max-overhang		The maximum overhang length [default: 1000] \n"
    "        --max-overhang-ratio	The maximum overhang to mapping length ratio [default: 0.8] \n"
    "        --print-read-cov	Print read id and coverage for each read to stdout; overwrites verbose flag \n"
    "        --print-gse-stat	Print genome size estimate statistics only \n"
    "        --print-new-paf		Print new paf file after filtering overlaps\n"	
    "\n"
    "Report bugs to https://github.com/simpsonlab/preqclr/issues"
    "\n";

    int rflag=0, nflag=0, pflag=0, gflag=0;
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
            if ( rflag == 1 ) {
                fprintf(stderr, "preqclr: multiple instances of option -r,--reads. \n\n");
                fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE, argv[0]); 
                exit(1);
            }
            rflag = 1;
            arg >> opt::reads_file;
            break;
        case 'n':
            if ( nflag == 1 ) {
                fprintf(stderr, "preqclr: multiple instances of option -n,--sample_name. \n\n");
                fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE, argv[0]);
                exit(1);
            }
            nflag = 1;
            arg >> opt::sample_name;
            break;
        case 'p':
            if ( pflag == 1 ) {
                fprintf(stderr, "preqclr: multiple instances of option -p,--paf. \n\n");
                fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE, argv[0]);
                exit(1);
            }
            pflag = 1;
            arg >> opt::paf_file;
            break;
        case 'g':
            if ( gflag == 1 ) {
                fprintf(stderr, "preqclr: multiple instances of option -g,--gfa. \n\n");
                fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE, argv[0]);
                exit(1);
            }
            gflag = 1;
            arg >> opt::gfa_file;
            break;
        case 'h':
            cout << PREQCLR_CALCULATE_USAGE_MESSAGE << endl;
            exit(0);
        case 'l':
            arg >> opt::rlen_cutoff;
            break;
        case 'm':
            arg >> opt::olen_cutoff;
            break;
        case OPT_KEEP_LOW_COV:
            opt::filter_low_cov = false;
            break;
        case OPT_KEEP_HIGH_COV:
            opt::filter_high_cov = false;
            break;
        case OPT_KEEP_DUPS:
            opt::keep_dups = true;
            break;
        case OPT_REMOVE_INT_MATCHES:
            opt::remove_internal_matches = true;
            break;
        case OPT_MAX_OVERHANG:
            arg >> opt::max_overhang;
            opt::remove_internal_matches = true;
            break;
        case OPT_MAX_OVERHANG_RATIO:
            arg >> opt::max_overhang_ratio;
            opt::remove_internal_matches = true;
            break;
        case OPT_REMOVE_CONTAINED:
            opt::remove_contained = true;
            break;
        case OPT_PRINT_READ_COV:
            opt::print_read_cov = true;
            break;
        case OPT_PRINT_GSE_STAT:
            opt::print_gse_stat = true;
            break;
        case OPT_KEEP_SELF_OVERLAPS:
            opt::keep_self_overlaps = true;
            break;
        case OPT_PRINT_NEW_PAF:
            opt::print_new_paf = true;
            break;
        case '?':
            // invalid option: getopt_long already printed an error message
            if (optopt == 'c') {
                fprintf (stderr, "preqclr: option -%c requires an argument.\n", optopt);
            } else if(isprint(optopt)) {
                fprintf(stderr, "preqclr: invalid option thus ignored\n");
            }
            break;
        }
    }
   if (optind < argc) {
        printf("Non-option argument: ");
        while (optind < argc)
            printf("%s ", argv[optind++]);
        printf("\n");
    }

    // overwrite verbose flag if --print-read-cov or --print-gse-cov in use
    if (opt::print_read_cov) {
        opt::verbose = false;
    }
    if (opt::print_gse_stat) {
        opt::verbose = false;
        opt::print_read_cov = false;
    }
    if (opt::print_new_paf) {
        opt::verbose=false;
        opt::print_gse_stat=false;
        opt::print_read_cov=false;
    }

    // check mandatory variables and assign defaults
    if ( rflag == 0 ) {
        fprintf(stderr, "./preqclr: missing -r,--reads option\n\n");
        fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE, argv[0]);
        exit(1);
    }
    if ( nflag == 0 ) {
        fprintf(stderr, "preqclr: missing -n,--sample_name option\n\n");
        fprintf(stderr, PREQCLR_CALCULATE_USAGE_MESSAGE);
        exit(1);
    }
    if ( pflag == 0 ) {
        fprintf(stderr, "preqclr: missing -p,--paf option\n\n");
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
    const char *c1 = opt::paf_file.c_str();
    paf_file_t *fp1;
    paf_rec_t r1;
    fp1 = paf_open(c1);
    if (!fp1) {
        fprintf(stderr, "ERROR: PAF file failed to open. Check to see if it exists, is readable, and is non-empty.\n\n");
        exit(1);
    }
    
    // we need to filter overlaps
    // self overlap: minimap2 reports an overlap between the same reads (i.e. when query read id == target read id in the PAF file)
    // duplicate overlap: the same two reads are reported to overlap (i.e. when the same query read id and target read id pair are 
    //                    seen in multiple lines in the PAF file)
    vector<int> badlines; // stores all the lines we do not want
    map<size_t, pair<int, int>> h; // stores all the hashed query read name + target read name pairs with line number and alignment length
    int ln1 = 0; // current line number
    map<string, sequence> paf_records;
    while (paf_read(fp1, &r1) >= 0) {
        string qname = r1.qn;
        string tname = r1.tn;
        unsigned int qlen = r1.ql;
        unsigned int qstart = r1.qs;
        unsigned int qend = r1.qe;
        unsigned int strand = r1.rev;
        unsigned int tlen = r1.tl;
        unsigned int tstart = r1.ts;
        unsigned int tend = r1.te;
        unsigned int match = r1.ml;
        unsigned int al = r1.bl;
        double al_id = (double)match/(double)al;

        //cout << qname << "\t" << tname << "\t" << r1.bl << "\t" << ln1 << "\t" << bad << "\n";

        // start filtering overlaps
        // remove self overlaps
        if ( (!opt::keep_self_overlaps) && (qname.compare(tname) == 0) ) { 
              //self-overlap: the same read
              badlines.push_back(ln1);
              ln1++;
              continue;
        }

        // remove overlaps with low match id, length below cutoff
        if ( al_id < opt::min_iden || match < opt::min_match || al < opt::olen_cutoff || qlen < opt::rlen_cutoff || tlen < opt::rlen_cutoff ) {
            badlines.push_back(ln1);
            ln1++;
            continue;
        }

        // remove overlaps with high indel error rate
        /*int omax = max(qend - qstart, tend - tstart); 
        int omin = min(qend - qstart, tend - tstart);
        if ( (1 - double(omin/omax)) < 0.3 ) {
           // cout << qstart << "\t" << qend << "\t" << tstart << "\t" << tend << "\n";
           badlines.push_back(ln1);
           ln1++;
           continue;
        }*/

        // remove internal matches
        if ( opt::remove_internal_matches ) {
            // calculate overhang region
            unsigned int qprefix_len = qstart;
            unsigned int qsuffix_len = qlen - qend - 1;
            unsigned int tprefix_len = tstart;
            unsigned int tsuffix_len = tlen - tend - 1;
            int left_clip = 0, right_clip = 0;
            if ( ( qstart != 0 ) && ( tstart != 0 ) ){
                if ( strand == 0 ) { // if both on same strand
                    left_clip += min(qprefix_len, tprefix_len);
                } else {
                    left_clip += min(qprefix_len, tsuffix_len);
                }
            }
            if ( ( qend != 0 ) && ( tend != 0 ) ){
                if ( strand == 0 ) { // if both on same strand
                    right_clip += min(qsuffix_len, tsuffix_len);
                } else {
                    right_clip += min(qsuffix_len, tprefix_len);
                }
            }
            double maplen = double(max( qend - qstart, tend - tstart ))*opt::max_overhang_ratio;
            int overhang = left_clip + right_clip;

            if( opt::remove_internal_matches && ( double(overhang) > min(opt::max_overhang, maplen)) ) {
                // filter overlaps with long overhang regions
                badlines.push_back(ln1);
                ln1++;
                continue;
            }
        }

        // remove duplicate overlaps
        if ( !opt::keep_dups ) {
            // create a hashkey with lexicographically smallest combination of read names
            size_t hashkey = min(hash<string>{}(qname + tname), hash<string>{}(tname + qname));
            // check if we've seen this overlap between these two reads before
            auto it = h.find(hashkey);
            if (it != h.end()) {
                // YES, duplicate detected
                // let's compare the length of overlaps. We want the longer overlap.
                int curr_aln_len = int(r1.bl);
                int curr_ln = ln1;
                int prev_aln_len = int(it->second.first);
                int prev_ln = int(it->second.second);
                if ( curr_aln_len > prev_aln_len ) {
                    // prev. overlap between these 2 reads is shorter, we use the current line instead
                    // prev. overlap's line number is recorded as "bad". it will be skipped in second pass
                    badlines.push_back(prev_ln);
                    h[hashkey] = make_pair(curr_aln_len, curr_ln);
                    ln1++;
                    continue;
                    //cout << qname << "\t" << tname << "\t" << prev_ln << "\n";
                } else {
                    badlines.push_back(curr_ln);
                    ln1++;
                    continue;
                    //cout << qname << "\t" << tname << "\t" <<curr_ln << "\n";
                }
            } else {
                // First time we've seen this pair
                int aln_len = int(r1.bl);
                h.insert(make_pair(hashkey, make_pair(aln_len, ln1)));
            }
        }

        // adjust read length: read length = the region of read with overlaps only
        // store region with overlap on read and init read in paf_records
        auto i = paf_records.find(qname);
        bool success = true;
        if ( i == paf_records.end() ) {
            // if read not found initialize in paf_records
            sequence qr;
            qr.set(qlen, 0, qstart, qend);
            // cout << qname << "\t" << qr.min_s << "\t" << qr.max_e << "\n";
            paf_records.insert(pair<string,sequence>(qname, qr));
        } else {
            // if read found, update the overlap info
            success = i->second.updateOvlpRgn(qstart, qend);
            // cout << qname << "\t" << i->second.min_s << "\t" << i->second.max_e << "\n";
        }
        auto j = paf_records.find(tname);
        if ( i == paf_records.end() ) {
            // if read not found initialize in paf_records
            sequence tr;
            tr.set(tlen, 0, tstart, tend);
            paf_records.insert(pair<string,sequence>(tname, tr));
        } else {
            // if read found, update the overlap info
            success = j->second.updateOvlpRgn(tstart, tend);
        }

        if ( !success ){
            badlines.push_back(ln1);
        } 
        // next overlap, read next line!
        ln1+=1;    
    }
    h.clear(); // free up memory
    sort(badlines.begin(), badlines.end());

    // PASS 2: read each line in PAF file that has NOT been noted in PASS 1 as a "bad" line
    const char *c2 = opt::paf_file.c_str();
    paf_file_t *fp2;
    paf_rec_t r2;
    fp2 = paf_open(c2);
    if (!fp2) {
        fprintf(stderr, "ERROR: PAF file failed to open. Check to see if it exists, is readable, and is non-empty.\n\n");
        exit(1);
    }
    int ln2 = 0;
    unsigned int iv = 0; // index in vector
    // the vector is sorted numerically
    // we can loop through the vector once by storing which is the next line to avoid
    // once we have reached this line, we can move on to the next bad line and
    // look out for that one while going through the next lines
    int bad;
    if ( !badlines.empty() ) {
        bad = int(badlines.at(iv)); // first bad line to watch out for
    } else {
        bad = ln1; // no badlines detected, set to last line
    }
    // read good lines in PAF
    while (paf_read(fp2, &r2) >= 0) { 
        if ( ln2 < bad ) {
            // read each line/overlap and save each column into variable
            string qname = r2.qn;
            string tname = r2.tn;
            auto i = paf_records.find(qname);
            auto j = paf_records.find(tname);
            unsigned int qlen = r2.ql;
            unsigned int qalen = i->second.max_e - i->second.min_s;
            unsigned int qstart = r2.qs;
            unsigned int qend = r2.qe;
            unsigned int strand = r2.rev;
            unsigned int tlen = r2.tl;
            unsigned int talen = j->second.max_e - j->second.min_s;
            unsigned int tstart = r2.ts;
            unsigned int tend = r2.te;
            if ( qalen > opt::rlen_cutoff && talen > opt::rlen_cutoff && double(tlen-talen)/tlen < 0.10 && double(qlen-qalen)/qlen < 0.10  ) {
                if ( opt::print_new_paf) {
                    string s = ( strand == 0 ) ? "-" : "+";
                    cout << qname << "\t" << qalen << "\t" << qstart << "\t" << qend << "\t" << s <<"\t" << tname << "\t" << talen << "\t" << tstart << "\t" << tend << "\t" << r2.ml << "\t"<< r2.bl << "\t255\n";
                }

                // calculate overlap length, we need to take into account minimap2's softclipping 
                // adjust to new read length (region with overlaps only)
                unsigned int qprefix_len = qstart - i->second.min_s;
                unsigned int qsuffix_len = i->second.max_e - qend;
                unsigned int tprefix_len = tstart - j->second.min_s;
                unsigned int tsuffix_len = j->second.max_e - tend;

                int left_clip = 0, right_clip = 0;
                if ( ( qstart != 0 ) && ( tstart !=0 )) {
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
                int overhang = left_clip + right_clip;
                int omax = max(qend - qstart, tend - tstart);
                int omin = min(qend - qstart, tend - tstart);
                int gaps = omax - omin;
                // calculate coverage per read               
                unsigned int qoverlap_len = abs(qend - qstart) - gaps + overhang;
                double qcov = double(qoverlap_len) / double(qalen);
                i->second.updateCov(qcov);
                unsigned int toverlap_len = abs(tend - tstart) - gaps + overhang;
                double tcov = double(toverlap_len) / double(talen); 
                j->second.updateCov(tcov);
            }
        } else if ( iv < badlines.size()-1 ) {
            // we have a bad line!
            // next bad line to look out for:
            iv+=1;
            bad = int(badlines.at(iv));
        }
        ln2+=1;
    }
    if ( opt::print_read_cov ) {
        for ( auto const& r : paf_records ) {
            sequence temp = r.second;
            cout << r.first << "\t" << temp.read_len << "\t" << temp.cov << "\n";
        }
    }

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
    return fq_records;
}

//
// Dust scoring scheme as given by:
// Morgulis A. "A fast and symmetric DUST implementation to Mask
// Low-Complexity DNA Sequences". J Comp Bio.
double calculateDustScore(const string& seq)
{
    map<string, int> scoreMap;

    // Cannot calculate dust scores on very short reads
    if(seq.size() < 5)
        return 0.0f;


    // Slide a 3-mer window over the sequence and insert the sequences into the map
    for(size_t i = 0; i < seq.size() - 5; ++i)
    {
        string fiveMer = seq.substr(i, 5);
        scoreMap[fiveMer]++;
    }

    // Calculate the score by summing the square of every element in the map
    double sum = 0;
    for (auto& iter : scoreMap) {
        int tc = iter.second;
        double score = (double)(tc * (tc - 1)) / 2.0f;
        sum += score;
    }
    return sum / (seq.size() - 4);
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
        int r_len = r.max_e - r.min_s;
        long double r_cov = r.cov;
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

    // get the index of the 25th and 75th percentile item
    int i25 = ceil(covs.size() * 0.25);
    int i75 = ceil(covs.size() * 0.75); 
    double IQR = covs[i75].first - covs[i25].first;
    double bd = IQR*1.5;
    double upperbound = round(double(covs[i75].first) + bd);
    double lowerbound = (round(double(covs[i25].first) - bd)>3.0) ? round(double(covs[i25].first) - bd) : 3.0;
    if ( !opt::filter_low_cov ) {
        lowerbound = 3.0;
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
    //double l = covs[0].first;
    //double u = covs[0].first + 0.25;
    double l = lowerbound;
    double u = l + 0.25;
    double curr_largest = -1000.0;
    double mode_cov = 0;
    int count = 0;
    unsigned int i = 0;
    // we want to only consider reads above the lowerbound
    while ( covs[i].first < l ) {
        i += 1;
    }
    while ( i < covs.size() ){ 
        // iterate through reads
        // look at reads that fall within current bin
        while ( covs[i].first >= l && covs[i].first < u ){ 
            // filter outliers: [Q25-IQR*1.5, Q75+IQR*1.5]
            if ( (covs[i].first >= lowerbound) && (covs[i].first <= upperbound) ){
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
        if (( count > curr_largest ) && ( u > lowerbound )) {
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
    if ( opt::print_gse_stat ) {
        cout <<  opt::sample_name << "\t" << mode_cov << "\t" <<  median_cov << "\t" << mean_read_len << "\t" << tot_reads  << "\t" << tot_reads_f << "\t" << sum_len << "\t" << sum_len_f << "\t"<< est_genome_size << "\t" <<  est_genome_size1 << "\n";
    }
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
