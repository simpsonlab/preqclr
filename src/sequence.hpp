//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Joanna Pineda (joanna.pineda@oicr.on.ca)
//---------------------------------------------------------
//
// preqc-lr sequence -- holds read information calculated from overlaps
//
#include <string>

using namespace std;

typedef std::string Seq;



class sequence
{
  public:
    unsigned long int read_len;
    double cov;
    void set_paf(string qname,string tname,unsigned int qlen,\
                unsigned int qstart, unsigned int qend,unsigned int strand,\
                unsigned int tlen,unsigned int tstart,unsigned int tend);
    void set(unsigned long int l, double c);
    void updateCov(double c);
    string qname;
    string tname;
    unsigned int qlen;
    unsigned int qstart;
    unsigned int qend;
    unsigned int strand;
    unsigned int tlen;
    unsigned int tstart;
    unsigned int tend;
    double dv;
};


Seq reverseComplement(const Seq& s);

