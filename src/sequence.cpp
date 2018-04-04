//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Joanna Pineda (joanna.pineda@oicr.on.ca)
//---------------------------------------------------------
//
// preqc-lr sequence -- holds read information calculated from overlaps
//
#include "sequence.hpp"
#include <string>
#include <iostream>
#include <algorithm>
#include <cassert>

using namespace std;

void sequence::set(unsigned long int l, double c)
{
    read_len = l;
    cov = c;
}

void sequence::updateCov( double c )
{
    cov += c;
}


void sequence::set_paf(string qname,string tname,unsigned int qlen,\
                unsigned int qstart, unsigned int qend,unsigned int strand,\
                unsigned int tlen,unsigned int tstart,unsigned int tend)
{

    this->qname = qname;
    this->tname = tname;
    this->qlen = qlen; 
    this->qstart = qstart;
    this->qend = qend;
    this->strand = strand;
    this->tlen = tlen;
    this->tstart = tstart;
    this->tend = tend;               

}

/* Return the complement of a given nucleotide. */
char complementBaseChar(char c)
{
    char rc;
    switch (toupper(c)) {
    case 'A':
        rc = 'T';
        break;
    case 'C':
        rc = 'G';
        break;
    case 'G':
        rc = 'C';
        break;
    case 'T':
        rc = 'A';
        break;
    case 'N':
        rc = 'N';
        break;
    default:
        cerr << "error: unexpected character: `" << c << "'\n";
        assert(false);
        abort();
    }
    return islower(c) ? tolower(rc) : rc;
}


/* Return the reverse complement of a given sequence. */
Seq reverseComplement(const Seq& s)
{
    Seq rc(s);
    reverse(rc.begin(), rc.end());
    transform(rc.begin(), rc.end(), rc.begin(),
              complementBaseChar);
    return rc;
}

