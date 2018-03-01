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

using namespace std;

void sequence::set(string i, unsigned long int l, double c)
{
    read_id = i;
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

