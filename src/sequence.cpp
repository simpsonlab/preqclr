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

void sequence::set(unsigned long int l, double c, int s, int e)
{
    read_len = l;
    cov = c;
    min_s = s;
    max_e = e;
}

void sequence::updateCov(double c )
{
    cov += c;
}

bool sequence::updateOvlpRgn(int s, int e) 
{
    if ( (s < (max_e + 200)) && (min_s < (e + 200))  ) {
        if ( e > max_e ) {
            max_e = e;
        }
        if ( s < min_s ) {
            min_s = s;
        }
        return true;
    } else {
        return false;
    }
}
