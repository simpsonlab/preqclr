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

void sequence::set(unsigned long int l, long double c)
{
    read_len = l;
    cov = c;
}

void sequence::updateCov( long double c )
{
    cov += c;
}
