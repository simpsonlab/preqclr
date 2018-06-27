//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Joanna Pineda (joanna.pineda@oicr.on.ca)
//---------------------------------------------------------
//
// preqc-lr sequence -- holds read information calculated from overlaps
//
#include <string>

using namespace std;

class sequence
{
  public:
    unsigned long int read_len;
    long double cov;
    void set(unsigned long int l, long double c);
    void updateCov(long double c);
};
