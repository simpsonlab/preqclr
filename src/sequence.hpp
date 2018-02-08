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
    string read_id;
    unsigned long int read_len;
    double cov;
    void set(string i, unsigned long int l, double c);
    void updateCov(double c);
};
