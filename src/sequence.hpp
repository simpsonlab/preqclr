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
    double cov;
    int min_s;
    int max_e;
    void set(unsigned long int l, double c, int s, int e);
    void updateCov(double c);
    bool updateOvlpRgn(int s, int e);
};

