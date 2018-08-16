//---------------------------------------------------------
// Copyright 2018 Ontario Institute for Cancer Research
// Written by Joanna Pineda (joanna.pineda@oicr.on.ca)
//---------------------------------------------------------
//
// preqclr contig -- holds contig information calculated from miniasm
//
#include <string>

using namespace std;

class contig
{
  public:
    int len, num_reads;
    void set(int l, int n);
};
