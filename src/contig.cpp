//---------------------------------------------------------
// Copyright 2017 Ontario Institute for Cancer Research
// Written by Joanna Pineda (joanna.pineda@oicr.on.ca)
//---------------------------------------------------------
//
// preqclr contig -- holds contig information calculated with miniasm
//
#include "contig.hpp"
#include <string>
#include <iostream>

using namespace std;

void contig::set(int l, int n)
{
    len = l;
    num_reads = n;
}

