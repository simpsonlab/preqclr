#include "read_seq.hpp"
#include <string>
#include <iostream>

using namespace std;

void read_seq::set(string i, unsigned long int l, double c)
{
    read_id = i;
    read_len = l;
    cov = c;
}

void read_seq::updateCov( double c )
{
    cov += c;
}
