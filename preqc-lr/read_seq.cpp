#include "read_seq.hpp"
#include <string>
#include <iostream>

using namespace std;

void read_seq::set(string i, int l, int lol)
{
    read_id = i;
    read_len = l;
    total_len_overlaps = lol;
    total_num_overlaps = 1;
}

void read_seq::updateOverlap(int new_ol_len)
{
    total_len_overlaps+=new_ol_len;
    total_num_overlaps+=1;
}
