#include <string>

using namespace std;

class read_seq
{
  public:
    string read_id;
    int read_len;
    int total_len_overlaps;
    int total_num_overlaps;
    void set(string i, int l, int lol);
    void updateOverlap(int new_ol_len);
};
