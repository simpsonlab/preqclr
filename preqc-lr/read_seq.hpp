#include <string>

using namespace std;

class read_seq
{
  public:
    string read_id;
    unsigned long int read_len;
    double cov;
    void set(string i, unsigned long int l, double c);
    void updateCov(double c);
};
