#include "BSEQE.h"
#include "list.h"
using namespace std;

int main() {
    double d = 0.1;
    vector<double> vec;
    string fn;
    for(double q = 1+d; q < 300; q+=d ) {
        vec.clear();
        for(int m=2; m<200; m++) {
            cout << m << endl;
            vec.push_back(derror(m, q, 0, 2, 20000/m));
        }
        fn = "data/m_"+to_string(q)+".list";
        list_write(vec, fn.c_str());
    }
    return 0;
}