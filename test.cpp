#include "BSEQE.h"
using namespace std;

int main() {
    for(double q=1; q<100; q+=1) {
        cout << q << " : " << H(q,1+(q-1)/2) << " " << D(q,1+(q-1)/2) << " " << H(q,1+(q-1)/2)*D(q,1+(q-1)/2) << endl;
    }
    for(double q=1; q<100; q+=1) {
        cout << q << " : " << H(q,1) << " " << D(q,1) << " " << H(q,1)*D(q,1) << endl;
    }

    return 0;
}