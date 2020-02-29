#include "BSEQE.h"
#include "list.h"
using namespace std;

int main() {
    for(int i=0; i<4; i++) {
        cout << "q=" << 2 << "; m=" << 100 << "; i=" << i << endl;
        cout << "0 offset: " << derror(100, 2, (double)i/4.0, 2, false, 10000) << endl;
        cout << "unfirom offset: " << derror(100, 2, (double)i/4.0, 2, true, 10000) << endl;
    }
    for(int i=0; i<4; i++) {
        cout << "q=" << 5 << "; m=" << 164 << "; i=" << i << endl;
        cout << "0 offset: " << derror(164, 5, (double)i/4.0, 2, false, 10000) << endl;
        cout << "unfirom offset: " << derror(164, 5, (double)i/4.0, 2, true, 10000) << endl;
    }
    for(int i=0; i<4; i++) {
        cout << "q=" << 9 << "; m=" << 208 << "; i=" << i << endl;
        cout << "0 offset: " << derror(208, 9, (double)i/4.0, 2, false, 10000) << endl;
        cout << "unfirom offset: " << derror(208, 9, (double)i/4.0, 2, true, 10000) << endl;
    }
    return 0;
}