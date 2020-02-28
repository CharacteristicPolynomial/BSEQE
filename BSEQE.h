#pragma once
#include <iostream>
#include <vector>
#include <random>
#include <chrono>
#include <cmath>
#include <cstdlib>
using namespace std;

double quant_fac(double q, const vector<double> x, double k, int ne=1000) {
    int m = x.size();
    // calculate the quantization factor
    // first construct a categarical random variable
    vector<double> ths(m); // so if a random number is between ths[m-1] and ths[m], then it is categray m (with ths[-1]=0)
    vector<double> tempv(m);
    double temp2 = 0;
    for(int j=0; j<m; j++) {
        temp2 += pow(1.0/q, x[j]);
        tempv[j] = temp2;
    }
    for(int j=0; j<m; j++) {
        ths[j] = tempv[j]/temp2*RAND_MAX;
    }

    int mb = 2*m; // number of summations
    double nomin =0;
    double denomin = 0;
    vector<double> ex(m); // x for experiment
    for(int i=0; i<ne; i++) {
        // cout << i << endl;
        double qj = 1;
        double mj1 = 1;
        double mj0 = 1;
        double phij = 1;
        fill(ex.begin(), ex.end(), 0);
        for (int j=0; j<mb; j++) {
            // qj for (q-1)^j/q^j
            // mj1 for {j+m-k-1 \choose j}
            // mj0 for {j+m-k \choose j}
            // phij for \phi(k;p)
            nomin += qj * mj1 * phij;
            denomin += qj * mj0 * phij;

            double draw = rand(); // draw your card!
            int di=-1;
            for(int t=0; t<m; t++) {
                if(draw < ths[t]) {
                    di = t;
                    break;
                }
            }
            if(di == -1) {
                // cerr << "what?" << endl;
                // exit(-1);
                di = m-1;
            }

            // update
            qj *= (q-1.0) / q;
            mj1 *= (double)(m-k+j+2)/(j+1);
            mj0 *= (double)(j+m-k+1)/(j+1);
            phij *= (double)(ex[di]+1)/(ex[di]+2);
            ex[di] += 1;
        }
    }
    return nomin / denomin;
}

double bseqe_q(const vector<double>& x, double q, double k) {
    // BSEQE g(q^{-x},q^{-x-1}) with prior parameter k
    // REQUIRE: x.size = m
    int m = x.size();
    double conte = 0; // the continuous estimator
    for(int j=0; j<m; j++) {
        conte += pow(1.0/q, x[j]);
    }
    conte = (double)(m+1-k)  / conte ;
    
    // the quantization factor
    // cout << nomin / denomin << endl;
    return conte * quant_fac(q,x,k);
}

double derror(int m, double q, double t, int k, int ne = 100000) {
    // REQUIRE: m>=1, q>1, t>=0, k>=1
    double lamb = 1+(q-1)*t;
    exponential_distribution<double> dist(lamb);
    default_random_engine gen;
    unsigned d = chrono::high_resolution_clock::now().time_since_epoch().count();
    gen.seed(d);
    
    vector<double> x(m);
    double ser = 0; // sum of errors
    for (int i=0; i<ne; i++) {
        for(int j=0; j<m; j++) {
            double roll = dist(gen);
            double xj = 0;
            if(roll > 1) {
                while(roll > 1) {
                    xj -= 1;
                    roll /= q;
                }
            } else {
                roll *= q;
                while(roll <= 1) {
                    xj += 1;
                    roll *= q;
                }
            }
            x[j] = xj;
        }
        double est = bseqe_q(x, q, k);
        // cout << i << " : " <<  est << endl;
        ser += (est/lamb - 1.0) * (est/lamb - 1.0);
    }
    return ser/ne;
}

double H(double q, double lamb) {
    // in unit of bits
    double ent = 0;
    double l = 10;
    int a = -l/2.0/log(q);
    int b = l/log(q);
    for(int i=a; i<b; i++) {
        double eqx1 = exp(-lamb*pow(1.0/q,i+1));
        double eqx = exp(-lamb*pow(1.0/q,i));
        ent -= (eqx1-eqx)*log2(eqx1-eqx);
    }
    return ent;
}

double D(double q, double lamb) {
    double d = 0;
    double l = 10;
    int a = -l/2.0/log(q);
    int b = l/log(q);
    for(int i=a; i<b; i++) {
        double qx1 = -lamb*pow(1.0/q,i+1);
        double qx = -lamb*pow(1.0/q,i);
        double eqx1 = exp(-lamb*pow(1.0/q,i+1));
        double eqx = exp(-lamb*pow(1.0/q,i));
        d += (qx1*eqx1-qx*eqx)*(qx1*eqx1-qx*eqx)/(eqx1-eqx);
    }
    return 1.0/d;
}