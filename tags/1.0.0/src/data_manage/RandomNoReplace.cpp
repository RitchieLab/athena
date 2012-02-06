#include "RandomNoReplace.h"

using namespace data_manage;


///
/// Algorithm 3.4.2S of Knuth's book Seminumeric Algorithm
///
void RandomNoReplace::SampleWithoutReplacement(int popSize, int sampSize, std::vector<int>& samples){
    // Use Knuth's variable names
    int& n = sampSize;
    int& N = popSize;

    int t = 0; // total input records dealt with
    int m = 0; // number of items selected so far
    double u;

    while (m < n){
        u = double(rand()) / RAND_MAX; // call a uniform(0,1) random number generator

        if ( (N - t)*u >= n - m ){
            t++;
        }
        else{
            samples[m] = t;
            t++; m++;
        }
    }
}
