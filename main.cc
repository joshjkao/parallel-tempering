#include "Replica.hh"
#include "PT.hh"
#include "util.hh"
#include "RNG.hh"
#include "pcg/pcg_random.hpp"

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <climits>
#include <chrono>
#include <fstream>
#include <float.h>
#include <omp.h>
#include <algorithm>
#include <string>

using namespace std;
using namespace std::chrono;


int main(int argc, char** argv) {
    // vector<double> T = {0.10,1.40,1.65,1.85,2.04,2.18,2.26,2.27,2.28,2.29,2.32,2.37,2.44,2.55,2.70,2.90,3.25,3.70,4.60,6.20,10.0}; // opt
    vector<double> T = {0.10,0.13,0.16,0.20,0.25,0.32,0.40,0.50,0.63,0.79,1.00,1.26,1.58,2.00,2.51,3.16,3.98,5.01,6.31,7.94,10.0};
    // vector<double> T = {.1,1.905,2.386,3.12,5};
    // vector<double> T = {.1, 10};
    // 0.2,0.320761,0.418838,0.525278,10
    int L = 20;

    cout << "starting optimization with T = " << T << "and L = " << L << endl;

    // convert to inverse temperature
    vector<double> B;
    for (auto& t: T) {
        B.push_back(pow(t, -1));
    }

    // start the timer
    auto start = high_resolution_clock::now();


    auto pt = PT(L, B, 1);
    auto p = pt.Start(1e5);
    cout << p ;
    
    
    





    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "\nexecution time: " << duration.count()/1e6 << " s\n";
}