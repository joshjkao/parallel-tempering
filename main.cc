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

using namespace std;
using namespace std::chrono;


int main(int argc, char** argv) {
    cout << "starting program" << endl;
    auto start = high_resolution_clock::now();

    // vector<double> T = {0.10,1.40,1.65,1.85,2.04,2.18,2.26,2.27,2.28,2.29,2.32,2.37,2.44,2.55,2.70,2.90,3.25,3.70,4.60,6.20,10.0}; // opt
    vector<double> T = {0.10,0.13,0.16,0.20,0.25,0.32,0.40,0.50,0.63,0.79,1.00,1.26,1.58,2.00,2.51,3.16,3.98,5.01,6.31,7.94,10.0};
    // vector<double> T = {10,4,2,1.1,.1};
    // vector<double> T = {1/0.1,1/0.379576,1/0.406646,1/0.518733,1/10.};
    // vector<double> T = {1/0.1,1,1/10.};
    // vector<double> B = {.1, 1, 2, 5};
    vector<double> B;
    for (auto& t: T) {
        B.push_back(pow(t, -1));
    }

    const int n_sweeps = 1e6;
    const int L = 20;

    int n_reps = B.size();

    vector<double> p;
    double tau = DBL_MAX;

    pcg32 engine(0);
    
    // for (int r = 0; r < 20; ++r) {
        auto reset = PT::construct_replicas(L, B, engine);
        p = PT::start(n_sweeps, reset, engine);
        tau = PT::expected_rt(p);
        PT::delete_replicas(reset);
        cout << tau << " ------ " << B;
        cout << p << endl;
        // int i = 0;
    //     #pragma omp parallel for
    //     for (i = 0; i < 200; ++i) {
    //         // setup the RNG
    //         mt19937_64 gen(0);

    //         // create a copy of the current best arrangement, only 1 thread at once
    //         vector<double> Bcpy;
    //         #pragma omp critical
    //         {
    //             Bcpy = B;
    //             gen.seed(RNG::uniform_int(engine));
    //         }

    //         // adjust 1 beta value in new array
    //         auto op = RNG::pick_op(gen);
    //         switch (op) {
    //             case 0: {
    //                 PT::insertion(Bcpy, gen);
    //                 break;
    //             }
    //             case 1: {
    //                 PT::deletion(Bcpy, gen);
    //                 break;
    //             }
    //             default: {
    //                 PT::internal_adj(Bcpy, gen);
    //                 break;
    //             }
    //         }

    //         // construct replica array
    //         auto reps = PT::construct_replicas(L, Bcpy, gen);
            
    //         // run the simulation and calculate trip time
    //         auto p_new = PT::start(n_sweeps, reps, gen);
    //         double tau_new = PT::expected_rt(p_new);

    //         // check if the new trip is better, only 1 thread at once
    //         #pragma omp critical
    //         {
    //             if (tau_new < tau) {
    //                 tau = tau_new;
    //                 p = p_new;
    //                 B = Bcpy;
    //                 cout << tau << endl;
    //             }
    //         }

    //         // delete replica array
    //         PT::delete_replicas(reps);
    //     }
    // }

    cout << endl << B << p << tau << endl;

    for (int j = 0; j < n_reps; ++j) {
        T[j] = 1/B[j];
    }

    write_to_file(T, p, "output.csv");

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "\nexecution time: " << duration.count()/1e6 << " s\n";
}