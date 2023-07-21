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
    // vector<double> T = {0.10,0.13,0.16,0.20,0.25,0.32,0.40,0.50,0.63,0.79,1.00,1.26,1.58,2.00,2.51,3.16,3.98,5.01,6.31,7.94,10.0};
    vector<double> T = {.1,.6, 1, 2, 5};

    vector<double> B;
    for (auto& t: T) {
        B.push_back(pow(t, -1));
    }

    const int n_sweeps = 1e6;
    const int L = 9;

    int n_reps = B.size();

    vector<double> p;
    double tau = DBL_MAX;


    unsigned int seed = time(0);
    pcg32 main_engine(seed);
    pcg32 local_engine(0);

    int failed_iterations = 0;
    
    int r;
    for (r = 0; r < 20; ++r) {
        auto reset = PT::construct_replicas(L, B, main_engine);
        p = PT::start(n_sweeps, reset, main_engine);
        tau = PT::expected_rt(p);
        PT::delete_replicas(reset);
        cout << tau << " ------ " << B;
        cout << p << endl;

        int i = 0;
        #pragma omp parallel for firstprivate(local_engine) num_threads(4)
        for (i = 0; i < 200; ++i) {

            // create a copy of the current best arrangement, only 1 thread at once
            vector<double> Bcpy;
            #pragma omp critical
            {
                Bcpy = B;
                local_engine.seed(RNG::uniform_int(main_engine));
            }

            // adjust 1 beta value in new array
            PT::internal_adj(Bcpy, local_engine);

            // construct replica array
            auto reps = PT::construct_replicas(L, Bcpy, local_engine);
            
            // run the simulation and calculate trip time
            auto p_new = PT::start(n_sweeps, reps, local_engine);
            double tau_new = PT::expected_rt(p_new);

            // delete replica array
            PT::delete_replicas(reps);

            // check if the new arrangement is better, only 1 thread at once
            #pragma omp critical
            {
                if (tau_new < tau) {
                    tau = tau_new;
                    p = p_new;
                    B = Bcpy;
                    cout << tau << endl;
                    failed_iterations = 0;
                }
                else ++failed_iterations;
            }
        }
        if (failed_iterations > i*200*.75) {
            cout << "terminated early\n";
            break;
        }
    }
    if (r == 20) cout << "reached maxed optimizing iterations\n";


    cout << endl << B << p << tau << endl;

    for (int j = 0; j < n_reps; ++j) {
        T[j] = 1/B[j];
    }

    write_to_file(T, p, "output.csv");

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "\nexecution time: " << duration.count()/1e6 << " s\n";
}