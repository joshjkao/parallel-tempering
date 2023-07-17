#include "Replica.hh"
#include "PT.hh"
#include "util.hh"

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <climits>
#include <chrono>
#include <fstream>
#include <float.h>
#include <omp.h>

using namespace std;
using namespace std::chrono;


int main(int argc, char** argv) {
    // srand(time(0));
    srand(0);
    auto start = high_resolution_clock::now();

    // vector<double> T = {0.10,1.40,1.65,1.85,2.04,2.18,2.26,2.27,2.28,2.29,2.32,2.37,2.44,2.55,2.70,2.90,3.25,3.70,4.60,6.20,10.0}; // opt
    // vector<double> T = {0.10,0.13,0.16,0.20,0.25,0.32,0.40,0.50,0.63,0.79,1.00,1.26,1.58,2.00,2.51,3.16,3.98,5.01,6.31,7.94,10.0};
    vector<double> T = {0.10, 1.1, 2, 4, 10};
    // vector<double> B = {.1, 1, 2, 5};
    vector<double> B;
    for (auto& t: T) {
        B.push_back(1/t);
    }

    int n_sweeps = 1e5;
    int L = 10;

    int n_reps = B.size();

    vector<double> p;
    double tau = DBL_MAX;
    
    for (int r = 0; r < 20; ++r) {
        vector<Replica*> reset(n_reps, nullptr);
        for (int j =0; j < n_reps; ++j) {
            reset[j] = new IsingFerromagnetReplica(L, B[j]);
        }
        tau = PT::expected_rt(PT::start(n_sweeps, reset));
        for (auto& res: reset) {
            delete res;
        }


        int i = 0;
        #pragma omp parallel for 
        for (i = 0; i < 200; ++i) {
            // create a copy of the current best arrangement, only 1 thread at once
            vector<double> Bcpy;
            #pragma omp critical
            {
            Bcpy = B;
            // cout << Bcpy;
            }

            // adjust 1 beta value in new array, update Bcpy
            int adj = rand()%(n_reps-2)+1;
            // cout << adj << " ";
            double Bl = Bcpy[adj-1];
            double Br = Bcpy[adj+1];
            Bcpy[adj] = (double)rand()/INT_MAX * (Br-Bl) + Bl;

            // construct replica array
            vector<Replica*> reps(n_reps, nullptr);
            for (int j =0; j < n_reps; ++j) {
                reps[j] = new IsingFerromagnetReplica(L, Bcpy[j]);
            }
            
            // run the simulation and calculate trip time
            auto p_new = PT::start(n_sweeps, reps);
            double tau_new = PT::expected_rt(p_new);

            // check if the new trip is better, only 1 thread at once
            #pragma omp critical
            {
            if (tau_new < tau) {
                tau = tau_new;
                p = p_new;
                B = Bcpy;
                cout << tau << endl;
            }
            }

            // delete replica array
            for (auto& rep: reps) {
                delete rep;
            }
        }
    }

    cout << endl << B << p << tau << endl;

    for (int j = 0; j < n_reps; ++j) {
        T[j] = 1/B[j];
    }

    write_to_file(T, p, "output.csv");

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "\nexecution time: " << duration.count()/1e6 << " s\n";
}