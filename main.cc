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
    // vector<double> T = {0.10,0.13,0.16,0.20,0.25,0.32,0.40,0.50,0.63,0.79,1.00,1.26,1.58,2.00,2.51,3.16,3.98,5.01,6.31,7.94,10.0};
    // vector<double> T = {.1,1.905,2.386,3.12,5};
    vector<double> T = {.1, 10};
    // 0.2,0.320761,0.418838,0.525278,10
    const int L = 10;

    cout << "starting optimization with T = " << T << "and L = " << L << endl;

    // convert to inverse temperature
    vector<double> B;
    for (auto& t: T) {
        B.push_back(pow(t, -1));
    }
    int n_reps = B.size();
    auto B_best = B;
    
    // return variables
    vector<double> p;
    double tau = INFINITY;

    // RNG setup
    unsigned int seed = time(0);
    pcg32 main_engine(seed);
    pcg32 local_engine(0);

    // max number of MC sweeps before timing out PT
    const int n_sweeps = 1e6;

    // max number of reps before terminating
    const int max_reps = 25;

    // number of optimization steps before recalculating tau
    const int n_adj = 200;

    // max number of resets before timing out optimization
    const int n_reset = 100;

    // highest fraction of old/new mavg to continue the loop
    const double fail_threshold = .01;

    // exponential mavg constant
    const double alpha = .20;

    // minimum number of consecutive failed updates to break the loop
    const int min_fails = 10;

    // ratio of failed to total updates to break the loop
    const double fail_ratio = 0.75;


    string filename = to_string(seed) + string(".txt");
    ofstream mavgout(filename);
    mavgout << "mavg" << endl;


    // multithreading
    int n_threads = 1;
    #ifdef USE_MULTITHREADING
    n_threads = omp_get_max_threads();
    #endif

    // start the timer
    auto start = high_resolution_clock::now();

    for (int reps = 0; reps < max_reps; ++reps) {
        auto B_this_n = B;
        auto tau_this_n = tau;
        auto p_this_n = p;

        double mavg = 0.;

        int stop_counter = 0;
        int r;

        // make continuous adjustments for this number of replicas
        for (r = 0; r < n_reset && stop_counter < min_fails + r*fail_ratio; ++r) {
            
            // recalculate expecte probabilities of the current best
            auto reset = PT::construct_replicas(L, B_this_n, main_engine);
            p_this_n = PT::start(n_sweeps, reset, main_engine);
            tau_this_n = PT::expected_rt(p_this_n);
            PT::delete_replicas(reset);
            cout << tau_this_n << " ------ " << B_this_n;
            cout << p_this_n << endl;

            // adjust and accept if better
            int i = 0;
            #pragma omp parallel for firstprivate(local_engine) num_threads(n_threads)
            for (i = 0; i < n_adj; ++i) {
                bool thread_init = false;

                // create a copy of the current best arrangement, thread protected
                vector<double> Bcpy;
                #pragma omp critical
                {
                    Bcpy = B_this_n;
                    if (!thread_init) {
                        local_engine.seed(RNG::uniform_int(main_engine));
                        thread_init = true;
                    }
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
                    if (tau_new < tau_this_n) {
                        tau_this_n = tau_new;
                        p_this_n = p_new;
                        B_this_n = Bcpy;
                    }
                }
            }
            if (tau_this_n == INFINITY) ++stop_counter;
            else if (mavg == 0.) mavg = tau_this_n;
            else {
                auto old_mavg = mavg;
                mavg = mavg*(1-alpha) + tau_this_n*alpha;
                cout << mavg << endl;
                mavgout << mavg << "            " << B_this_n;
                if ((mavg - old_mavg)/old_mavg < fail_threshold) ++stop_counter;
                else stop_counter = 0;
            }
        }
        if (r == n_reset) cout << "reached max optimizing iterations\n";
        if (tau_this_n <= tau) {
            tau = tau_this_n;
            B = B_best = B_this_n;
            p = p_this_n;
            PT::insertion(B, main_engine);
        }
        else {
            break;
        }
    }

    
    cout << endl << B_best << p << tau << endl;

    for (int j = 0; j < n_reps; ++j) {
        T[j] = 1/B[j];
    }

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "\nexecution time: " << duration.count()/1e6 << " s\n";
}