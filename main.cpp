#include "Replica.h"
#include "PT.h"
#include "RNG.h"
#include "util.h"

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <climits>
#include <chrono>
#include <fstream>
#include <float.h>
// #include <omp.h>
#include <algorithm>
#include <string>

using namespace std;
using namespace std::chrono;

void original_optimize(int L, vector<double> B, const char* filename);
void optimize(int L, vector<double> B, const char* filename);

void test(int L, vector<double> B);


int main(int argc, char** argv) {

    // parse cmd line arguments
    int L = 7;
    // if (argc < 2 || atoi(argv[1]) < 3) {
    //     cout << "error with command line args\n";
    //     //return 1;
    // } else {
    //     L = atoi(argv[1]);
    // }

    // start the timer
    auto start = high_resolution_clock::now();

    // perform optimization
    //vector<double> B = {10, .1};
    //optimize(L, B, "output20.txt");

    // vector<double> B = {10,0.662364,0.538213,0.48103,0.434498,0.396992,0.340488,0.268908,0.190814,0.1};
    // test(L, B);
    //vector<double> B = {10,0.662364,0.538213,0.48103,0.434498,0.396992,0.340488,0.268908,0.190814,0.1};
    vector<double> B(10, 0);
    B[0] = 10;
    B[B.size() - 1] = 0.1;
    double range = B[0] - B[B.size()-1];
    double step = range / (B.size() - 1);
    for (int i = 1; i < B.size()-1; ++i) {
        B[i] = B[i - 1] - step;
    }
    //optimize(L, B, "output20.txt");

    // test(L, B);

    original_optimize(L, B, "");

    // stop the timer
    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "\nexecution time: " << duration.count()/1e6 << " s\n";
}

void original_optimize(int L, vector<double> B, const char* filename) {
    cout << "starting original optimization algorithm with B = " << B << "and L = " << L << endl;

    unsigned int seed = 0;
    RNGEngine main_eng(seed);

    auto best_reps = PT(L, B, 1);

    cout << "burn in" << endl;
    auto p = best_reps.Start(1e5);

    auto best_tau = PT::Expected_RT(p);
    cout << "the initial round trip time is " << best_tau << endl;

    auto new_reps = best_reps;
#pragma omp parallel for firstprivate(new_reps)
    for (int i = 0; i < 20000; ++i) {
#pragma omp critical
        {
            new_reps = best_reps;
            //cout << "attempting an adjustment " << i << endl;
            new_reps.Adjustment(main_eng);
            //cout << "the new temperature set is " << new_reps.GetBetas() << endl;
        }

        auto p_new = new_reps.Start(1e5);
        auto tau_new = PT::Expected_RT(p_new);

#pragma omp critical
        {
            //cout << "the new round trip time is " << tau_new << endl;
            if (tau_new < best_tau || best_tau == INFINITY) {
                best_reps = new_reps;
                best_tau = tau_new;
                cout << "the best round trip time is now " << best_tau << endl;
            }
            else {
                //cout << "the best round trip time is still " << best_tau << endl;
            }
        }
    }

    cout << "the best round trip time was " << best_tau << endl;
    cout << "the best temperature set was " << best_reps.GetBetas() << endl;
}

void optimize(int L, vector<double> B, const char* filename) {
    cout << "starting optimization with B = " << B << "and L = " << L << endl;

    unsigned int seed = 0;
    RNGEngine main_eng(seed);          // SHARED, ALWAYS USED FOR TEMP SET UPDATES

    auto best_reps = PT(L, B, 1);  // PT OBJECTS CONTAIN THEIR OWN RNG ENGINE
                                   // COPY CONSTR. COPIES THE CURRENT STATE OF THE ENGINE

    auto best_tau = PT::Expected_RT(best_reps.Start(1e5));

    vector<double> rt_list;
    vector<vector<double>> beta_list;
    vector<vector<double>> p_list;

    // MAXIMUM NO. OF INSERTIONS AND DELETIONS
    for (int i = 0; i < 30; ++i) {
        auto best_reps_n = best_reps;
        if (best_reps_n.size() < 3 || best_tau == INFINITY) best_reps_n.Insertion(main_eng);
        else {
            if (RNG::zero_one_int(main_eng)) {
                best_reps_n.Insertion(main_eng);
            }
            else {
                best_reps_n.Deletion(main_eng);
            }
        }

        // RECORD BEST STATS FOR THIS  N ONLY
        auto best_tau_n = INFINITY;
        auto local_reps = best_reps_n;
        cout << local_reps.GetBetas();

        // DOUBLE INTERVAL WITH EACH ITERATION OF OUTER LOOP
        int step = 100;
        for (int k = 0; k < 10; ++k) {
            // BEFORE EACH CHECK, RECALCULATE EXPECTED VALUE OF THIS N
            best_tau_n = PT::Expected_RT(best_reps_n.Start(1e5));
            bool changed = false;

            // MINIMUM INTERVAL PLUS STEP THAT DOUBLES EACH TIME
            int check = 1000 + step;
            step *= 2;

            // COPY OF local_reps IS INITIALIZED FOR EACH THREAD
            #pragma omp parallel for firstprivate(local_reps)
            for (int j = 0; j < check; ++j) {
                #pragma omp critical 
                {   
                    // CRITICAL: READ CURRENT BEST SET
                    local_reps = best_reps_n;
                    local_reps.Adjustment(main_eng);
                }
                // PARALLEL SECTION
                auto p_new = local_reps.Start(1e8);
                auto tau_new = PT::Expected_RT(p_new);
                #pragma omp critical
                {
                    // CRITICAL: UPDATE CURRENT BEST TEMP SET, IF BETTER
                    if (tau_new <= best_tau_n) {
                        changed = true;
                        cout << best_tau_n << " --- " << local_reps.GetBetas();
                        best_tau_n = tau_new;
                        best_reps_n = local_reps;

                        rt_list.push_back(best_tau_n);
                        beta_list.push_back(best_reps_n.GetBetas());
                        p_list.push_back(p_new);
                    }
                }
            }
            // IF A STEP/ITERATION FAILS TO LOWER BEST TAU FOR THIS N, STOP
            if (!changed || best_tau_n == INFINITY) break;
        }
        // UPDATE BEST TAU BASED ON RESULTS FROM THIS N
        // (OUTER LOOP NEVER STOPS EARLY)
        if (best_tau_n <= best_tau || best_tau == INFINITY) {
            best_tau = best_tau_n;
            best_reps = best_reps_n;
        }
        if (best_tau_n == INFINITY) {
            i--;
        }
    }

    cout << best_reps.GetBetas();
    ofstream ofile(filename);
    ofile << seed << endl;
    for (unsigned int i = 0; i < rt_list.size(); ++i) {
        ofile << rt_list[i] << " --- " << beta_list[i];
    }
    ofile << best_reps.GetBetas();
    ofile.close();
}

void test(int L, vector<double> B) {
    cout << "starting a test with B = " << B << "and L = " << L << endl;

    unsigned int seed = 0;
    RNGEngine main_eng(seed);          // SHARED, ALWAYS USED FOR TEMP SET UPDATES
    unsigned int nupdates = 1;

    auto best_reps = PT(L, B, 1);  // PT OBJECTS CONTAIN THEIR OWN RNG ENGINE
    // COPY CONSTR. COPIES THE CURRENT STATE OF THE ENGINE

    auto local_reps = best_reps;

    // COPY OF local_reps IS INITIALIZED FOR EACH THREAD
#pragma omp parallel for firstprivate(local_reps) num_threads(2)
    for (int j = 0; j < nupdates; ++j) {
#pragma omp critical 
        {
            // CRITICAL: READ CURRENT BEST SET
            local_reps = best_reps;
            local_reps.Adjustment(main_eng);
        }
        // PARALLEL SECTION
        auto p_new = local_reps.Start(1e8);
        auto tau_new = PT::Expected_RT(p_new);
#pragma omp critical
        {
            // CRITICAL: 
            cout << "finished run " << j << endl;
        }
    }
}