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
    int L = 16;
    
    if (argc > 1) {
	for (int i = 0; i < argc; ++i) {
            if (string(argv[i]) == "-m") { 
		cout << argv[i+1] << endl;
	    } else if (string(argv[i]) == "-L") {
		L = atoi(argv[i+1]);
	    }
	}
    }

    // start the timer
    auto start = high_resolution_clock::now();

    // perform optimization
    //vector<double> B = {10, .1};
    //optimize(L, B, "output20.txt");

    // vector<double> B = {10,0.662364,0.538213,0.48103,0.434498,0.396992,0.340488,0.268908,0.190814,0.1};
    // test(L, B);
    //vector<double> B = {10,0.662364,0.538213,0.48103,0.434498,0.396992,0.340488,0.268908,0.190814,0.1};
    vector<double> B(20, 0);
    B[0] = 10;
    B[B.size() - 1] = 0.1;
    double range = B[0] - B[B.size()-1];
    double step = range / (B.size() - 1);
    for (int i = 1; i < B.size()-1; ++i) {
        B[i] = B[i - 1] - step;
    }
    //B = {10,9.47895,8.18324,8.04597,7.41457,7.36569,7.30659,6.35263,1.2814,0.713888,0.624601,0.551903,0.497851,0.468332,0.4278,0.383205,0.331821,0.268845,0.19389,0.1};
    //B = {10,9.57445,8.7345,8.01931,6.58886,0.934614,0.82281,0.724389,0.63988,0.574662,0.53832,0.500207,0.462507,0.434108,0.40361,0.366129,0.31381,0.25219,0.180653,0.1};
    //B = {10,7.05052,4.66507,3.27823,1.13212,0.791818,0.669847,0.608387,0.55638,0.523979,0.484614,0.461364,0.436778,0.417395,0.389406,0.349967,0.296762,0.231269,0.169787,0.1};
    //B = {10,4.67737,2.72574,0.847269,0.715726,0.635995,0.617965,0.564841,0.528914,0.501875,0.479639,0.453024,0.428162,0.402818,0.371865,0.332404,0.287669,0.226995,0.167552,0.1};
    //B = {10,5.76714,0.83096,0.801046,0.697532,0.650352,0.598702,0.54773,0.514966,0.480195,0.454756,0.429781,0.407565,0.381762,0.344282,0.306928,0.259876,0.211268,0.156177,0.1};
    //B = {10,1.32284,0.836715,0.702703,0.638139,0.599297,0.562789,0.525261,0.488333,0.455047,0.437187,0.415893,0.392294,0.365957,0.336062,0.301478,0.254297,0.204025,0.151595,0.1};
    //B = {10,0.991885,0.773716,0.685872,0.623665,0.572046,0.529449,0.49306,0.47219,0.448207,0.425505,0.408578,0.389556,0.363422,0.330791,0.293416,0.253985,0.204948,0.149088,0.1};
    //B = {10,0.891955,0.732779,0.677932,0.622229,0.565548,0.531036,0.494341,0.46424,0.439804,0.420873,0.400526,0.377164,0.351881,0.322231,0.283375,0.240621,0.194289,0.151904,0.1};
    //B = {10,0.803415,0.706419,0.656094,0.590807,0.551518,0.515002,0.485782,0.458216,0.434074,0.417434,0.395171,0.367197,0.341951,0.314868,0.283553,0.241135,0.190548,0.144848,0.1};
    //B = {10,0.803415,0.682838,0.617702,0.565577,0.534245,0.504565,0.479108,0.45916,0.437047,0.420797,0.399228,0.378461,0.351065,0.319009,0.277878,0.236072,0.190952,0.144976,0.1};
    //B = {10,0.803415,0.686698,0.615372,0.582873,0.530332,0.495937,0.467841,0.44434,0.425736,0.405442,0.387065,0.363037,0.336003,0.29853,0.266927,0.232783,0.187775,0.145409,0.1};
    //B = {10,0.803415,0.691273,0.616145,0.566909,0.525436,0.492405,0.464054,0.442107,0.422412,0.403158,0.385745,0.361657,0.334925,0.297267,0.264227,0.227735,0.182597,0.138412,0.1};
    //B = {10,0.803415,0.687169,0.618126,0.560575,0.522178,0.489235,0.461404,0.439098,0.421128,0.40301,0.380716,0.353594,0.324425,0.290869,0.258222,0.22392,0.183659,0.142824,0.1};
    //B = {10,0.808525,0.682604,0.607103,0.555901,0.516019,0.485856,0.46094,0.438792,0.421139,0.401899,0.379797,0.352124,0.322329,0.289901,0.255275,0.217891,0.18061,0.139694,0.1};
    B = {10,0.808525,0.690239,0.608825,0.55478,0.517458,0.483245,0.458075,0.436438,0.417221,0.396702,0.374917,0.348612,0.318465,0.286626,0.253612,0.216804,0.17905,0.140645,0.1};

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

    unsigned int seed = time(0);
    RNGEngine main_eng(seed);

    auto best_reps = PT(L, B, 1);

    cout << "burn in" << endl;
    auto p = best_reps.Start(1.6e7);

    auto best_tau = PT::Expected_RT(p);
    cout << "the initial round trip time is " << best_tau << endl;
    cout << "initial swap probabilities were " << p << endl;

    auto new_reps = best_reps;
    int redo = 1000;
    unsigned int counter = 0;
#pragma omp parallel for firstprivate(new_reps)
    for (int i = 0; i < 40000; ++i) {
#pragma omp critical
        {
	    counter++;
	    new_reps = best_reps;
            //cout << "attempting an adjustment " << i << endl;
	    if (counter % redo != 0) {            
		new_reps.Adjustment(main_eng);
	    } else {
		cout << "redo iteration, retest the current best" << endl;
	    }
            //cout << "the new temperature set is " << new_reps.GetBetas() << endl;
        }

        auto p_new = new_reps.Start(1.6e7);
        auto tau_new = PT::Expected_RT(p_new);

#pragma omp critical
        {
            //cout << "the new round trip time is " << tau_new << endl;
            if (tau_new < best_tau || best_tau == INFINITY) {
                best_reps = new_reps;
                best_tau = tau_new;
                cout << "the best round trip time is now " << best_tau << endl;
		cout << best_reps.GetBetas() << endl;	    
            }
            else if (counter % redo == 0) {
                //cout << "the best round trip time is still " << best_tau << endl;
                best_tau = tau_new;
	        }
            if (counter % redo == 0) {
		        cout << "n_up: " << new_reps.up();
                cout << "n_dwn: " << new_reps.down();
                cout << p_new << endl;
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

    auto best_tau = PT::Expected_RT(best_reps.Start(1.6e7));

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
