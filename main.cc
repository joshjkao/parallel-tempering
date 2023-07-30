#include "Replica.h"
#include "PT.h"
#include "util.h"
#include "RNG.h"
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
    vector<double> T = {1./10, 1./2, 1/.1};
    // 0.2,0.320761,0.418838,0.525278,10
    int L = 10;

    cout << "starting optimization with T = " << T << "and L = " << L << endl;

    // convert to inverse temperature
    vector<double> B;
    for (auto& t: T) {
        B.push_back(pow(t, -1));
    }

    // start the timer
    auto start = high_resolution_clock::now();

    pcg32 main_eng(time(0));

    auto best_reps = PT(L, B, 1);
    auto best_tau = PT::Expected_RT(best_reps.Start(1e5));

    vector<double> rt_list;
    vector<vector<double>> beta_list;
    vector<vector<double>> p_list;

    for (int i = 0; i < 3; ++i) {
        auto best_reps_n = best_reps;
        best_reps_n.Insertion(main_eng);
        auto best_tau_n = INFINITY;
        auto local_reps = best_reps_n;
        for (int k = 0; k < 10; ++k) {
            best_tau_n = PT::Expected_RT(best_reps_n.Start(1e5));
            #pragma omp parallel for firstprivate(local_reps)
            for (int j = 0; j < 1000; ++j) {
                #pragma omp critical 
                {
                    local_reps = best_reps_n;
                    local_reps.Adjustment(main_eng);
                }
                auto p_new = local_reps.Start(1e7);
                auto tau_new = PT::Expected_RT(p_new);
                #pragma omp critical
                {
                    if (tau_new < best_tau_n) {
                        cout << best_tau_n << " --- " << local_reps.GetBetas() << endl;
                        best_tau_n = tau_new;
                        best_reps_n = local_reps;

                        rt_list.push_back(best_tau_n);
                        beta_list.push_back(best_reps_n.GetBetas());
                        p_list.push_back(p_new);
                    }
                }
            }
        }
        if (best_tau_n <= best_tau) {
            best_tau = best_tau_n;
            best_reps = best_reps_n;
        }
        else {
            break;
        }
    }

    cout << best_reps.GetBetas();


    ofstream ofile("output.csv");
    for (unsigned int i = 0; i < rt_list.size(); ++i) {
        ofile << rt_list[i] << " --- " << beta_list[i] << endl;
    }
    ofile.close();



    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "\nexecution time: " << duration.count()/1e6 << " s\n";
}