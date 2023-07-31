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


void optimize(int L, vector<double> B, const char* filename);

void parse(const char* filename, vector<int>& L_list, vector<vector<double>>& B_list);

int main(int argc, char** argv) {
    vector<double> T = {1./10, 1/.1};

    // convert to inverse temperature
    vector<double> B;
    for (auto& t: T) {
        B.push_back(pow(t, -1));
    }

    // start the timer
    auto start = high_resolution_clock::now();

    // optimize(5, B, "output1.txt");
    // optimize(10, B, "output2.txt");
    optimize(15, B, "output15.txt");

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "\nexecution time: " << duration.count()/1e6 << " s\n";
}


void optimize(int L, vector<double> B, const char* filename) {
    cout << "starting optimization with B = " << B << "and L = " << L << endl;

    unsigned int seed = time(0);
    pcg32 main_eng(seed);

    auto best_reps = PT(L, B, 1);
    auto best_tau = PT::Expected_RT(best_reps.Start(1e5));

    vector<double> rt_list;
    vector<vector<double>> beta_list;
    vector<vector<double>> p_list;

    for (int i = 0; i < 25; ++i) {
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
                        cout << best_tau_n << " --- " << local_reps.GetBetas();
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
    ofstream ofile(filename);
    ofile << seed << endl;
    for (unsigned int i = 0; i < rt_list.size(); ++i) {
        ofile << rt_list[i] << " --- " << beta_list[i];
    }
    ofile << best_reps.GetBetas();
    ofile.close();
}

void parse(const char* filename, vector<int>& L_list, vector<vector<double>>& B_list) {
    ofstream ifile(filename);
    
}
