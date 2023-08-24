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
    // vector<double> T = {1./10, 1/.1};

    // // convert to inverse temperature
    // vector<double> B;
    // for (auto& t: T) {
    //     B.push_back(pow(t, -1));
    // }

    // vector<double> B = {10,0.619754,0.50124,0.435079,0.386917,0.327503,0.25645,0.182552,0.1};
    vector<double> B = {10, .1};

    // start the timer
    auto start = high_resolution_clock::now();

    // optimize(3, B, "output3.txt");
    optimize(5, B, "output5.txt");
    // optimize(10, B, "output10.txt");
    // optimize(13, B, "output13.txt");
    // optimize(15, B, "output15.txt");
    // optimize(20, B, "output20.txt");

    auto end = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(end - start);
    cout << "\nexecution time: " << duration.count()/1e6 << " s\n";
}


void optimize(int L, vector<double> B, const char* filename) {
    cout << "starting optimization with B = " << B << "and L = " << L << endl;

    unsigned int seed = 0;
    pcg32 main_eng(seed);

    auto best_reps = PT(L, B, 1);
    auto best_tau = PT::Expected_RT(best_reps.Start(1e5));

    vector<double> rt_list;
    vector<vector<double>> beta_list;
    vector<vector<double>> p_list;

    for (int i = 0; i < 20; ++i) {
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
        auto best_tau_n = INFINITY;
        auto local_reps = best_reps_n;
        cout << local_reps.GetBetas();
        double mavg = -1;
        int step = 100;
        for (int k = 0; k < 10; ++k) {
            best_tau_n = PT::Expected_RT(best_reps_n.Start(1e5));
            int check = 1000 + step;
            step *= 2;
            #pragma omp parallel for firstprivate(local_reps)
            for (int j = 0; j < check; ++j) {
                #pragma omp critical 
                {
                    local_reps = best_reps_n;
                    local_reps.Adjustment(main_eng);
                }
                auto p_new = local_reps.Start(1e7);
                auto tau_new = PT::Expected_RT(p_new);
                #pragma omp critical
                {
                    if (tau_new <= best_tau_n) {
                        cout << best_tau_n << " --- " << local_reps.GetBetas();
                        best_tau_n = tau_new;
                        best_reps_n = local_reps;

                        rt_list.push_back(best_tau_n);
                        beta_list.push_back(best_reps_n.GetBetas());
                        p_list.push_back(p_new);
                    }
                }
            }
            if (mavg == -1) mavg = best_tau_n;
            else {
                double old_mavg = mavg;
                mavg = (1-.2)*mavg + .2*best_tau_n;
                if (abs(mavg-old_mavg)/mavg < .05) break;
                else if (mavg == INFINITY && old_mavg == INFINITY) break;
            }
        }
        if (best_tau_n <= best_tau || best_tau == INFINITY) {
            best_tau = best_tau_n;
            best_reps = best_reps_n;
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
