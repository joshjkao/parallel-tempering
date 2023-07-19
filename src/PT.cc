#include "PT.hh"
#include "util.hh"
#include "RNG.hh"
#include "pcg/pcg_random.hpp"

#include <cmath>
#include <climits>
#include <random>
#include <iostream>
#include <algorithm>
#include <float.h>


std::vector<Replica*> PT::construct_replicas(int L, const std::vector<double>& B, pcg32& gen) {
    std::vector<Replica*> ret;
    for (auto& b: B) {
        ret.push_back(new IsingFerromagnetReplica(L, b, &gen));
    }
    return ret;
}


void PT::delete_replicas(std::vector<Replica*>& reps) {
    for (auto& r: reps) {
        delete r;
    }
    return;
}


std::vector<double> PT::start(int n_sweeps, std::vector<Replica*>& reps, pcg32& gen) {
    std::vector<double> p_acc(reps.size()-1, 0.);
    std::vector<double> variance_acc(reps.size()-1, 0.);

    int step_size = 1e3;

    for (int i = 0; i < n_sweeps/step_size; ++i) {
        // std::cout << "Event " << i*step_size << " starts\n"; 
        std::vector<double> mean_i(reps.size()-1, 0.);
        for (int j = 0; j < step_size; ++j) {
            // #pragma omp parallel for schedule(static)
            for (auto& rep: reps) {
                rep->Update();
            }

            for (unsigned int j = 0; j < reps.size()-1; ++j) {
                double dB = reps[j+1]->B - reps[j]->B;
                double dE = reps[j+1]->cost - reps[j]->cost;
                // double A = std::min(1., exp(dB*dE));
                if (dB*dE > 0. || RNG::zero_one_double(gen) < exp(dB*dE)) {
                    std::swap(reps[j], reps[j+1]);
                    std::swap(reps[j]->B, reps[j+1]->B);
                    ++p_acc[j];
                    ++mean_i[j];
                }
            }
        }
        std::vector<double> mean_total;
        for (unsigned int j = 0; j < reps.size()-1; ++j) {
            mean_total.push_back(p_acc[j]/((i+1)*step_size));
            mean_i[j] /= step_size;
            variance_acc[j] += pow(mean_total[j]-mean_i[j], 2);
        }

    }

    for (auto& p: p_acc) {
        p/=n_sweeps;
    }
    // std::cout << "WARNING: reached maxed MC iterations" << std::endl;
    return p_acc;
}


std::vector<double> PT::verify_ave_energy(int n_sweeps, std::vector<Replica*>& reps) {
    std::vector<int> sum_f(reps.size(), 0.);

    int print_progress = 1e4;

    for (int i = 0; i < n_sweeps; ++i) {

        if (i%print_progress == 0) std::cout << "Event " << i << " starts\n"; 

        #pragma omp parallel for
        for (unsigned int j = 0; j < reps.size(); ++j) {
            reps[j]->Update();
            sum_f[j] += reps[j]->cost;
        }

    }

    std::vector<double> avg_energy(reps.size(), 0.);
    for (unsigned int j = 0; j < sum_f.size(); ++j) {
        avg_energy[j] = (double)sum_f[j]/n_sweeps;
    }

    return avg_energy;
}


double PT::expected_rt(std::vector<double> p) {
    int n = p.size()+1;
    double sum = 0;
    for (auto& s: p) {
        sum += 1/s;
    }
    return n*(n-1)*sum;
}

void PT::internal_adj(std::vector<double>& B, pcg32& gen) {
    int adj = RNG::uniform_int(gen)%(B.size()-2)+1;
    double Bl = B[0];
    double Br = B[B.size()-1];

    B[adj] = RNG::zero_one_double(gen) * (Br-Bl) + Bl;
    sort(B.begin(), B.end());
}


void PT::insertion(std::vector<double>& B, pcg32& gen) {
    if (B.size() < 3) return;
    int ins = RNG::uniform_int(gen)%(B.size()-1)+1;
    double Bl = B[ins-1];
    double Br = B[ins];
    double Bnew = RNG::zero_one_double(gen) * (Br-Bl) + Bl;
    B.insert(B.begin()+ins, Bnew);
}


void PT::deletion(std::vector<double>& B, pcg32& gen) {
    int rem = RNG::uniform_int(gen)%(B.size()-2)+1;
    B.erase(B.begin()+rem);
}
    