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
    std::vector<unsigned int> count_acc(reps.size()-1, 0.);
    std::vector<double> old_p(reps.size()-1, 0.);
    std::vector<double> p(reps.size()-1, 0.);

    int step_size = 100;

    std::vector<double> mean_i(reps.size()-1, 0.);
    for (int i = 1; i < n_sweeps+1; ++i) {
        for (auto& rep: reps) {
            rep->Update();
        }

        for (unsigned int j = 0; j < reps.size()-1; ++j) {
            double dB = reps[j+1]->B - reps[j]->B;
            double dE = reps[j+1]->cost - reps[j]->cost;
            if (dB*dE > 0. || RNG::zero_one_double(gen) < exp(dB*dE)) {
                std::swap(reps[j], reps[j+1]);
                std::swap(reps[j]->B, reps[j+1]->B);
                ++count_acc[j];
            }
        }

        if (!(i%step_size)) {
            bool terminate = true;
            for (unsigned int k = 0; k < p.size(); ++k) {
                p[k] = (double)count_acc[k]/i;
                if (abs(p[k]-old_p[k])/p[k] > .001) terminate = false;
                old_p[k] = p[k];
            }
            if (terminate) {
                return p;
            }
        }
    }

    for (unsigned int k = 0; k < p.size(); ++k) {
        p[k] = (double)count_acc[k]/n_sweeps;
    }
    std::cout << "reached max MC iterations\n";
    return p;
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
    