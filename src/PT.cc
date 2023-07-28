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


PT::PT(int L, const std::vector<double>& B, unsigned int seed):
eng(seed)
{
    for (auto& b: B) {
        reps.push_back(new IsingFerromagnetReplica(L, b));
    }
    for (auto& rep: reps) {
        rep->Init(eng);
    }
}

PT::PT(const PT& other) {
    eng = pcg32(other.eng);
    for (auto& orep: other.reps) {
        reps.push_back(new IsingFerromagnetReplica(*orep));
    }
}

PT& PT::operator=(const PT& rhs) {
    eng = rhs.eng;
    for (auto& r: reps) {
        delete r;
    }
    reps.clear();
    for (auto& orep: rhs.reps) {
        reps.push_back(new IsingFerromagnetReplica(*orep));
    }
    return *this;
}

PT::~PT() {
    for (auto& rep: reps) {
        delete rep;
    }
}


std::vector<double> PT::Start(unsigned int n_sweeps) {
    std::vector<unsigned int> count_acc(reps.size()-1, 0.);
    std::vector<double> old_p(reps.size()-1, 0.);
    std::vector<double> p(reps.size()-1, 0.);

    int step_size = 100;

    std::vector<double> mean_i(reps.size()-1, 0.);
    for (unsigned int i = 1; i < n_sweeps+1; ++i) {
        for (auto& rep: reps) {
            rep->Update(eng);
        }

        for (unsigned int j = 0; j < reps.size()-1; ++j) {
            double dB = reps[j+1]->B - reps[j]->B;
            double dE = reps[j+1]->cost - reps[j]->cost;
            if (dB*dE > 0. || RNG::zero_one_double(eng) < exp(dB*dE)) {
                std::swap(reps[j], reps[j+1]);
                std::swap(reps[j]->B, reps[j+1]->B);
                ++count_acc[j];
            }
        }

        // if (!(i%step_size)) {
        //     bool terminate = true;
        //     for (unsigned int k = 0; k < p.size(); ++k) {
        //         p[k] = (double)count_acc[k]/i;
        //         if (abs(p[k]-old_p[k])/p[k] > .0005) terminate = false;
        //         old_p[k] = p[k];
        //     }
        //     if (terminate) {
        //         return p;
        //     }
        // }
    }

    for (unsigned int k = 0; k < p.size(); ++k) {
        p[k] = (double)count_acc[k]/n_sweeps;
    }
    std::cout << "reached max MC iterations\n";
    return p;
}


double PT::Expected_RT() {
    return 0.;
}



void PT::Adjustment() {
    return;
}


void PT::Insertion() {
    return;
}


void PT::Deletion() {
    return;
}
    