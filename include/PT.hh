#pragma once
#include "Replica.hh"


class PT {
    public:
    static std::vector<double> start(int n_sweeps, std::vector<Replica*>& reps);
    static std::vector<double> verify_ave_energy(int n_sweeps, std::vector<Replica*>& reps);
    static double expected_rt(std::vector<double> p);
};