#pragma once
#include "Replica.hh"


class PT {
    public:

    static std::vector<Replica*> construct_replicas(int L, const std::vector<double>& B, pcg32& gen);
    static void delete_replicas(std::vector<Replica*>& reps);

    static std::vector<double> start(int n_sweeps, std::vector<Replica*>& reps, pcg32& gen);
    static double expected_rt(std::vector<double> p);

    static std::vector<double> verify_ave_energy(int n_sweeps, std::vector<Replica*>& reps);

    static void internal_adj(std::vector<double>& B, pcg32& gen);
    static void insertion(std::vector<double>& B, pcg32& gen);
    static void deletion(std::vector<double>& B, pcg32& gen);
    
};