#pragma once
#include "Replica.h"


class PT {
    public:
    PT(int L, const std::vector<double>& B, unsigned int seed);
    PT(const PT& other);
    PT& operator=(const PT& rhs);
    ~PT();
    void Seed(unsigned int seed);

    void Adjustment(pcg32& main_eng);
    void Insertion(pcg32& main_eng);
    void Deletion(pcg32& main_eng);

    std::vector<double> Start(unsigned int n_sweeps);
    static double Expected_RT(const std::vector<double>& p);
    std::vector<double> GetBetas();

    std::vector<IsingFerromagnetReplica*> reps;
    pcg32 eng;
};