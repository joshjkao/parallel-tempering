#pragma once
#include "Replica.hh"


class PT {
    public:
    PT(int L, const std::vector<double>& B, unsigned int seed);
    PT(const PT& other);
    PT& operator=(const PT& rhs);
    ~PT();

    void Adjustment();
    void Insertion();
    void Deletion();

    std::vector<double> Start(unsigned int n_sweeps);
    double Expected_RT();

    std::vector<IsingFerromagnetReplica*> reps;
    pcg32 eng;
    
};