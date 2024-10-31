#pragma once
#include "Replica.h"


class PT {
    public:
        static long double Expected_RT(const std::vector<double>& p);

        PT(int L, const std::vector<double>& B, unsigned int seed);
        PT(const PT& other);
        PT& operator=(const PT& rhs);
        ~PT();

        void Seed(unsigned int seed);

        void Adjustment(pcg32& main_eng);
        void Insertion(pcg32& main_eng);
        void Deletion(pcg32& main_eng);

        std::vector<double> Start(unsigned int n_sweeps);
        std::vector<double> GetBetas();
        inline size_t size() {return reps.size();}

        std::vector<IsingFerromagnetReplica*> reps;
        pcg32 eng;
};