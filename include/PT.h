#pragma once
#include "Replica.h"

typedef IsingFerromagnetReplica ReplicaType;

class PT {
    public:
        static long double Expected_RT(const std::vector<double>& p);

        PT(int L, const std::vector<double>& B, unsigned int seed);

        void Seed(unsigned int seed);

        void Adjustment(RNGEngine& main_eng);
        void Insertion(RNGEngine& main_eng);
        void Deletion(RNGEngine& main_eng);

        std::vector<double> Start(unsigned int n_sweeps);
        std::vector<double> GetBetas();
        inline size_t size() { return reps.size(); }
        std::vector<unsigned int> up() { return h_up; }
        std::vector<unsigned int> down() { return h_down; }

        std::vector<ReplicaType> reps;
        std::vector<unsigned int> h_up;
        std::vector<unsigned int> h_down;
        RNGEngine eng;
};