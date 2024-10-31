#pragma once
#include <vector>
#include <random>
#include "RNG.h"


class Replica {
    public:
        Replica() = default;
        Replica(const Replica& other) = default;
        virtual ~Replica() = default;

        virtual void Init(RNGEngine& eng) = 0;
        virtual int Cost() = 0;
        virtual void Update(RNGEngine& eng) = 0;

        virtual void Dump() {}

        int cost;
        double B;
};


class IsingFerromagnetReplica: public Replica {
    public:
        IsingFerromagnetReplica(int L_, double B_);
        ~IsingFerromagnetReplica();

        void Init(RNGEngine& eng);
        int Cost();
        void Update(RNGEngine& eng);

        void Dump();

        unsigned int L;
        std::vector<std::vector<int>> lattice;

    private:
        inline int arithmod(int i);
        inline double int_exp(int x);
        std::vector<int> pick_order;
        std::vector<int> Ldecr;
        std::vector<int> Lincr;
        
};


int IsingFerromagnetReplica::arithmod(int i) {
    return (i + L) % L;
}