#pragma once
#include <vector>
#include <random>
#include "pcg/pcg_random.hpp"

class Replica {
    public:
        Replica() = default;
        virtual ~Replica() = default;

        virtual void Init() = 0;
        virtual int Cost() = 0;
        virtual void Update() = 0;

        virtual void Dump() {}

        int cost;
        double B;
};

class IsingFerromagnetReplica: public Replica {
    public:
        IsingFerromagnetReplica(int L_, double B_, pcg32* gen_);
        ~IsingFerromagnetReplica();

        void Init();
        int Cost();
        void Update();

        void Dump();

        unsigned int L;
        std::vector<std::vector<int>> lattice;

    private:
        inline int arithmod(int i);
        inline double int_exp(int x);
        std::vector<int> pick_order;
        std::vector<int> Ldecr;
        std::vector<int> Lincr;

        pcg32* gen;
};

int IsingFerromagnetReplica::arithmod(int i) {
    return (i + L) % L;
}