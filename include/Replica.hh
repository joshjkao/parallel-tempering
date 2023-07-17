#pragma once
#include <vector>

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
        IsingFerromagnetReplica(int L_, double B_);
        ~IsingFerromagnetReplica();

        void Init();
        int Cost();
        void Update();

        void Dump();

        unsigned int L;
        std::vector<std::vector<int>> lattice;    

    private:
        inline int arithmod(int i);
        std::vector<int> pick_order;
        std::vector<int> Ldecr;
        std::vector<int> Lincr;
};

int IsingFerromagnetReplica::arithmod(int i) {
    return (i + L) % L;
}