#include "Replica.hh"
#include <random>
#include <iostream>
#include <algorithm>
#include <cmath>
#include <climits>


IsingFerromagnetReplica::IsingFerromagnetReplica(int L_, double B_):
Replica(), L(L_)
{
    lattice = std::vector< std::vector<int> >(L_, std::vector<int>(L_, 0));
    pick_order = std::vector<int>(L_*L_, 0);
    Ldecr = std::vector<int>(L, 0);
    Lincr = std::vector<int>(L, 0);
    cost = 0;
    Init();
    B = B_;
}

IsingFerromagnetReplica::~IsingFerromagnetReplica() {

}

void IsingFerromagnetReplica::Init() {
    for (unsigned int i = 0; i < L; ++i) {
        Lincr[i] = (i+1)%L;
        Ldecr[i] = arithmod(i-1);
    }

    unsigned int i = 0;
    for (auto& row: lattice) {
        for (auto& s: row) {
            s = rand() % 2;
            if (s == 0) {
                s = -1;
            }
            pick_order[i] = i;
            ++i;
        }
    }

    Cost();
}

int IsingFerromagnetReplica::Cost() {
    int sum = 0;
    for (unsigned int i = 0; i < L; ++i) {
        for (unsigned int j = 0; j < L; ++j) {
            sum += (lattice)[i][j] * ( (lattice)[i][Lincr[j]] + (lattice)[Lincr[i]][j] );
        }
    }
    cost = -sum;
    return cost;
}

void IsingFerromagnetReplica::Update() {
// Perform sweep to evolve replica
    std::random_shuffle(pick_order.begin(), pick_order.end());
    for (auto& x: pick_order) {
        int i = x/L;
        int j = x%L;
        int contr = (lattice)[i][j] * ((lattice)[i][Lincr[j]] + (lattice)[i][Ldecr[j]] + (lattice)[Ldecr[i]][j] + (lattice)[Lincr[i]][j]);
        double A = std::min(1., exp(-B * (2*contr)));
        if ((double)rand()/(double)INT_MAX < A) {
            cost += 2*contr;
            (lattice)[i][j] *= -1;
        }
    }
}

void IsingFerromagnetReplica::Dump() {
    for (auto& row: (lattice)) {
        for (auto& s: row) {
            std::cout << s << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl << "cost: " << cost << std::endl;
}
