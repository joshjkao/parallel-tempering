#include "util.hh"

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <climits>
#include <chrono>
#include <fstream>
#include <float.h>
#include <omp.h>
#include <string>


std::ostream& operator<<(std::ostream& os, const std::vector<double>& v) {
    for (auto& x: v) {
        os << x << ",";
    }
    os << std::endl;
    return os;
}

std::ostream& operator<<(std::ostream& os, const std::vector<int>& v) {
    for (auto& x: v) {
        os << x << ",";
    }
    os << std::endl;
    return os;
}

std::vector<double> divide(const std::vector<double>& v, double d) {
    std::vector<double> ret;
    for (auto& i: v) {
        ret.push_back(i/d);
    }
    return ret;
}


int cost_helper(const std::vector<std::vector<int>>& lat) {
    int sum = 0;
    for (unsigned int i = 0; i < lat.size(); ++i) {
        for (unsigned int j = 0; j < lat.size(); ++j) {
            sum += lat[i][j] * ( lat[i][(j+1)%lat.size()] + lat[(i+1)%lat.size()][j] );
        }
    }
    return -sum;
}
void helper(std::vector<std::vector<int>>& lat, double B, double& num, double& den, unsigned int i, unsigned int j) {
    if (i == lat.size()) {
        double cost = (double)cost_helper(lat);
        // cout << cost_helper(lat) << endl;
        num += cost * exp(-B * cost);
        // cout << num << " ";
        den += exp(-B * cost);
        // cout << den << endl;
        return;
    }
    else {
        lat[i][j] = -1;
        if (j < lat.size()-1) helper(lat, B, num, den, i, j+1);
        else helper(lat, B, num, den, i+1, 0);
        lat[i][j] = 1;
        if (j < lat.size()-1) helper(lat, B, num, den, i, j+1);
        else helper(lat, B, num, den, i+1, 0);
        return;
    }
}
double ave_thermal_energy(int L, double B) {
    // ( sum_c f(c) exp(-beta f(c)))/( sum_c exp(-beta f(c)))

    std::vector<std::vector<int>> lat(L, std::vector<int>(L, -1));
    double num = 0;
    double den = 0;

    helper(lat, B, num, den, 0, 0);
    return num/den;
}

void write_to_file(const std::vector<double>& T, const std::vector<double>& p_ave, std::string filename) {
    std::ofstream ofile(filename);
    ofile << "T,swap_prob";
    for (unsigned int i = 0; i < p_ave.size(); ++i) {
        ofile << std::endl << T[i] << "," << p_ave[i];
    }
    ofile.close();
}