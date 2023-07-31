#include "util.h"

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

void write_to_file(const std::vector<double>& T, const std::vector<double>& p_ave, std::string filename) {
    std::ofstream ofile(filename);
    ofile << "T,swap_prob";
    for (unsigned int i = 0; i < p_ave.size(); ++i) {
        ofile << std::endl << T[i] << "," << p_ave[i];
    }
    ofile.close();
}