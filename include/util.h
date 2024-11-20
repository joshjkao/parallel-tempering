#pragma once

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <climits>
#include <chrono>
#include <fstream>
#include <float.h>


std::vector<double> divide(const std::vector<double>& v, double d);

void write_to_file(const std::vector<double>& T, const std::vector<double>& p_ave, std::string);

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& v) {
    for (auto& x: v) {
        os << x << ",";
    }
    os << std::endl;
    return os;
}