#pragma once

#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <climits>
#include <chrono>
#include <fstream>
#include <float.h>
#include <omp.h>


std::ostream& operator<<(std::ostream& os, const std::vector<double>& v);
std::ostream& operator<<(std::ostream& os, const std::vector<int>& v);
std::vector<double> divide(const std::vector<double>& v, double d);

double ave_thermal_energy(int L, double B);

int cost_helper(const std::vector<std::vector<int>>& lat);

void helper(std::vector<std::vector<int>>& lat, double B, double& num, double& den, unsigned int i, unsigned int j);

void write_to_file(const std::vector<double>& T, const std::vector<double>& p_ave, std::string);