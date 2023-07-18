#include "RNG.hh"
#include <random>

std::uniform_int_distribution<int> RNG::zero_one_int(0,1);
std::uniform_real_distribution<double> RNG::zero_one_double(0,1);
std::uniform_int_distribution<int> RNG::uniform_int(0, RAND_MAX);
std::uniform_int_distribution<int> RNG::pick_op(0, 99);