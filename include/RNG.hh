#pragma once
#include <random>
#include "pcg/pcg_random.hpp"

class RNG {
    public:
    static std::uniform_int_distribution<int> zero_one_int;
    static std::uniform_real_distribution<double> zero_one_double;
    static std::uniform_int_distribution<int> uniform_int;
    static std::uniform_int_distribution<int> pick_op;
};