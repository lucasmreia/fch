///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#pragma once

#include "test_case_template.h"

template <>
struct TestCaseImpl<TestCase::HYPERSHPERE_5TO1> {
    static constexpr size_t NDIM = 5;
    static constexpr size_t KDIM = 1;

    static constexpr std::array<double, NDIM> DOMAIN_MIN{
        -1.1,
        -1.1,
        -1.1,
        -1.1,
        -1.1
    };
    static constexpr std::array<double, NDIM> DOMAIN_MAX{
        1.1,
        1.1,
        1.1,
        1.1,
        1.1
    };
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{
        20,
        20,
        20,
        20,
        20,
    };

    static constexpr std::array<double, NDIM> FIRST_POINT{
        1,
        0,
        0,
        0,
        0
    };

    static void func(const double *x, double *f) {
        f[0] = x[0] * x[0] +
               x[1] * x[1] +
               x[2] * x[2] +
               x[3] * x[3] +
               x[4] * x[4] - 1.0;
    }
};
