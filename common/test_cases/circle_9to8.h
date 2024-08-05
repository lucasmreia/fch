///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#pragma once

#include "test_case_template.h"

template <>
struct TestCaseImpl<TestCase::CIRCLE_9TO8> {
    static constexpr size_t NDIM = 9;
    static constexpr size_t KDIM = 8;

    static constexpr std::array<double, NDIM> DOMAIN_MIN{
        -1.1,
        -1.1,
        -1.1,
        -1.1,
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
        1.1,
        1.1,
        1.1,
        1.1,
        1.1
    };
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{
        11,
        11,
        11,
        11,
        11,
        11,
        11,
        11,
        11,
    };

    static constexpr std::array<double, NDIM> FIRST_POINT{
        1,
        1.000099995000500,
        1.000149988751687,
        1.000199980003999,
        1.000249968757810,
        1.000299955013495,
        1.000349938771428,
        1.000399920031984,
        0.01
    };

    static void func(const double *x, double *f) {
        f[0] = x[0] * x[0] + x[8] * x[8] - 1.0001;
        f[1] = x[1] * x[1] - 1.0002;
        f[2] = x[2] * x[2] - 1.0003;
        f[3] = x[3] * x[3] - 1.0004;
        f[4] = x[4] * x[4] - 1.0005;
        f[5] = x[5] * x[5] - 1.0006;
        f[6] = x[6] * x[6] - 1.0007;
        f[7] = x[7] * x[7] - 1.0008;
    }
};
