///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#pragma once

#include "test_case_template.h"

template <>
struct TestCaseImpl<TestCase::KLEIN_5TO3> {
    static constexpr size_t NDIM = 5;
    static constexpr size_t KDIM = 3;

    static constexpr std::array<double, NDIM> DOMAIN_MIN{
        -5,
        -5,
        -2,
        0,
        0
    };
    static constexpr std::array<double, NDIM> DOMAIN_MAX{
        5,
        5,
        2,
        2 * M_PI,
        2 * M_PI
    };
    static constexpr std::array<size_t, NDIM> DOMAIN_DIV{
        30,
        30,
        30,
        30,
        30,
    };

    static constexpr std::array<double, NDIM> FIRST_POINT{
        /*-3,   //,*/ -2.88283384,
        /*0,    //,*/ 0.08969913,
        /*0,    //,*/ -0.05940595,
        /*M_PI, //,*/ 3.11048777,
        /*M_PI  //,*/ 3.19927956
    };

    static void func(const double *x, double *f) {
        constexpr double a = 3;
        const double cos_half_x3 = std::cos(x[3] / 2);
        const double sin_x4 = std::sin(x[4]);
        const double sin_half_x3 = std::sin(x[3] / 2);
        const double sin_2_x4 = std::sin(2 * x[4]);
        const double aux = a + cos_half_x3 * sin_x4 - sin_half_x3 * sin_2_x4;
        f[0] = x[0] - aux * std::cos(x[3]);
        f[1] = x[1] - aux * std::sin(x[3]);
        f[2] = x[2] - sin_half_x3 * sin_x4 + cos_half_x3 * sin_2_x4;
    }
};
