///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#pragma once

#include "test_case_template.h"

template <>
struct TestCaseImpl<TestCase::S1_K_TIMES> {
    static constexpr size_t KDIM = 4;
    static constexpr size_t NDIM = KDIM * 2;

    static constexpr std::array<double, NDIM> DOMAIN_MIN = []() {
        std::array<double, NDIM> domain_min{};
        domain_min.fill(-1.1);
        return domain_min;
    }.operator()();

    static constexpr std::array<double, NDIM> DOMAIN_MAX = []() {
        std::array<double, NDIM> domain_max{};
        domain_max.fill(1.1);
        return domain_max;
    }.operator()();

    static constexpr std::array<size_t, NDIM> DOMAIN_DIV = []() {
        std::array<size_t, NDIM> domain_div{};
        domain_div.fill(5);
        return domain_div;
    }.operator()();

    static constexpr std::array<double, NDIM> FIRST_POINT = []() {
        std::array<double, NDIM> first_point{};
        first_point.fill(0);
        constexpr double STEP = 1. / static_cast<double>(KDIM + 1);
        for (size_t i = 0; i < KDIM; ++i) {
            const double x = STEP * static_cast<double>(i + 1);
            // doing binary search to find y = sqrt(y) = sqrt(1-x*x).
            double y = 1. - (x * x);
            double u = 1, d = 0, s = 0.5;
            bool finished = false;
            const double thresh = 1e-11;
            while (!finished) {
                s = (u + d) / 2.0;
                if ((s * s) > y) {
                    if (((s * s) - y) < thresh) {
                        finished = true;
                    } else {
                        u = s;
                    }
                } else {
                    if ((y - (s * s)) < thresh) {
                        finished = true;
                    } else {
                        d = s;
                    }
                }
            }
            y = s;
            first_point[2 * i] = x;
            first_point[2 * i + 1] = y;
        }
        return first_point;
    }.operator()();

    static void func(const double *x, double *f) {
        for (size_t i = 0; i < KDIM; ++i) {
            f[i] = x[2 * i] * x[2 * i] + x[2 * i + 1] * x[2 * i + 1] - 1.;
        }
    }
};
