///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#pragma once

#include "test_case_template.h"

template <>
struct TestCaseImpl<TestCase::HYPERSHPERE_N_TO_K> {
    static constexpr size_t NDIM = 7;
    static constexpr size_t KDIM = 2;

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
        domain_div.fill(10);
        return domain_div;
    }.operator()();

    static constexpr std::array<double, NDIM> ZERO = []() {
        // calculating the reference values for the hyperplane cuts for the coordinates that do not form the sphere.
        std::array<double, NDIM> zero{};
        zero.fill(0);
        for (size_t i = (NDIM - KDIM + 1); i < NDIM; ++i) {
            zero[i] = DOMAIN_MIN[i] + ((DOMAIN_MAX[i] - DOMAIN_MIN[i]) / DOMAIN_DIV[i]) * ((DOMAIN_DIV[i] >> 1) + 0.1 + 0.8 * (i - (NDIM - KDIM + 1)) / std::max(KDIM - 2.0, 1.0));
        }
        return zero;
    }.operator()();

    static constexpr std::array<double, NDIM> FIRST_POINT = []() {
        std::array<double, NDIM> first_point{};
        // filling the first point of the hypersphere.
        // first_point[1] -- first_point[NDIM - KDIM] are assigned directly.
        // first_point[0] is calculated to complete the radius of the sphere.
        first_point.fill(0);
        first_point[0] = 1.;
        for (size_t i = 1; i < (NDIM - KDIM + 1); ++i) {
            first_point[i] = 0.01 * static_cast<double>(i);  // because of the implementation, i will always be smaller than 100.
            first_point[0] -= (first_point[i] * first_point[i]);
        }
        // doing binary search to find the sqrt(first_point[0]).
        double u = 1, d = 0, s = 0.5;
        bool finished = false;
        const double thresh = 1e-11;
        while (!finished) {
            s = (u + d) / 2.0;
            if ((s * s) > first_point[0]) {
                if (((s * s) - first_point[0]) < thresh) {
                    finished = true;
                } else {
                    u = s;
                }
            } else {
                if ((first_point[0] - (s * s)) < thresh) {
                    finished = true;
                } else {
                    d = s;
                }
            }
        }
        first_point[0] = s;
        // filling the other coordinates.
        for (size_t i = (NDIM - KDIM + 1); i < NDIM; ++i) {
            first_point[i] = ZERO[i];
        }
        return first_point;
    }.operator()();

    static void func(const double *x, double *f) {
        f[0] = 1.;
        for (size_t i = 0; i < (NDIM - KDIM + 1); ++i) {
            f[0] -= x[i] * x[i];
        }
        for (size_t i = 1; i < KDIM; ++i) {
            f[i] = x[NDIM - KDIM + i] - ZERO[NDIM - KDIM + i];
        }
    }
};
