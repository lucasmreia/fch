///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#pragma once

// Includes.
#include <array>
#include <cmath>
#include <cstddef>

#include "test_cases/all_cases.h"

// Definitions.
#ifndef SELECTED_CASE
#define SELECTED_CASE KLEIN_5TO3
#endif

// #define PRINT_DEBUG
#define PRECALCULATED_PERTURBATION
#define PTA_DOOR_IN_DOOR_OUT

// Error thresholds.
constexpr double PERTURBATION_AMPLITUDE = 1e-8;    // amplitude of the perturbation applied to the vertices.
constexpr double DIAGONAL_ZERO_THRESHOLD = 1e-11;  // threshold used to check if all diagonal elements of a matrix are zero.
constexpr double LAMBDA_ZERO_THRESHOLD = 0;        // threshold used to check if all lambda are greater than zero.

// Parameters.
constexpr size_t NDIM = TestCaseImpl<TestCase::SELECTED_CASE>::NDIM;
constexpr size_t KDIM = TestCaseImpl<TestCase::SELECTED_CASE>::KDIM;

// Range (domain).
constexpr std::array<double, NDIM> DOMAIN_MIN = TestCaseImpl<TestCase::SELECTED_CASE>::DOMAIN_MIN;
constexpr std::array<double, NDIM> DOMAIN_MAX = TestCaseImpl<TestCase::SELECTED_CASE>::DOMAIN_MAX;
constexpr std::array<size_t, NDIM> DOMAIN_DIV = TestCaseImpl<TestCase::SELECTED_CASE>::DOMAIN_DIV;

// First point for the continuation methods.
constexpr std::array<double, NDIM> FIRST_POINT = TestCaseImpl<TestCase::SELECTED_CASE>::FIRST_POINT;

// Implicit function evaluation.
void func(const double *x, double *f);
