///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#pragma once

#if defined(__BMI2__)

#include <immintrin.h>

#include <cinttypes>

typedef uint32_t P_MASK;

uint32_t pmask(uint32_t mask);

uint64_t pdep(uint64_t src, uint32_t mask);

uint64_t pext(uint64_t src, uint32_t mask);

#else

#include <cinttypes>

#define N_BITS (6)

typedef struct {
    uint64_t mask;
    uint64_t ppp_bit[N_BITS];
} zp7_masks_64_t;

typedef zp7_masks_64_t P_MASK;

zp7_masks_64_t pmask(uint32_t mask);

uint64_t pdep(uint64_t src, const zp7_masks_64_t &mask);

uint64_t pext(uint64_t src, const zp7_masks_64_t &mask);

#endif