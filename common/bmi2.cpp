///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#include "bmi2.h"

#if defined(__BMI2__)

uint32_t pmask(const uint32_t mask) {
    return mask;
}

uint64_t pdep(const uint64_t src, const uint32_t mask) {
    return _pdep_u64(src, mask);
}

uint64_t pext(const uint64_t src, const uint32_t mask) {
    return _pext_u64(src, mask);
}

#else

#include "zp7.h"

zp7_masks_64_t pmask(const uint32_t mask) {
    return zp7_ppp_64(mask);
}

uint64_t pdep(const uint64_t src, const zp7_masks_64_t &mask) {
    return zp7_pdep_pre_64(src, &mask);
}

uint64_t pext(const uint64_t src, const zp7_masks_64_t &mask) {
    return zp7_pext_pre_64(src, &mask);
}

#endif