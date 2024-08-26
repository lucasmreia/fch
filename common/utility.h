///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#pragma once

#include <bit>
#include <bitset>
#include <cstdio>
#include <limits>
#include <random>
#include <type_traits>
#include <typeinfo>

#include "definitions.h"

// Limiting the dimension of the domain.
// In all methods, it must be able to address all 2^NDIM vertices of a hypercube.
// This cannot exceed the maximum size of size_t.
// The amount of RAM will run out well before that.
// In the PTA, it must be able to do ((1ULL << (NDIM + 1)) - 1) to create the last vector.
// I sum 2 because if I shift by (NDIM + 1) there has to be two extra bits at the end.
// So (NDIM + 2) must be within the amount of bits available.
static_assert((NDIM + 2) < (sizeof(size_t) * 8));

constexpr bool checkFirstPoint() {
    // if the initial cube is out of range, throw exception.
    for (size_t i = 0; i < NDIM; ++i) {
        if ((FIRST_POINT[i] < DOMAIN_MIN[i]) || (FIRST_POINT[i] > DOMAIN_MAX[i])) {
            return false;
        }
    }
    return true;
}
static_assert(checkFirstPoint(), "first point out of range");

// Dimension of the manifold.
constexpr size_t MDIM = NDIM - KDIM;

// Cell length at a given dimension (step).
constexpr std::array<double, NDIM> createDomainStep() {
    std::array<double, NDIM> step{};
    for (size_t i = 0; i < NDIM; ++i) {
        step[i] = (DOMAIN_MAX[i] - DOMAIN_MIN[i]) / static_cast<double>(DOMAIN_DIV[i]);
    }
    return step;
}
constexpr std::array<double, NDIM> DOMAIN_STEP = createDomainStep();

// Required bit width for the divisions of each coordinate of the domain.
constexpr std::array<size_t, NDIM> createDomainDivBitWidth() {
    std::array<size_t, NDIM> bit_width{};
    for (size_t i = 0; i < NDIM; ++i) {
        bit_width[i] = std::bit_width(DOMAIN_DIV[i]);
    }
    return bit_width;
}
constexpr std::array<size_t, NDIM> DOMAIN_DIV_BIT_WIDTH = createDomainDivBitWidth();

// Total required bit width for the divisions of the domain.
constexpr size_t createDomainDivTotalBitWidth() {
    size_t total_bit_width = 0;
    for (size_t i = 0; i < NDIM; ++i) {
        total_bit_width += DOMAIN_DIV_BIT_WIDTH[i];
    }
    return total_bit_width;
}
constexpr size_t DOMAIN_DIV_TOTAL_BIT_WIDTH = createDomainDivTotalBitWidth();

// Max required bit width for one of the divisions of the domain.
constexpr size_t createDomainDivMaxBitWidth() {
    size_t max_bit_width = 0;
    for (size_t i = 0; i < NDIM; ++i) {
        max_bit_width = std::max(max_bit_width, DOMAIN_DIV_BIT_WIDTH[i]);
    }
    return max_bit_width;
}
constexpr size_t DOMAIN_DIV_MAX_BIT_WIDTH = createDomainDivMaxBitWidth();

// Required bit width for the index of each vector of the inductive notation.
constexpr size_t INDUCTIVE_INDEX_BIT_WIDTH = std::bit_width(NDIM - 1);

// Integer type used by the coords.
using UINT_COORD = std::conditional<(DOMAIN_DIV_MAX_BIT_WIDTH > 32), uint64_t, std::conditional<(DOMAIN_DIV_MAX_BIT_WIDTH > 16), uint32_t, std::conditional<(DOMAIN_DIV_MAX_BIT_WIDTH > 8), uint16_t, uint8_t>::type>::type>::type;

// Integer type used by the vector indices.
using UINT_EIDX = std::conditional<(INDUCTIVE_INDEX_BIT_WIDTH > 32), uint64_t, std::conditional<(INDUCTIVE_INDEX_BIT_WIDTH > 16), uint32_t, std::conditional<(INDUCTIVE_INDEX_BIT_WIDTH > 8), uint16_t, uint8_t>::type>::type>::type;

// Integer type used by the vectors in binary format.
using UINT_EBIN = std::conditional<((NDIM + 1) > 32), uint64_t, std::conditional<((NDIM + 1) > 16), uint32_t, std::conditional<((NDIM + 1) > 8), uint16_t, uint8_t>::type>::type>::type;

// calculate which cube the first point belongs to.
constexpr std::array<UINT_COORD, NDIM> getFirstPointHypercubeCoord() {
    std::array<UINT_COORD, NDIM> grid_coord{};
    for (size_t i = 0; i < NDIM; ++i) {
        grid_coord[i] = static_cast<UINT_COORD>(std::max(0., FIRST_POINT[i] - DOMAIN_MIN[i]) / DOMAIN_STEP[i]);
    }
    return grid_coord;
}

// Utility class: bitset with functions for writing and reading to/from a binary file.
template <size_t BIT_SIZE>
class BitLabel : public std::bitset<BIT_SIZE> {
   public:
    static constexpr size_t bitSize() {
        return BIT_SIZE;
    }

    static constexpr size_t byteSize() {
        return (BIT_SIZE + 7) / 8;
    }

    void write(FILE *fout) const {
        constexpr size_t BYTE_SIZE = byteSize();
        std::bitset<BIT_SIZE> aux = (*this), temp;
        uint8_t byte;
        // writing byte by byte, starting from the least significant.
        for (size_t i = 0; i < BYTE_SIZE; ++i) {
            temp = aux;
            temp &= 0xFF;
            byte = temp.to_ulong();
            fwrite(&byte, 1, 1, fout);
            aux >>= 8;
        }
    }

    void read(FILE *fin) {
        constexpr size_t BYTE_SIZE = byteSize();
        std::bitset<BIT_SIZE> aux;
        uint8_t byte;
        this->reset();
        // reading byte by byte, starting from the least significant.
        for (size_t i = 0; i < BYTE_SIZE; ++i) {
            fread(&byte, 1, 1, fin);
            aux.reset();
            aux |= byte;
            aux <<= (8ULL * i);
            (*this) |= aux;
        }
    }

    bool operator<(const BitLabel<BIT_SIZE> &other) const {
        for (int64_t i = (BIT_SIZE - 1); i >= 0; --i) {
            if ((*this)[i] ^ other[i]) return other[i];
        }
        return false;
    }
};

namespace std {
template <size_t BIT_SIZE>
struct hash<BitLabel<BIT_SIZE>> {
    size_t operator()(const BitLabel<BIT_SIZE> &x) const {
        return std::hash<std::bitset<BIT_SIZE>>()(x);
    }
};
}  // namespace std

// Functions to convert bitset to/from coord.
void coordToBitset(const std::array<UINT_COORD, NDIM> &coord, BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH> &b);

void bitsetToCoord(const BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH> &b, std::array<UINT_COORD, NDIM> &coord);

// Function to add random perturbation to a vertex of the grid.
void addPerturbationVertex(const std::array<UINT_COORD, NDIM> &grid_coord, size_t local_idx,
                           std::array<double, NDIM> &vertex);

// Oracle to check if k-simples crosses the manifold, and compute the manifold intersecrtion vertex.
bool solve(const std::array<std::array<double, NDIM>, KDIM + 1> &vert_simplex_k,
           std::array<double, NDIM> &mf_vertex);
