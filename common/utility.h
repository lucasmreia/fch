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

#include "definitions.h"

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

constexpr bool checkDomainDiv() {
    // if size_t can't store the grid size, it throws an exception.
    size_t aux = std::numeric_limits<size_t>::max();
    aux /= DOMAIN_DIV[0];
    for (size_t i = 1; i < NDIM; ++i) {
        if (aux < DOMAIN_DIV[i]) {
            return false;
        }
        aux /= DOMAIN_DIV[i];
    }
    return true;
}

static_assert(checkDomainDiv(), "domain has too many divisions");

// Dimension of the manifold.
constexpr size_t MDIM = NDIM - KDIM;

// Cell length at a given dimension (step).
constexpr std::array<double, NDIM> createDomainRange() {
    std::array<double, NDIM> range{};
    for (size_t i = 0; i < NDIM; ++i) {
        range[i] = (DOMAIN_MAX[i] - DOMAIN_MIN[i]) / static_cast<double>(DOMAIN_DIV[i]);
    }
    return range;
}

constexpr std::array<double, NDIM> DOMAIN_RANGE = createDomainRange();

// considering only the inner part of the grid.
constexpr std::array<size_t, NDIM + 1> createDomainGridBasis() {
    std::array<size_t, NDIM + 1> basis{};
    basis[0] = 1;
    for (size_t i = 1; i <= NDIM; ++i) {
        basis[i] = basis[i - 1] * DOMAIN_DIV[i - 1];
    }
    return basis;
}

constexpr std::array<size_t, NDIM + 1> DOMAIN_GRID_BASIS = createDomainGridBasis();
constexpr size_t BIT_WIDTH_GRID_ENUM = std::bit_width(DOMAIN_GRID_BASIS[NDIM]);

// considering that the grid has 1 extra point at the end of each dimension (outer boundary).
constexpr std::array<size_t, NDIM + 1> createDomainGridBasis1() {
    std::array<size_t, NDIM + 1> basis{};
    basis[0] = 1;
    for (size_t i = 1; i <= NDIM; ++i) {
        basis[i] = basis[i - 1] * (DOMAIN_DIV[i - 1] + 1);
    }
    return basis;
}

constexpr std::array<size_t, NDIM + 1> DOMAIN_GRID_BASIS1 = createDomainGridBasis1();
constexpr size_t BIT_WIDTH_GRID1_ENUM = std::bit_width(DOMAIN_GRID_BASIS1[NDIM]);

template <size_t SIZE>
size_t coordToEnum(const size_t *coord, const size_t *basis) {
    size_t enm = 0;
    for (size_t i = 0; i < SIZE; ++i) {
        enm += coord[i] * basis[i];
    }
    return enm;
}

template <size_t SIZE>
size_t coordToEnum(const std::array<size_t, SIZE> &coord, const std::array<size_t, SIZE + 1> &basis) {
    return coordToEnum<SIZE>(coord.data(), basis.data());
}

template <size_t SIZE>
void enumToCoord(const size_t enm, const size_t *basis, size_t *coord) {
    size_t aux = enm;
    for (size_t i = 0; i < SIZE; ++i) {
        coord[SIZE - 1 - i] = aux / basis[SIZE - 1 - i];
        aux %= basis[SIZE - 1 - i];
    }
}

template <size_t SIZE>
void enumToCoord(const size_t enm, const std::array<size_t, SIZE + 1> &basis, std::array<size_t, SIZE> &coord) {
    enumToCoord<SIZE>(enm, basis.data(), coord.data());
}

// calculate which cube the first point belongs to.
constexpr std::array<size_t, NDIM> getFirstPointHypercubeCoord() {
    std::array<size_t, NDIM> grid_coord{};
    for (size_t i = 0; i < NDIM; ++i) {
        grid_coord[i] = static_cast<size_t>(std::max(0., FIRST_POINT[i] - DOMAIN_MIN[i]) / DOMAIN_RANGE[i]);
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

// Function to add random perturbation to a vertex of the grid.
void addPerturbationVertex(const std::array<size_t, NDIM> &grid_coord, size_t local_idx,
                           std::array<double, NDIM> &vertex);

// Oracle to check if k-simples crosses the manifold, and compute the manifold intersecrtion vertex.
bool solve(const std::array<std::array<double, NDIM>, KDIM + 1> &vert_simplex_k,
           std::array<double, NDIM> &mf_vertex);
