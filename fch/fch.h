///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#pragma once

#include <stdexcept>

#include "../common/utility.h"

// Uses only the NDIM lsb because the msb are always zero.
class CanonicalFaceLabel : public BitLabel<NDIM> {
};

// Bitmask to mark the cofaces yet to be traversed.
class CofacesBitmask : public BitLabel<2 * NDIM> {
};

template <size_t DIM>
class BitsetSimplex : public BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH + INDUCTIVE_INDEX_BIT_WIDTH * DIM> {
};

namespace std {
template <size_t DIM>
struct hash<BitsetSimplex<DIM>> {
    size_t operator()(const BitsetSimplex<DIM> &x) const {
        return std::hash<BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH + INDUCTIVE_INDEX_BIT_WIDTH * DIM>>()(x);
    }
};
}  // namespace std

// Equivalent to the canonical permutahedron representation, but it is valid only for
// DIM-simplices inside hypercube DIM-faces.
template <size_t DIM>
struct CanonicalSimplex {
    std::array<UINT_COORD, NDIM> grid_coord{};
    std::array<UINT_EIDX, DIM> e_idx{};  // vector indices, not the binary representation.

    void toBitset(BitsetSimplex<DIM> &b) const {
        b.reset();
        b |= grid_coord[0];
        for (size_t i = 1; i < NDIM; ++i) {
            b <<= DOMAIN_DIV_BIT_WIDTH[i];
            b |= grid_coord[i];
        }
        b <<= INDUCTIVE_INDEX_BIT_WIDTH;
        b |= e_idx[0];
        for (size_t i = 1; i < DIM; ++i) {
            b <<= INDUCTIVE_INDEX_BIT_WIDTH;
            b |= e_idx[i];
        }
    }

    static std::array<BitsetSimplex<DIM>, NDIM> getMaskCoord() {
        std::array<BitsetSimplex<DIM>, NDIM> mask{};
        for (size_t i = 0; i < NDIM; ++i) {
            mask[i].reset();
            mask[i] |= ((1ULL << DOMAIN_DIV_BIT_WIDTH[i]) - 1ULL);
        }
        return mask;
    }

    void fromBitset(const BitsetSimplex<DIM> &b) {
        BitsetSimplex<DIM> aux1 = b, aux2;
        for (size_t i = 0; i < DIM; ++i) {
            aux2 = aux1;
            aux2 &= ((1ULL << INDUCTIVE_INDEX_BIT_WIDTH) - 1ULL);
            e_idx[DIM - 1 - i] = aux2.to_ullong();
            aux1 >>= INDUCTIVE_INDEX_BIT_WIDTH;
        }
        for (size_t i = 0; i < NDIM; ++i) {
            aux2 = aux1;
            aux2 &= ((1ULL << DOMAIN_DIV_BIT_WIDTH[NDIM - 1 - i]) - 1ULL);
            grid_coord[NDIM - 1 - i] = aux2.to_ullong();
            aux1 >>= DOMAIN_DIV_BIT_WIDTH[NDIM - 1 - i];
        }
    }

    void toCanonicalFaceLabel(CanonicalFaceLabel &label) const {
        label.set();
        for (size_t i = 0; i < DIM; ++i) {
            label.reset(e_idx[i]);
        }
    }

    void toVertices(std::array<std::array<double, NDIM>, DIM + 1> &vertices) const {
        // in the canonical notation, the coord is the minimal point of the simplex.
        std::array<UINT_COORD, NDIM> grid_coord_aux = grid_coord;
        for (size_t i = 0; i < NDIM; ++i) {
            vertices[0][i] = DOMAIN_MIN[i] + DOMAIN_STEP[i] * static_cast<double>(grid_coord_aux[i]);
        }
        size_t local_vertex = 0;
        addPerturbationVertex(grid_coord, local_vertex, vertices[0]);
        for (size_t i = 1; i <= DIM; ++i) {
            grid_coord_aux[e_idx[i - 1]] += 1;
            local_vertex |= (1ULL << e_idx[i - 1]);
            for (size_t j = 0; j < NDIM; ++j) {
                vertices[i][j] = DOMAIN_MIN[j] + DOMAIN_STEP[j] * static_cast<double>(grid_coord_aux[j]);
            }
            addPerturbationVertex(grid_coord, local_vertex, vertices[i]);
        }
    }

    void toVertLabelsAndVertices(std::array<size_t, DIM + 1> &vert_labels,
                                 std::array<std::array<double, NDIM>, DIM + 1> &vertices) const {
        // in the canonical notation, the coord is the minimal point of the simplex.
        std::array<UINT_COORD, NDIM> grid_coord_aux = grid_coord;
        for (size_t i = 0; i < NDIM; ++i) {
            vertices[0][i] = DOMAIN_MIN[i] + DOMAIN_STEP[i] * static_cast<double>(grid_coord_aux[i]);
        }
        size_t local_vertex = 0;
        vert_labels[0] = local_vertex;
        addPerturbationVertex(grid_coord, local_vertex, vertices[0]);
        for (size_t i = 1; i <= DIM; ++i) {
            grid_coord_aux[e_idx[i - 1]] += 1;
            local_vertex |= (1ULL << e_idx[i - 1]);
            vert_labels[i] = local_vertex;
            for (size_t j = 0; j < NDIM; ++j) {
                vertices[i][j] = DOMAIN_MIN[j] + DOMAIN_STEP[j] * static_cast<double>(grid_coord_aux[j]);
            }
            addPerturbationVertex(grid_coord, local_vertex, vertices[i]);
        }
    }
};

template <size_t DIM>
void simplexFromVertLabels(const std::array<UINT_COORD, NDIM> &grid_coord,
                           const std::array<size_t, DIM + 1> &vert_labels,
                           CanonicalSimplex<DIM> &simplex) {
    // assigning the reference vertex.
    simplex.grid_coord = grid_coord;
    for (size_t j = 0; j < NDIM; ++j) {
        if (vert_labels[0] & (1ULL << j)) {
            simplex.grid_coord[j] += 1;
        }
    }

    // assigning the vectors.
    size_t e_bin;
    for (size_t j = 0; j < DIM; ++j) {
        e_bin = vert_labels[j + 1] - vert_labels[j];
        if (std::popcount(e_bin) != 1) {
            throw std::runtime_error("std::popcount(e_bin) must be 1");
        }
        simplex.e_idx[j] = std::countr_zero(e_bin);
    }
}

void writeOutputCell(const size_t g,
                     const BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH> &bitset_coord,
                     const CanonicalFaceLabel &hcoface_label,
                     const BitsetSimplex<KDIM> &bitset_simplex_in,
                     const std::array<double, NDIM> &mf_vert_in,
                     const BitsetSimplex<KDIM> &bitset_simplex_out,
                     const std::array<double, NDIM> &mf_vert_out,
                     FILE *fout);

void readOutputCell(FILE *fin,
                    BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH> &bitset_coord,
                    CanonicalFaceLabel &hcoface_label,
                    BitsetSimplex<KDIM> &bitset_simplex_in,
                    std::array<double, NDIM> &mf_vert_in,
                    BitsetSimplex<KDIM> &bitset_simplex_out,
                    std::array<double, NDIM> &mf_vert_out);
