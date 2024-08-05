///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#pragma once

#include <iostream>

#include "../common/utility.h"

// saves in canonical notation, doesn't need the last vector ((KDIM+1)-1).
template <size_t DIM>
class BitsetPermutahedron : public BitLabel<BIT_WIDTH_GRID1_ENUM + (DIM * NDIM)> {
};

namespace std {
template <size_t DIM>
struct hash<BitsetPermutahedron<DIM>> {
    size_t operator()(const BitsetPermutahedron<DIM> &x) const {
        return std::hash<BitLabel<BIT_WIDTH_GRID1_ENUM + (DIM * NDIM)>>()(x);
    }
};
}  // namespace std

// saves in canonical notation, doesn't need the last vector ((KDIM+1)-1).
typedef BitsetPermutahedron<KDIM> BitsetPermutahedronSFace;

// saves in canonical notation, doesn't need the last vector ((KDIM+2)-1).
typedef BitsetPermutahedron<KDIM + 1> BitsetPermutahedronSCoface;

template <size_t DIM>
struct PermutahedronNotation {
    std::array<size_t, NDIM> grid_coord{};  // position of the reference vertex in the grid.
    std::array<size_t, DIM + 1> e{};        // k simplex: k+1 vectors. this is binary notation.
    // unlike hypercube vectors, here an element can be the sum of
    // more than one vector.

    virtual ~PermutahedronNotation() = default;

    virtual void print() const {
        std::cout << "(";
        for (size_t i = 0; i < NDIM; ++i) {
            std::cout << grid_coord[i] << " ";
        }
        std::cout << ") ";
        for (size_t i = 0; i < (DIM + 1); ++i) {
            std::cout << e[i] << " ";
        }
        std::cout << std::endl;
    }

    void fromBitset(const BitsetPermutahedron<DIM> &b) {
        BitsetPermutahedron<DIM> aux = b, temp;
        grid_coord.fill(0);
        e.fill(0);
        size_t last_e = (1ULL << (NDIM + 1)) - 1;
        // vectors are LSB and grid coords are MSB.
        // reads from last vector (LSB) to first vector (MSB).
        for (size_t i = 0; i < DIM; ++i) {
            temp = aux;
            temp &= ((1ULL << NDIM) - 1);
            e[DIM - 1 - i] = temp.to_ullong();
            aux >>= NDIM;
            last_e -= e[DIM - 1 - i];
        }
        e[DIM] = last_e;
        // grid coords are left.
        const size_t grid_enm = aux.to_ullong();
        enumToCoord(grid_enm, DOMAIN_GRID_BASIS1, grid_coord);
    }

    void toBitset(BitsetPermutahedron<DIM> &b) const {
        b.reset();
        // writes from grid to vectors, pushing the grid to the left (MSB).
        // in the end, vectors will be LSB and grid will be MSB.
        const size_t grid_coord_enm = coordToEnum(grid_coord, DOMAIN_GRID_BASIS1);
        b |= grid_coord_enm;
        // writes from first vector to last vector, pushing the first vectors to the left (MSB).
        // in the end, last vectors will be LSB and first vectors will be MSB.
        for (size_t i = 0; i < DIM; ++i) {
            b <<= NDIM;
            b |= e[i];
        }
        // in the end, vectors will be LSB and grid coords will be MSB.
    }

    void toCanonical(PermutahedronNotation<DIM> &canonical) const {
        size_t last_idx = -1, idx;
        for (size_t i = 0; i < (DIM + 1); ++i) {
            if (e[i] & (1ULL << NDIM)) {
                last_idx = i;
                break;
            }
        }
        if (last_idx == static_cast<size_t>(-1)) {
            std::cout << "error" << std::endl;
            print();
            throw std::runtime_error("face_to_canonical: last_idx cant be -1");
        }
        if (last_idx == DIM) {
            canonical = *this;
        } else {
            canonical.grid_coord = grid_coord;
            size_t acc_e = 0;  // accumulate the vectors that will be rotated.
            for (size_t i = 0; i < (DIM + 1); ++i) {
                idx = (last_idx + 1 + i) % (DIM + 1);
                if (idx > last_idx) {
                    acc_e |= e[idx];  // accumulating the e to move the grid.
                }
                canonical.e[i] = e[idx];
            }
            //-- doing: face_canonical.grid -= face.e[idx];
            for (size_t j = 0; j < NDIM; ++j) {
                if (acc_e & (1ULL << j)) {
                    --canonical.grid_coord[j];
                }
                if (acc_e & (1ULL << NDIM)) {
                    ++canonical.grid_coord[j];
                }
            }
            //--
        }
    }

    // simplex must be in the canonical notation for this to work properly.
    void toVertices(std::array<std::array<double, NDIM>, DIM + 1> &vertices) const {
        std::array<size_t, NDIM> aux_grid_coord = grid_coord;
        size_t curr_vidx = 0;
        for (size_t i = 0; i < (DIM + 1); ++i) {
            for (size_t j = 0; j < NDIM; ++j) {
                vertices[i][j] = DOMAIN_MIN[j] + static_cast<double>(aux_grid_coord[j]) * DOMAIN_RANGE[j];
                if (e[i] & (1ULL << j)) {
                    ++aux_grid_coord[j];
                }
            }
            // add perturbation to the vertices.
            addPerturbationVertex(grid_coord, curr_vidx, vertices[i]);
            curr_vidx |= e[i];
        }
    }
};

struct SCoface : public PermutahedronNotation<KDIM + 1> {
    void print() const override {
        std::cout << "face: ";
        PermutahedronNotation<KDIM + 1>::print();
    }
};

struct SFace : public PermutahedronNotation<KDIM> {
    void print() const override {
        std::cout << "coface: ";
        PermutahedronNotation<KDIM>::print();
    }
};

bool isSFaceInsideGrid(const SFace &face);

bool isSCofaceInsideGrid(const SCoface &coface);
