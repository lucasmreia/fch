///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#pragma once

#include "../common/bmi2.h"
#include "../common/utility.h"
#include "../gcch/gch.h"

// I have to perform this check because I do pdep with the label, and I also use to_ullong().
static_assert((2 * NDIM) < (sizeof(size_t) * 8));

struct HFace {
    std::array<size_t, NDIM> grid_coord{};
    Label2N label{};

    void toCanonical(HFace &canonical) const {
        canonical.grid_coord = grid_coord;
        canonical.label = label;
        for (size_t i = 0; i < NDIM; ++i) {
            if (canonical.label[i + NDIM] && (canonical.grid_coord[i] < (DOMAIN_DIV[i] - 1))) {
                canonical.label.reset(i + NDIM);
                canonical.label.set(i);
                canonical.grid_coord[i] += 1;
            }
        }
    }
};

constexpr bool checkInductiveBasis() {
    // if size_t can't store the grid size, it throws an exception.
    size_t aux = std::numeric_limits<size_t>::max();
    aux /= NDIM;
    for (size_t i = 1; i < NDIM; ++i) {
        if (aux < NDIM) {
            return false;
        }
        aux /= NDIM;
    }
    return true;
}

static_assert(checkInductiveBasis(), "inductive basis cant be stored");

constexpr std::array<size_t, NDIM + 1> createInductiveBasis() {
    std::array<size_t, NDIM + 1> basis{};
    basis[0] = 1;
    for (size_t i = 1; i <= NDIM; ++i) {
        basis[i] = basis[i - 1] * NDIM;
    }
    return basis;
}

constexpr std::array<size_t, NDIM + 1> INDUCTIVE_BASIS = createInductiveBasis();

template <size_t DIM>
class BitsetSimplex : public BitLabel<BIT_WIDTH_GRID_ENUM + (2 * NDIM) + std::bit_width(INDUCTIVE_BASIS[DIM])> {
};

namespace std {
template <size_t DIM>
struct hash<BitsetSimplex<DIM>> {
    size_t operator()(const BitsetSimplex<DIM> &x) const {
        return std::hash<BitLabel<BIT_WIDTH_GRID_ENUM + (2 * NDIM) + std::bit_width(INDUCTIVE_BASIS[DIM])>>()(x);
    }
};
}  // namespace std

template <size_t DIM>
struct Simplex {
    HFace hface{};
    std::array<size_t, DIM> e_idx{};  // vector indices, not the binary representation.

    void toCanonical(Simplex &canonical) const {
        canonical.e_idx = e_idx;
        hface.toCanonical(canonical.hface);
    }

    void toBitset(BitsetSimplex<DIM> &b) const {
        b.reset();
        b |= coordToEnum(hface.grid_coord, DOMAIN_GRID_BASIS);
        b <<= 2 * NDIM;
        b |= hface.label.to_ullong();
        b <<= std::bit_width(INDUCTIVE_BASIS[DIM]);
        b |= coordToEnum<DIM - 1>(e_idx.data(), INDUCTIVE_BASIS.data());
    }

    void fromBitset(const BitsetSimplex<DIM> &b) {
        BitsetSimplex<DIM> aux;
        aux = b;
        aux &= (1ULL << std::bit_width(INDUCTIVE_BASIS[DIM])) - 1;
        enumToCoord<DIM - 1>(aux.to_ullong(), INDUCTIVE_BASIS.data(), e_idx.data());
        aux = b;
        aux >>= std::bit_width(INDUCTIVE_BASIS[DIM]);
        aux &= (1ULL << (2 * NDIM)) - 1;
        hface.label.reset();
        hface.label |= aux.to_ullong();
        aux = b;
        aux >>= std::bit_width(INDUCTIVE_BASIS[DIM]) + (2 * NDIM);
        enumToCoord(aux.to_ullong(), DOMAIN_GRID_BASIS, hface.grid_coord);
        size_t fixed_coord, free_coord;
        hface.label.getFixedAndFreeCoord(fixed_coord, free_coord);
        for (size_t i = 0; i < (DIM - 1); ++i) {
            free_coord ^= (1ULL << e_idx[i]);
        }
        e_idx[DIM - 1] = std::countr_zero(free_coord);
    }

    void toVertices(std::array<std::array<double, NDIM>, DIM + 1> &vertices) const {
        const size_t first_vertex = (hface.label >> NDIM).to_ullong();
        std::array<size_t, NDIM> grid_coord_aux = hface.grid_coord;
        for (size_t i = 0; i < NDIM; ++i) {
            if (first_vertex & (1ULL << i)) {
                grid_coord_aux[i] += 1;
            }
            vertices[0][i] = DOMAIN_MIN[i] + DOMAIN_RANGE[i] * static_cast<double>(grid_coord_aux[i]);
        }
        size_t local_vertex = first_vertex;
        addPerturbationVertex(hface.grid_coord, local_vertex, vertices[0]);
        for (size_t i = 1; i <= DIM; ++i) {
            grid_coord_aux[e_idx[i - 1]] += 1;
            local_vertex |= (1ULL << e_idx[i - 1]);
            for (size_t j = 0; j < NDIM; ++j) {
                vertices[i][j] = DOMAIN_MIN[j] + DOMAIN_RANGE[j] * static_cast<double>(grid_coord_aux[j]);
            }
            addPerturbationVertex(hface.grid_coord, local_vertex, vertices[i]);
        }
    }
};

// I assume that bitIdx is a fixed coordinated.
// "Mirroring" actually means taking representations of the same simplex in different hypercubes,
// taking a neighboring hypercube in the direction of one of the fixed coordinates.
// since the direction of a fixed coordinate was chosen, the neighboring hypercube also shares this simplex.
// returns true if the representation exists (hypercube is inside the grid).
// bit_idx = direction of mirroring.
bool mirrorHFace(const std::array<size_t, NDIM> &in_grid, const Label2N &in_label,
                 size_t bit_idx,
                 std::array<size_t, NDIM> &out_grid, Label2N &out_label);

// vertices, Vert_simplex_k_1 and extra_vertex_idx will be modified as you pivot.
template <size_t DIM>
void traverseHCoface(const std::array<size_t, NDIM> &grid_coord,
                     std::array<std::array<double, NDIM>, KDIM + 2> &vertices,
                     std::array<size_t, DIM + 2> &vert_labels,
                     size_t &extra_vertex_idx,
                     std::array<double, NDIM> &mf_vertex) {
    // progressively solve the system and pivot until you find the output simplex.
    // create the system.
    Eigen::Matrix<double, DIM + 1, DIM + 1> A{};
    Eigen::Matrix<double, NDIM, DIM + 1> V{};
    std::array<size_t, DIM + 1> columns{};  // label of the vertex assigned to each column.

    // calculating the positions of the vertices and evalVertices.
    // creating the matrices of the A*x=b system.
    for (int64_t j = 0, v_idx = 0; j < (DIM + 1); ++j, ++v_idx) {
        if (v_idx == extra_vertex_idx) {
            ++v_idx;
        }
        columns[j] = v_idx;
        std::memcpy(&V(0, j), vertices[v_idx].data(), NDIM * sizeof(double));
        A(0, j) = 1;
        func(vertices[v_idx].data(), &A(1, j));
    }
    Eigen::Vector<double, DIM + 1> b_lambda = Eigen::Vector<double, DIM + 1>::Zero();
    b_lambda(0) = 1;

    // assigning b_mu.
    Eigen::Vector<double, DIM + 1> b_mu{};
    b_mu(0) = 1;
    func(vertices[extra_vertex_idx].data(), &b_mu(1));

    // variables for solving the system.
    Eigen::Vector<double, DIM + 1> lambda{};
    Eigen::PartialPivLU<Eigen::Matrix<double, DIM + 1, DIM + 1>> lu_decomp;

    std::array<size_t, NDIM> grid_coord_aux{};

    while (true) {
        // solve the system to calculate lambda.
        lu_decomp = A.partialPivLu();
        lambda = lu_decomp.solve(b_lambda);

        const int64_t pivot_column_idx = getPivotColumnIdx<DIM>(lu_decomp, lambda, b_mu);
        if ((columns[pivot_column_idx] == 0) || (columns[pivot_column_idx] == (DIM + 1))) {
            // end of pivoting.
            // create the system with the extra vertex in place of the one to be replaced.
            std::swap(columns[pivot_column_idx], extra_vertex_idx);
            std::memcpy(&V(0, pivot_column_idx), vertices[columns[pivot_column_idx]].data(),
                        NDIM * sizeof(double));
            func(vertices[columns[pivot_column_idx]].data(), &A(1, pivot_column_idx));
            // solve the system to calculate lambda.
            lu_decomp = A.partialPivLu();
            lambda = lu_decomp.solve(b_lambda);

            // calculating the approximated position where the manifold intersects the simplex and return.
            Eigen::Map<Eigen::Vector<double, NDIM>> center{ mf_vertex.data() };
            center = V * lambda;
            return;
        }

        const size_t previous_v_idx = columns[pivot_column_idx] - 1;
        const size_t next_v_idx = columns[pivot_column_idx] + 1;
        const size_t dif_e = vert_labels[next_v_idx] - vert_labels[previous_v_idx];
        const size_t pivoted_v_label = vert_labels[columns[pivot_column_idx]] ^ dif_e;

        // assign a new pivoted vertex.
        vert_labels[columns[pivot_column_idx]] = pivoted_v_label;

        // calculate the position of the new pivoted vertex.
        for (size_t i = 0; i < NDIM; ++i) {
            grid_coord_aux[i] = grid_coord[i];
            if (pivoted_v_label & (1ULL << i)) {
                grid_coord_aux[i] += 1;
            }
            vertices[columns[pivot_column_idx]][i] =
                DOMAIN_MIN[i] + static_cast<double>(grid_coord_aux[i]) * DOMAIN_RANGE[i];
        }

        // add perturbation.
        addPerturbationVertex(grid_coord, pivoted_v_label, vertices[columns[pivot_column_idx]]);

        // replace the pivot_column_idx column with the extra_vertex_idx.
        std::swap(columns[pivot_column_idx], extra_vertex_idx);

        // updating b_mu with the new pivoted vertex.
        func(vertices[extra_vertex_idx].data(), &b_mu(1));

        // updating the column with the replaced vertex.
        std::memcpy(&V(0, pivot_column_idx), vertices[columns[pivot_column_idx]].data(),
                    NDIM * sizeof(double));
        func(vertices[columns[pivot_column_idx]].data(), &A(1, pivot_column_idx));
    }
}

template <size_t DIM>
void computeSimplexLabelAndInductive(const std::array<size_t, DIM + 1> &vert_labels,
                                     Label2N &hface_label,
                                     std::array<size_t, DIM> &e_idx) {
    // binary representation of the face's simplex vectors.
    std::array<size_t, DIM> e_hface{};

    // assigning part of the face of the k-simplex.
    hface_label.reset();
    for (size_t j = 0; j < NDIM; ++j) {
        if (vert_labels[0] & (1ULL << j)) {
            hface_label.set(j + NDIM);
        }
    }

    // assigning the vectors of the k-simplex.
    size_t fixed_coord = (1ULL << NDIM) - 1;
    for (size_t j = 0; j < DIM; ++j) {
        e_hface[j] = vert_labels[j + 1] - vert_labels[j];
        fixed_coord -= e_hface[j];
    }

    // assigning the vectors of the face.
    for (size_t j = 0; j < DIM; ++j) {
        if (std::popcount(e_hface[j]) != 1) {
            throw std::runtime_error("wrong number of vectors");
        }
        e_idx[j] = std::countr_zero(e_hface[j]);
    }

    // assigning the rest of the label of the face.
    for (size_t j = 0; j < NDIM; ++j) {
        if ((fixed_coord & (1ULL << j)) && !hface_label[j + NDIM]) {
            hface_label.set(j);
        }
    }
}

void writeOutputCell(size_t hcoface_g,
                     const Label2N &hcoface_label,
                     size_t hface_g_in,
                     const std::array<size_t, KDIM + 1> &simplex_vert_label_in,
                     const std::array<double, NDIM> &mf_vert_in,
                     size_t hface_g_out,
                     const std::array<size_t, KDIM + 1> &simplex_vert_label_out,
                     const std::array<double, NDIM> &mf_vert_out,
                     FILE *fout);

void readOutputCell(FILE *fin,
                    Label2N &hcoface_label,
                    size_t &hface_g_in,
                    std::array<size_t, KDIM + 1> &simplex_vert_label_in,
                    std::array<double, NDIM> &mf_vert_in,
                    size_t &hface_g_out,
                    std::array<size_t, KDIM + 1> &simplex_vert_label_out,
                    std::array<double, NDIM> &mf_vert_out);
