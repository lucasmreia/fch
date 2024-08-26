///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#include "utility.h"

#include <Eigen/Dense>
#include <cstring>

// Functions to convert bitset to/from coord.
void coordToBitset(const std::array<UINT_COORD, NDIM> &coord, BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH> &b) {
    b.reset();
    b |= coord[0];
    for (size_t i = 1; i < NDIM; ++i) {
        b <<= DOMAIN_DIV_BIT_WIDTH[i];
        b |= coord[i];
    }
}

void bitsetToCoord(const BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH> &b, std::array<UINT_COORD, NDIM> &coord) {
    BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH> aux1 = b, aux2;
    for (size_t i = 0; i < NDIM; ++i) {
        aux2 = aux1;
        aux2 &= ((1ULL << DOMAIN_DIV_BIT_WIDTH[NDIM - 1 - i]) - 1ULL);
        coord[NDIM - 1 - i] = aux2.to_ullong();
        aux1 >>= DOMAIN_DIV_BIT_WIDTH[NDIM - 1 - i];
    }
}

#ifdef PRECALCULATED_PERTURBATION
// Precalculated perturbation.
std::vector<double> createPerturbation() {
    std::vector<double> pert((1ULL << NDIM) * NDIM);
    std::seed_seq seed({ 0 });
    std::default_random_engine eng{ seed };
    std::uniform_real_distribution<double> urd(-1, 1);
    for (size_t i = 0; i < (1ULL << NDIM); ++i) {
        for (size_t j = 0; j < NDIM; ++j) {
            pert[i * NDIM + j] = urd(eng) * PERTURBATION_AMPLITUDE;
        }
    }
    return pert;
}
const std::vector<double> pert = createPerturbation();

// Function to add random perturbation to a vertex of the grid.
void addPerturbationVertex(const std::array<UINT_COORD, NDIM> &grid_coord, const size_t local_idx,
                           std::array<double, NDIM> &vertex) {
    // Getting the position of the vertex in the grid (integer index for each dimension).
    size_t idx = local_idx;
    for (size_t j = 0; j < NDIM; ++j) {
        // mod 2 is used to get the mirrored coords based on the first hypercube.
        idx ^= ((grid_coord[j] % 2) << j);
    }
    // Using the position in the grid as seed to generate the random perturbation.
    for (size_t j = 0; j < NDIM; ++j) {
        vertex[j] += pert[idx * NDIM + j];
    }
}
#else
// Function to add random perturbation to a vertex of the grid.
void addPerturbationVertex(const std::array<UINT_COORD, NDIM> &grid_coord, const size_t local_idx,
                           std::array<double, NDIM> &vertex) {
    // Getting the position of the vertex in the grid (integer index for each dimension).
    std::array<UINT_COORD, NDIM> seed_values = grid_coord;
    for (size_t j = 0; j < NDIM; ++j) {
        if (local_idx & (1ULL << j)) {
            ++seed_values[j];
        }
        seed_values[j] %= 2;  // mod 2 is used to get the mirrored coords based on the first hypercube.
    }
    // Using the position in the grid as seed to generate the random perturbation.
    std::seed_seq seed(seed_values.begin(), seed_values.end());
    std::default_random_engine eng{ seed };
    std::uniform_real_distribution<double> urd(-1, 1);
    for (size_t j = 0; j < NDIM; ++j) {
        const double pert = urd(eng) * PERTURBATION_AMPLITUDE;
        vertex[j] += pert;
    }
}
#endif

// Oracle to check if k-simples crosses the manifold, and compute the manifold intersecrtion vertex.
bool solve(const std::array<std::array<double, NDIM>, KDIM + 1> &vert_simplex_k,
           std::array<double, NDIM> &mf_vertex) {
    // assigning the positions of the vertices and evalVertices.
    // creating the matrices of the system A*lambda=b, V*lambda=mfVert.
    // eigen stores as column major.
    Eigen::Matrix<double, KDIM, KDIM> A{};
    Eigen::Vector<double, KDIM> b{};
    Eigen::Matrix<double, NDIM, KDIM + 1> V{};
    func(vert_simplex_k[0].data(), b.data());
    b *= -1;
    std::memcpy(&V(0, 0), vert_simplex_k[0].data(), NDIM * sizeof(double));
    for (int64_t j = 1; j < (KDIM + 1); ++j) {
        std::memcpy(&V(0, j), vert_simplex_k[j].data(), NDIM * sizeof(double));
        func(vert_simplex_k[j].data(), &A(0, j - 1));  // eigen has to be column major for this to work. check it.
        A.col(j - 1) += b;
    }

    // calculating the intersection of the manifold and the simplex.
    const Eigen::PartialPivLU<Eigen::Matrix<double, KDIM, KDIM>> lu_decomp = A.partialPivLu();
    const Eigen::Matrix<double, KDIM, KDIM> upper = lu_decomp.matrixLU().triangularView<Eigen::UpLoType::Upper>();

    for (int64_t i = 0; i < KDIM; ++i) {
        if (std::abs(upper(i, i)) < DIAGONAL_ZERO_THRESHOLD) {
            // not linearly independent.
            return false;
        }
    }

    // only works if A is linearly independent.
    Eigen::Vector<double, KDIM + 1> lambda = Eigen::Vector<double, KDIM + 1>::Zero();
    lambda.tail(KDIM) = lu_decomp.solve(b);
    lambda(0) = 1. - lambda.sum();
    // barycentric coordinates criterion.
    if (lambda.minCoeff() >= LAMBDA_ZERO_THRESHOLD) {
        // calculating the approximated position of the point where manifold intersects the simplex.
        Eigen::Map<Eigen::Vector<double, NDIM>> center{ mf_vertex.data() };
        center = V * lambda;
        return true;
    }

    return false;
}
