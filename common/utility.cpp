///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#include "utility.h"

#include <Eigen/Dense>
#include <cstring>

// Function to add random perturbation to a vertex of the grid.
void addPerturbationVertex(const std::array<size_t, NDIM> &grid_coord, const size_t local_idx,
                           std::array<double, NDIM> &vertex) {
    // Getting the position of the vertex in the grid (integer index for each dimension).
    std::array<size_t, NDIM> seed_values = grid_coord;
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

// Oracle to check if k-simples crosses the manifold, and compute the manifold intersecrtion vertex.
bool solve(const std::array<std::array<double, NDIM>, KDIM + 1> &vert_simplex_k,
           std::array<double, NDIM> &mf_vertex) {
    // calculating the positions of the vertices and evalVertices.
    // creating the matrices of the system A*lambda=b, V*lambda=mfVert.
    // eigen stores as column major.
    Eigen::Matrix<double, KDIM + 1, KDIM + 1> A{};
    Eigen::Matrix<double, NDIM, KDIM + 1> V{};
    for (int64_t j = 0; j < (KDIM + 1); ++j) {
        std::memcpy(&V(0, j), vert_simplex_k[j].data(), NDIM * sizeof(double));
        A(0, j) = 1;
        func(vert_simplex_k[j].data(), &A(1, j));  // eigen has to be column major for this to work. check it.
    }
    Eigen::Vector<double, KDIM + 1> b = Eigen::Vector<double, KDIM + 1>::Zero();
    b(0) = 1;

    // calculating the intersection of the manifold and the simplex.
    const Eigen::PartialPivLU<Eigen::Matrix<double, KDIM + 1, KDIM + 1>> lu_decomp = A.partialPivLu();
    const Eigen::Matrix<double, KDIM + 1, KDIM + 1> upper = lu_decomp.matrixLU().triangularView<Eigen::UpLoType::Upper>();

    for (int64_t i = 0; i < KDIM + 1; ++i) {
        if (std::abs(upper(i, i)) < DIAGONAL_ZERO_THRESHOLD) {
            // not linearly independent.
            return false;
        }
    }

    // only works if A is linearly independent.
    const Eigen::Vector<double, KDIM + 1> lambda = lu_decomp.solve(b);
    // barycentric coordinates criterion.
    if (lambda.minCoeff() >= 0) {
        // calculating the approximated position of the point where manifold intersects the simplex.
        Eigen::Map<Eigen::Vector<double, NDIM>> center{ mf_vertex.data() };
        center = V * lambda;
        return true;
    }

    return false;
}
