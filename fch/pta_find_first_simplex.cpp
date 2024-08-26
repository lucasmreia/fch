///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#include "pta_find_first_simplex.h"

#include <Eigen/Dense>
#include <algorithm>
#include <argparse/argparse.hpp>
#include <array>
#include <bit>
#include <chrono>
#include <iostream>
#include <unordered_set>

#include "../common/bmi2.h"
#include "../pta/permutahedron.h"

#ifndef PTA_DOOR_IN_DOOR_OUT
bool createFacesFromCoface(std::unordered_set<BitsetPermutahedronSFace> &faces_to_be_processed,
                           std::unordered_set<BitsetPermutahedronSFace> &faces_already_processed,
                           const SCoface &coface,
                           CanonicalSimplex<KDIM> &found_simplex) {
    // print_coface(coface);
    SFace face;
    BitsetPermutahedronSFace bitset_face;
    std::array<std::array<double, NDIM>, KDIM + 1> vertices{};
    std::array<double, NDIM> mf_vertex{};
    size_t idx;
    // combining the vectors two by two, always taking two in sequence.
    for (size_t i = 0; i < (KDIM + 1); ++i) {  // choosing which coface vector to combine with the next.
        face.grid_coord = coface.grid_coord;
        idx = 0;
        for (size_t j = 0; j < (KDIM + 1); ++j) {  // filling the face.
            if (j != i) {
                face.e[j] = coface.e[idx];
                ++idx;
            } else {
                face.e[j] = coface.e[idx] + coface.e[idx + 1];
                idx += 2;
            }
        }
        //-- error check.
        if (!(face.e[KDIM] & (1ULL << NDIM))) {
            std::cout << "error" << std::endl;
            face.print();
            throw std::runtime_error("last face vector must contain e_{NDIM}");
        }
        //--
        if (isSFaceInsideGrid(face)) {
            face.toBitset(bitset_face);
            if (!faces_already_processed.contains(bitset_face)) {
                face.toVertices(vertices);
                bool contains_manifold = solve(vertices, mf_vertex);
                if (contains_manifold) {
                    faces_to_be_processed.insert(bitset_face);
                    faces_already_processed.insert(bitset_face);
                    // returning the simplex found if it is aligned with the grid.
                    if (std::popcount(face.e[KDIM]) == (NDIM + 1 - KDIM)) {
                        found_simplex.grid_coord = face.grid_coord;
                        for (size_t ii = 0; ii < KDIM; ++ii) {
                            if (std::popcount(face.e[ii]) != 1) {
                                throw std::runtime_error("std::popcount(face.e[ii]) must be 1");
                            }
                            found_simplex.e_idx[ii] = std::countr_zero(face.e[ii]);
                        }
                        return true;
                    }
                }
            }
        }
    }
    // combining the last vector with the first.
    // I have to change the reference vertex: advance 1 vector.
    face.grid_coord = coface.grid_coord;
    for (size_t j = 0; j < NDIM; ++j) {   // coface.grid += face.e[0], only considers canonical notation.
        if (coface.e[0] & (1ULL << j)) {  // e[0] will never have the negative vector (e_n).
            ++face.grid_coord[j];
        }
    }
    for (size_t i = 0; i < (KDIM + 1); ++i) {
        face.e[i] = coface.e[i + 1];
    }
    face.e[KDIM] += coface.e[0];

    //-- error check.
    if (!(face.e[KDIM] & (1ULL << NDIM))) {
        std::cout << "error" << std::endl;
        face.print();
        throw std::runtime_error("last face vector must contain e_{NDIM}");
    }
    //--
    if (isSFaceInsideGrid(face)) {
        face.toBitset(bitset_face);
        if (!faces_already_processed.contains(bitset_face)) {
            face.toVertices(vertices);
            bool contains_manifold = solve(vertices, mf_vertex);
            if (contains_manifold) {
                faces_to_be_processed.insert(bitset_face);
                faces_already_processed.insert(bitset_face);
                // returning the simplex found if it is aligned with the grid.
                if (std::popcount(face.e[KDIM]) == (NDIM + 1 - KDIM)) {
                    found_simplex.grid_coord = face.grid_coord;
                    for (size_t ii = 0; ii < KDIM; ++ii) {
                        if (std::popcount(face.e[ii]) != 1) {
                            throw std::runtime_error("std::popcount(face.e[ii]) must be 1");
                        }
                        found_simplex.e_idx[ii] = std::countr_zero(face.e[ii]);
                    }
                    return true;
                }
            }
        }
    }
    return false;
}

bool createCofacesFromFace(std::unordered_set<BitsetPermutahedronSFace> &faces_to_be_processed,
                           std::unordered_set<BitsetPermutahedronSFace> &faces_already_processed,
                           const SFace &face,
                           CanonicalSimplex<KDIM> &found_simplex) {
    // print_face(face);
    SCoface coface, coface_canonical;
    coface.grid_coord = face.grid_coord;
    size_t pcount, idx = 0, part_a, part_b;
    P_MASK pdep_mask;
    for (size_t i = 0; i < (KDIM + 1); ++i) {  // looping through the face vectors to see which one can be split.
        pcount = std::popcount(face.e[i]);     // popcount gives the number of vectors added together.
        if (pcount <= 1) {                     // if it can't be divided, just assign it.
            coface.e[idx] = face.e[i];
            ++idx;
        } else {
            pdep_mask = pmask(face.e[i]);  // the mask for pdep is the vector itself.
            // looping from 000...0001 to 111...1110 to split the groups of vectors.
            for (size_t j = 1; j <= ((1ULL << pcount) - 2); ++j) {
                part_a = j;  // separating in two groups, group a and group b.
                part_b = part_a ^ ((1ULL << pcount) - 1);
                coface.e[idx] = pdep(part_a, pdep_mask);         // new vertex added to create the (k+1)-simplex.
                coface.e[idx + 1] = pdep(part_b, pdep_mask);     // vertex that already existed in the k-simplexo.
                for (size_t m = 1; m < ((KDIM + 1) - i); ++m) {  // filling in the following elements, only assigns.
                    coface.e[(idx + 1) + m] = face.e[i + m];
                }
                coface.toCanonical(coface_canonical);
                //-- error check.
                if (!(coface_canonical.e[KDIM + 1] & (1ULL << NDIM))) {
                    std::cout << "error" << std::endl;
                    coface.print();
                    throw std::runtime_error("last coface vector must contain e_{NDIM}");
                }
                //--
                if (isSCofaceInsideGrid(coface_canonical)) {
                    if (createFacesFromCoface(faces_to_be_processed,
                                              faces_already_processed,
                                              coface_canonical,
                                              found_simplex)) {
                        return true;
                    }
                }
            }
            coface.e[idx] = face.e[i];  // now just assign this one to split the next ones.
            ++idx;
        }
    }
    return false;
}
#else
bool traverseCofacesFromFace(std::unordered_set<BitsetPermutahedronSFace> &faces_to_be_processed,
                             std::unordered_set<BitsetPermutahedronSFace> &faces_already_processed,
                             const SFace &in_face,
                             CanonicalSimplex<KDIM> &found_simplex) {
    // print_face(face);
    SCoface coface, coface_canonical;
    coface.grid_coord = in_face.grid_coord;
    size_t pcount, idx = 0, part_a, part_b;

    Eigen::Matrix<double, KDIM, KDIM> A{};
    Eigen::Matrix<double, NDIM, KDIM + 1> V{};
    std::array<std::array<double, NDIM>, KDIM + 2> vertices_simplex_k_1{};
    size_t coface_vert_idx;

    Eigen::Vector<double, KDIM + 1> lambda{};
    Eigen::Vector<double, KDIM + 1> mu{};
    Eigen::Vector<double, KDIM> b_lambda{};
    Eigen::Vector<double, KDIM> b_mu{};
    Eigen::PartialPivLU<Eigen::Matrix<double, KDIM, KDIM> > lu_decomp;

    SFace out_face;
    BitsetPermutahedronSFace out_bitset_face;

    std::array<size_t, KDIM + 1> vlabel_simplex_k{};
    std::array<size_t, KDIM + 2> vlabel_simplex_k_1{};
    std::array<size_t, KDIM + 1> vlabel_simplex_k_aux{};

    std::array<double, NDIM> mf_vertex{};
    Eigen::Map<Eigen::Vector<double, NDIM> > center{ mf_vertex.data() };

    P_MASK pdep_mask;

    for (size_t i = 0; i < (KDIM + 1); ++i) {  // looping through the face vectors to see which one can be split.
        pcount = std::popcount(in_face.e[i]);  // popcount gives the number of vectors added together.
        if (pcount <= 1) {                     // if it can't be divided, just assign it.
            coface.e[idx] = in_face.e[i];
            ++idx;
        } else {
            pdep_mask = pmask(in_face.e[i]);  // the mask for pdep is the vector itself.
            // looping from 000...0001 to 111...1110 to split the groups of vectors.
            for (size_t j = 1; j <= ((1ULL << pcount) - 2); ++j) {
                part_a = j;  // separating in two groups, group a and group b.
                part_b = part_a ^ ((1ULL << pcount) - 1);
                coface.e[idx] = pdep(part_a, pdep_mask);         // new vertex added to create the (k+1)-simplex.
                coface.e[idx + 1] = pdep(part_b, pdep_mask);     // vertex that already existed in the k-simplexo.
                for (size_t m = 1; m < ((KDIM + 1) - i); ++m) {  // filling in the following elements, only assigns.
                    coface.e[(idx + 1) + m] = in_face.e[i + m];
                }
                coface.toCanonical(coface_canonical);
                //-- error check.
                if (!(coface_canonical.e[KDIM + 1] & (1ULL << NDIM))) {
                    std::cout << "error" << std::endl;
                    coface.print();
                    throw std::runtime_error("last coface vector must contain e_{NDIM}");
                }
                //--
                if (isSCofaceInsideGrid(coface_canonical)) {
                    // create the vlabels of both the face and the coface, and the vertices of the coface.
                    in_face.toVertLabels(vlabel_simplex_k);
                    coface_canonical.toVertLabelsAndVertices(vlabel_simplex_k_1, vertices_simplex_k_1);

                    // the reference grid_coord of the vlabels of both the face and the coface must be the same.
                    // adjusting the vlabels of the face if needed.
                    if (coface_canonical.grid_coord != in_face.grid_coord) {
                        size_t off = 0;
                        for (size_t ii = 0; ii < NDIM; ++ii) {
                            // face grid_coord must always be greater or equal to coface grid_coord.
                            if (in_face.grid_coord[ii] > coface_canonical.grid_coord[ii]) {
                                off += 1ULL << ii;
                            } else if (in_face.grid_coord[ii] < coface_canonical.grid_coord[ii]) {
                                throw std::runtime_error("coface_canonical.grid must be less or equal face.grid");
                            }
                        }
                        for (size_t ii = 0; ii < (KDIM + 1); ++ii) {
                            vlabel_simplex_k[ii] += off;
                        }
                    }

                    // identify which coface vertex is not in the face.
                    coface_vert_idx = 0;
                    for (size_t face_vert_idx = 0; face_vert_idx < (KDIM + 1); ++face_vert_idx) {
                        if (vlabel_simplex_k_1[coface_vert_idx] == vlabel_simplex_k[face_vert_idx]) {
                            ++coface_vert_idx;
                        } else {
                            break;
                        }
                    }

                    // creating the systems of equations used by the door-in/door-out principle.
                    // assigning the positions of the vertices and evalVertices.
                    for (int64_t jj = 0, v_idx = 0; jj < (KDIM + 1); ++jj, ++v_idx) {
                        if (v_idx == coface_vert_idx) {
                            ++v_idx;
                        }
                        std::memcpy(&V(0, jj), vertices_simplex_k_1[v_idx].data(), NDIM * sizeof(double));
                        if (jj == 0) {
                            func(vertices_simplex_k_1[v_idx].data(), b_lambda.data());
                            b_lambda *= -1;
                        } else {
                            func(vertices_simplex_k_1[v_idx].data(), &A(0, jj - 1));
                            A.col(jj - 1) += b_lambda;
                        }
                    }

                    // assigning b_mu.
                    func(vertices_simplex_k_1[coface_vert_idx].data(), b_mu.data());
                    b_mu += b_lambda;

                    // solve system to calculate both lambda e mu.
                    lu_decomp = A.partialPivLu();

                    lambda(0) = 0;
                    lambda.tail(KDIM) = lu_decomp.solve(b_lambda);
                    lambda(0) = 1. - lambda.sum();

                    mu(0) = 0;
                    mu.tail(KDIM) = lu_decomp.solve(b_mu);
                    mu(0) = 1. - mu.sum();

                    // calcualte min{lambda/mu} and the index of the vertex that must be replaced.
                    double min_lambda_mu = 0;
                    int64_t replaced_column_idx = -1;
                    for (int64_t ii = 0; ii < (KDIM + 1); ++ii) {
                        if (mu(ii) > 0 && ((replaced_column_idx < 0) || (min_lambda_mu > lambda(ii) / mu(ii)))) {
                            replaced_column_idx = ii;
                            min_lambda_mu = lambda(ii) / mu(ii);
                        }
                    }
                    //-- error check.
                    if (replaced_column_idx < 0) {
                        std::cout << "error replaced_column_idx: " << replaced_column_idx << std::endl;
                        throw std::runtime_error("replaced_column_idx must be >= 0");
                    }
                    //--

                    // determine the output face.
                    vlabel_simplex_k_aux = vlabel_simplex_k;
                    vlabel_simplex_k_aux[replaced_column_idx] = vlabel_simplex_k_1[coface_vert_idx];
                    std::sort(vlabel_simplex_k_aux.begin(), vlabel_simplex_k_aux.end());

                    out_face.grid_coord = coface_canonical.grid_coord;
                    // ajust the out face grid_coord to go back to canonical.
                    if (vlabel_simplex_k_aux[0] > 0) {
                        for (size_t ii = 0; ii < NDIM; ++ii) {
                            if (vlabel_simplex_k_aux[0] & (1ULL << ii)) {
                                out_face.grid_coord[ii] += 1;
                            }
                        }
                        size_t off = vlabel_simplex_k_aux[0];
                        for (size_t ii = 0; ii < (KDIM + 1); ++ii) {
                            vlabel_simplex_k_aux[ii] -= off;
                        }
                    }

                    out_face.e.fill(0);
                    out_face.e[KDIM] = (1ULL << (NDIM + 1)) - 1;
                    for (int64_t ii = 0; ii < KDIM; ++ii) {
                        out_face.e[ii] = vlabel_simplex_k_aux[ii + 1] - vlabel_simplex_k_aux[ii];
                        out_face.e[KDIM] -= out_face.e[ii];
                    }

                    //-- error check.
                    if (!(out_face.e[KDIM] & (1ULL << NDIM))) {
                        out_face.print();
                        throw std::runtime_error("last out_face vector must contain e_{NDIM}");
                    }
                    //--
                    if (isSFaceInsideGrid(out_face)) {
                        out_face.toBitset(out_bitset_face);
                        if (!faces_already_processed.contains(out_bitset_face)) {
                            // replace the vertex in the matrix.
                            // updating the column with the new vertex.
                            std::memcpy(&V(0, replaced_column_idx), vertices_simplex_k_1[coface_vert_idx].data(),
                                        NDIM * sizeof(double));

                            // updating the system of equations.
                            if (replaced_column_idx > 0) {
                                // if it is not the first vertex, just update the replaced column.
                                func(vertices_simplex_k_1[coface_vert_idx].data(), &A(0, replaced_column_idx - 1));
                                A.col(replaced_column_idx - 1) += b_lambda;
                            } else {
                                // if it is the first vertex, recreate all matrices.
                                func(&V(0, 0), b_lambda.data());
                                b_lambda *= -1;
                                for (int64_t jj = 1; jj < (KDIM + 1); ++jj) {
                                    func(&V(0, jj), &A(0, jj - 1));  // eigen has to be column major for this to work. check it.
                                    A.col(jj - 1) += b_lambda;
                                }
                            }

                            // get output lambda.
                            lu_decomp = A.partialPivLu();
                            lambda(0) = 0;
                            lambda.tail(KDIM) = lu_decomp.solve(b_lambda);
                            lambda(0) = 1. - lambda.sum();

                            // check if output vertex was found.
                            if (!(lambda.minCoeff() >= LAMBDA_ZERO_THRESHOLD)) {
                                throw std::runtime_error("door-in,door-out could not find output k-simplex");
                            }

                            // add simplex to already processed/to be processed.
                            faces_to_be_processed.insert(out_bitset_face);
                            faces_already_processed.insert(out_bitset_face);
                            // returning the simplex found if it is aligned with the grid.
                            if (std::popcount(out_face.e[KDIM]) == (NDIM + 1 - KDIM)) {
                                found_simplex.grid_coord = out_face.grid_coord;
                                for (size_t ii = 0; ii < KDIM; ++ii) {
                                    if (std::popcount(out_face.e[ii]) != 1) {
                                        throw std::runtime_error("std::popcount(face.e[ii]) must be 1");
                                    }
                                    found_simplex.e_idx[ii] = std::countr_zero(out_face.e[ii]);
                                }
                                return true;
                            }
                        }
                    }
                }
            }
            coface.e[idx] = in_face.e[i];  // now just assign this one to split the next ones.
            ++idx;
        }
    }
    return false;
}
#endif

void addFirstPoint(std::unordered_set<BitsetPermutahedronSFace> &faces_to_be_processed,
                   std::unordered_set<BitsetPermutahedronSFace> &faces_already_processed) {
    // calculate which n-hypercube the first point belongs to.
    constexpr std::array<UINT_COORD, NDIM> GRID = getFirstPointHypercubeCoord();

    // coordinates of the vertices of a simplex, to check for intersection with the manifold.
    std::array<std::array<double, NDIM>, KDIM + 1> vertices{};
    std::array<double, NDIM> mf_vertex{};

    // determine, in the permutahedral representation, which n-simplex the point belongs to.
    std::array<double, NDIM + 1> aux{};
    std::array<size_t, NDIM + 1> e_idx{};
    for (size_t i = 0; i < NDIM; ++i) {
        // taking the percentage of each coordinate of the point inside the unit hypercube.
        aux[i] = (FIRST_POINT[i] - (DOMAIN_MIN[i] + static_cast<double>(GRID[i]) * DOMAIN_STEP[i])) / DOMAIN_STEP[i];
        e_idx[i] = i;
    }
    aux[NDIM] = 0;
    e_idx[NDIM] = NDIM;
    // printing aux to help with debugging if no point is found.
    std::cout << "FIRST_POINT coords inside the unitary hypercube: (";
    for (size_t i = 0; i < NDIM; ++i) {
        std::cout << aux[i];
        if (i < (NDIM - 1)) {
            std::cout << ", ";
        }
    }
    std::cout << ")" << std::endl;
    std::cout << "if two coords are equal, might not find the manifold (may get the wrong n-simplex)" << std::endl;
    // argsort: rearranges e_idx by sorting aux in decreasing order.
    std::stable_sort(e_idx.begin(), e_idx.end() - 1,  // excludes the last one.
                     [&aux](int left, int right) -> bool {
                         // sort indices according to corresponding array element.
                         return aux[left] > aux[right];
                     });
    // now e_idx represents the indices of the ordered vectors of the n-simplex that contains the first point.
    // converting it to the binary representation.
    for (size_t i = 0; i < (NDIM + 1); ++i) {
        e_idx[i] = (1ULL << e_idx[i]);
    }

    // now decompose the n-simplex in k-simplices and mark the faces to be processed.
    // n-simplex has n+1 vectors. k-simplex has k+1 vectors. I need to combine (n+1)-(k+1)=n-k vectors to
    // reduce the dimension.
    // if I represent the vertices of the cube as indices, I can do nchoosek(n+1, k+1), and then go back to the
    // permutahedral representation.
    // calculating the vertices of the n-simplex.
    std::array<size_t, NDIM + 1> simplex_n_orig{};
    simplex_n_orig[0] = 0;
    for (size_t i = 1; i < (NDIM + 1); ++i) {
        simplex_n_orig[i] = simplex_n_orig[i - 1] | e_idx[i - 1];
    }

    std::array<size_t, NDIM + 1> simplex_n{};

    // for the case when the first simplex is not found, try to pivot to neighboring vertices.
    for (size_t piv_idx = 0; piv_idx < NDIM; ++piv_idx) {
        simplex_n = simplex_n_orig;
        // first time doesn't pivot. next times pivot on piv_idx.
        if (piv_idx > 0) {
            simplex_n[piv_idx] ^= (simplex_n[piv_idx + 1] - simplex_n[piv_idx - 1]);
        }

        // now I need to do nchoosek(n+1, k+1).
        // making n choose k (actually (n+1) choose (k+1)) and taking the face.
        SFace face, face_canonical;
        BitsetPermutahedronSFace bitset_face;
        std::array<size_t, KDIM + 1> idx{};
        for (size_t i = 0; i < (KDIM + 1); ++i) {
            idx[i] = i;
        }
        while (idx[0] <= ((NDIM + 1) - (KDIM + 1))) {
            // assigning the initial vertex of k-simplex.
            face.grid_coord = GRID;
            for (size_t j = 0; j < NDIM; ++j) {
                if (simplex_n[idx[0]] & (1ULL << j)) {
                    ++face.grid_coord[j];
                }
            }

            // assigning the k-simplex vectors.
            size_t last_e = (1ULL << (NDIM + 1)) - 1;
            for (size_t j = 0; j < KDIM; ++j) {
                face.e[j] = simplex_n[idx[j + 1]] - simplex_n[idx[j]];
                last_e -= face.e[j];
            }
            face.e[KDIM] = last_e;

            face.toCanonical(face_canonical);

            if (isSFaceInsideGrid(face_canonical)) {
                face_canonical.toVertices(vertices);
                bool contains_manifold = solve(vertices, mf_vertex);
                if (contains_manifold) {
                    face_canonical.toBitset(bitset_face);
                    faces_to_be_processed.insert(bitset_face);
                    return;
                }
            }

            int64_t to_be_increased = (KDIM + 1) - 1;
            bool reset = false;
            while ((to_be_increased >= 0) && (++idx[to_be_increased] > ((NDIM + 1) - (KDIM + 1) + to_be_increased))) {
                --to_be_increased;
                reset = true;
            }
            if (reset && (to_be_increased >= 0)) {
                for (size_t i = to_be_increased + 1; i < (KDIM + 1); ++i) {
                    idx[i] = idx[i - 1] + 1;
                }
            }
        }
    }
}

bool find_first_simplex(CanonicalSimplex<KDIM> &found_simplex) {
    std::unordered_set<BitsetPermutahedronSFace> faces_to_be_processed;
    std::unordered_set<BitsetPermutahedronSFace> faces_already_processed;

    addFirstPoint(faces_to_be_processed, faces_already_processed);

    BitsetPermutahedronSFace bitset_face;
    SFace face;

    while (!faces_to_be_processed.empty()) {
        bitset_face = *faces_to_be_processed.begin();
        faces_to_be_processed.erase(faces_to_be_processed.begin());

        // convert the label of the k-simplex to coordinates.
        face.fromBitset(bitset_face);

        // calculate the cofaces.
        // then, for each coface, calculate the faces.
#ifndef PTA_DOOR_IN_DOOR_OUT
        if (createCofacesFromFace(faces_to_be_processed,
                                  faces_already_processed,
                                  face,
                                  found_simplex)) {
            return true;
        }
#else
        if (traverseCofacesFromFace(faces_to_be_processed,
                                    faces_already_processed,
                                    face,
                                    found_simplex)) {
            return true;
        }
#endif
    }
    return false;
}