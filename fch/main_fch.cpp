///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#include <Eigen/Dense>
#include <algorithm>
#include <argparse/argparse.hpp>
#include <chrono>
#include <deque>
#include <iostream>
#include <memory>
#include <unordered_set>

#include "fch.h"
#include "pta_find_first_simplex.h"

class DequeUniqueElements {
   private:
    std::deque<std::set<BitsetSimplex<KDIM>>::const_iterator> q;
    std::set<BitsetSimplex<KDIM>> s;

   public:
    void push_back(const BitsetSimplex<KDIM> &simplex) {
        const auto p = s.insert(simplex);  // std::pair<std::set<Foo>::iterator, bool>
        if (p.second) {
            q.push_back(p.first);
        }
    }

    BitsetSimplex<KDIM> pop_front() {
        const auto iter = q.front();
        const BitsetSimplex<KDIM> f = *iter;
        s.erase(iter);
        q.pop_front();
        return f;
    }

    size_t size() {
        return q.size();
    }

    bool empty() {
        return q.empty();
    }
};

class UnorderedDequeUniqueElements : public std::unordered_set<BitsetSimplex<KDIM>> {
   public:
    void push_back(const BitsetSimplex<KDIM> &simplex) {
        insert(simplex);
    }

    BitsetSimplex<KDIM> pop_front() {
        const auto iter = begin();
        const BitsetSimplex<KDIM> f = *iter;
        erase(iter);
        return f;
    }
};

std::unordered_map<BitsetSimplex<KDIM>, CofacesBitmask> cofaces_yet_to_be_processed;
DequeUniqueElements tau_to_be_processed;

size_t n_saved_cells = 0;

void computeOutputTauAndSigma(const std::array<UINT_COORD, NDIM> &grid_coord,
                              const std::array<size_t, KDIM + 2> &sigma_vert_labels,
                              const size_t extra_vertex_idx,
                              CanonicalSimplex<KDIM> &tau,
                              CanonicalSimplex<KDIM + 1> &sigma) {
    // --- error check.
    if (!std::is_sorted(sigma_vert_labels.begin(), sigma_vert_labels.end())) {
        std::cout << "sigma_vertices must be sorted" << std::endl;
        throw std::runtime_error("sigma_vertices not sorted");
    }
    // ---

    std::array<size_t, KDIM + 1> tau_vert_labels{};
    for (size_t i = 0, v_idx = 0; i < (KDIM + 1); ++i, ++v_idx) {
        if (v_idx == extra_vertex_idx) {
            ++v_idx;
        }
        tau_vert_labels[i] = sigma_vert_labels[v_idx];
    }

    simplexFromVertLabels(grid_coord,
                          tau_vert_labels,
                          tau);

    simplexFromVertLabels(grid_coord,
                          sigma_vert_labels,
                          sigma);
}

// Get the index of the column that has to be replaced.
// I already pass lu_decomp and lambda to reuse those values.
int64_t getReplacedColumnIdx(const Eigen::PartialPivLU<Eigen::Matrix<double, KDIM, KDIM>> &lu_decomp,
                             const Eigen::Vector<double, KDIM + 1> &lambda,
                             const Eigen::Vector<double, KDIM> &b_mu) {
    // solve the system to calculate mu.
    Eigen::Vector<double, KDIM + 1> mu = Eigen::Vector<double, KDIM + 1>::Zero();
    mu.tail(KDIM) = lu_decomp.solve(b_mu);
    mu(0) = 1. - mu.sum();

    // calcualte min{lambda/mu} and the index of the vertex that must be replaced.
    double min_lambda_mu = 0;
    int64_t idx_k = -1;
    for (int64_t ii = 0; ii < (KDIM + 1); ++ii) {
        if (mu(ii) > 0 && ((idx_k < 0) || (min_lambda_mu > lambda(ii) / mu(ii)))) {
            idx_k = ii;
            min_lambda_mu = lambda(ii) / mu(ii);
        }
    }
    //-- error check.
    if (idx_k < 0) {
        std::cout << "error idx_k: " << idx_k << std::endl;
        throw std::runtime_error("idx_k must be >= 0");
    }
    //--
    return idx_k;
}

// sigma_vert_labels, sigma_vertices and extra_vertex_idx will be modified as you pivot.
void traverseHCoface(const std::array<UINT_COORD, NDIM> &grid_coord,
                     std::array<size_t, KDIM + 2> &sigma_vert_labels,
                     std::array<std::array<double, NDIM>, KDIM + 2> &sigma_vertices,
                     size_t &extra_vertex_idx,
                     std::array<double, NDIM> &mf_vert_in,
                     std::array<double, NDIM> &mf_vertex_out) {
    // progressively solve the system and pivot until you find the output simplex.
    // create the system.
    Eigen::Matrix<double, KDIM, KDIM> A{};
    Eigen::Matrix<double, NDIM, KDIM + 1> V{};
    std::array<size_t, KDIM + 1> columns{};  // label of the vertex assigned to each column.

    Eigen::Vector<double, KDIM> b_lambda{};

    // calculating the positions of the vertices and evalVertices.
    // creating the matrices of the A*x=b system.
    for (int64_t j = 0, v_idx = 0; j < (KDIM + 1); ++j, ++v_idx) {
        if (v_idx == extra_vertex_idx) {
            ++v_idx;
        }
        columns[j] = v_idx;
        std::memcpy(&V(0, j), sigma_vertices[v_idx].data(), NDIM * sizeof(double));
        if (j == 0) {
            func(sigma_vertices[v_idx].data(), b_lambda.data());
            b_lambda *= -1;
        } else {
            func(sigma_vertices[v_idx].data(), &A(0, j - 1));
            A.col(j - 1) += b_lambda;
        }
    }

    // assigning b_mu.
    Eigen::Vector<double, KDIM> b_mu{};
    func(sigma_vertices[extra_vertex_idx].data(), b_mu.data());
    b_mu += b_lambda;

    // variables for solving the system.
    Eigen::Vector<double, KDIM + 1> lambda{};
    Eigen::PartialPivLU<Eigen::Matrix<double, KDIM, KDIM>> lu_decomp;

    // solve the system to calculate lambda for the input k-simplex.
    lu_decomp = A.partialPivLu();
    lambda(0) = 0;
    lambda.tail(KDIM) = lu_decomp.solve(b_lambda);
    lambda(0) = 1. - lambda.sum();

    if (!(lambda.minCoeff() >= LAMBDA_ZERO_THRESHOLD)) {
        throw std::runtime_error("could not find input k-simplex");
    }

    // calculating the approximated position where the manifold intersects the input simplex.
    {
        Eigen::Map<Eigen::Vector<double, NDIM>> center{ mf_vert_in.data() };
        center = V * lambda;
    }

    std::array<size_t, NDIM> grid_coord_aux{};

    while (true) {
        const int64_t replaced_column_idx = getReplacedColumnIdx(lu_decomp, lambda, b_mu);

        // replace the column with the extra_vertex_idx.
        std::swap(columns[replaced_column_idx], extra_vertex_idx);

        // updating the column with the new vertex.
        std::memcpy(&V(0, replaced_column_idx), sigma_vertices[columns[replaced_column_idx]].data(),
                    NDIM * sizeof(double));

        // updating the system of equations.
        if (replaced_column_idx > 0) {
            // if it is not the first vertex, just update the replaced column.
            func(sigma_vertices[columns[replaced_column_idx]].data(), &A(0, replaced_column_idx - 1));
            A.col(replaced_column_idx - 1) += b_lambda;
        } else {
            // if it is the first vertex, recreate all matrices.
            func(&V(0, 0), b_lambda.data());
            b_lambda *= -1;
            for (int64_t j = 1; j < (KDIM + 1); ++j) {
                func(&V(0, j), &A(0, j - 1));  // eigen has to be column major for this to work. check it.
                A.col(j - 1) += b_lambda;
            }
        }

        // solve the system to calculate lambda.
        lu_decomp = A.partialPivLu();
        lambda(0) = 0;
        lambda.tail(KDIM) = lu_decomp.solve(b_lambda);
        lambda(0) = 1. - lambda.sum();

        if (!(lambda.minCoeff() >= LAMBDA_ZERO_THRESHOLD)) {
            throw std::runtime_error("door-in,door-out could not find output k-simplex");
        }

        // pivoting the (k+1)-simplex on extra_vertex_idx to find the new extra vertex.
        // however, if it wants to pivot on the first or the last vertex of the (k+1)-simplex, if means the manifold
        // has left the (k+1)-hypercube, so the traversal ends.
        if ((extra_vertex_idx == 0) || (extra_vertex_idx == (KDIM + 1))) {
            // end of the traversal inside the coface.
            // calculating the approximated position where the manifold intersects the output k-simplex and return.
            Eigen::Map<Eigen::Vector<double, NDIM>> center{ mf_vertex_out.data() };
            center = V * lambda;
            return;
        }

        // pivoting the (k+1)-simplex on extra_vertex_idx to find the new extra vertex.
        const size_t previous_v_idx = extra_vertex_idx - 1;
        const size_t next_v_idx = extra_vertex_idx + 1;
        const size_t dif_e = sigma_vert_labels[next_v_idx] - sigma_vert_labels[previous_v_idx];
        const size_t pivoted_v_label = sigma_vert_labels[extra_vertex_idx] ^ dif_e;

        // assign a new pivoted vertex.
        sigma_vert_labels[extra_vertex_idx] = pivoted_v_label;

        // calculate the position of the new pivoted vertex.
        for (size_t i = 0; i < NDIM; ++i) {
            grid_coord_aux[i] = grid_coord[i];
            if (pivoted_v_label & (1ULL << i)) {
                grid_coord_aux[i] += 1;
            }
            sigma_vertices[extra_vertex_idx][i] =
                DOMAIN_MIN[i] + static_cast<double>(grid_coord_aux[i]) * DOMAIN_STEP[i];
        }

        // add perturbation.
        addPerturbationVertex(grid_coord, pivoted_v_label, sigma_vertices[extra_vertex_idx]);

        // updating b_mu with the new pivoted vertex.
        func(sigma_vertices[extra_vertex_idx].data(), b_mu.data());
        b_mu += b_lambda;
    }
}

void createCosimplicesFromSimplex(const BitsetSimplex<KDIM> &bitset_tau_in,
                                  const CanonicalSimplex<KDIM> &tau_in,
                                  FILE *fout) {
    std::array<double, NDIM> mf_vert_in{};

    CanonicalSimplex<KDIM + 1> sigma_in{}, sigma_out{};
    BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH> bitset_sigma_coord{};

    CanonicalSimplex<KDIM> tau_out{};
    std::array<double, NDIM> mf_vert_out{};
    BitsetSimplex<KDIM> bitset_tau_out;

    std::array<size_t, KDIM + 2> sigma_vert_labels{};
    std::array<std::array<double, NDIM>, KDIM + 2> sigma_vertices{};
    CanonicalFaceLabel label_tau_out, label_sigma_out;
    size_t bitmask_traversed_out;
    size_t bit_idx_out;

    CofacesBitmask yet_to_be_processed;
    {
        auto iter_in = cofaces_yet_to_be_processed.find(bitset_tau_in);
        yet_to_be_processed = iter_in->second;
        cofaces_yet_to_be_processed.erase(iter_in);
    }

    std::array<UINT_EIDX, NDIM - KDIM> tau_fixed_coord{};
    {
        CanonicalFaceLabel label;
        tau_in.toCanonicalFaceLabel(label);
        size_t counter = 0;
        for (size_t i = 0; i < NDIM; ++i) {
            if (label[i]) {
                tau_fixed_coord[counter++] = i;
            }
        }
    }

    for (size_t i = 0; i < (NDIM - KDIM); ++i) {
        const size_t unfix_coord_idx = tau_fixed_coord[i];
        if (yet_to_be_processed[unfix_coord_idx]) {
            yet_to_be_processed.reset(unfix_coord_idx);
            // checking if the coface is inside the domain.
            if (tau_in.grid_coord[unfix_coord_idx] < DOMAIN_DIV[unfix_coord_idx]) {
                // computing coface_in.
                sigma_in.grid_coord = tau_in.grid_coord;
                std::memcpy(sigma_in.e_idx.data(), tau_in.e_idx.data(), KDIM * sizeof(UINT_EIDX));
                sigma_in.e_idx[KDIM] = unfix_coord_idx;
                // call TraverseCoface and assign mfVert_out.
                size_t extra_vertex_idx = KDIM + 1;
                sigma_out = sigma_in;
                sigma_out.toVertLabelsAndVertices(sigma_vert_labels, sigma_vertices);
                traverseHCoface(sigma_out.grid_coord,
                                sigma_vert_labels,
                                sigma_vertices,
                                extra_vertex_idx,
                                mf_vert_in,
                                mf_vert_out);
                // compute tau out and sigma out in the canonical notation.
                computeOutputTauAndSigma(sigma_in.grid_coord,
                                         sigma_vert_labels,
                                         extra_vertex_idx,
                                         tau_out,
                                         sigma_out);
                // compute their labels.
                tau_out.toCanonicalFaceLabel(label_tau_out);
                sigma_out.toCanonicalFaceLabel(label_sigma_out);
                // get the bitmask indicating the coface traversed from tau out.
                bitmask_traversed_out = label_tau_out.to_ullong();
                bitmask_traversed_out ^= label_sigma_out.to_ullong();
                bit_idx_out = std::countr_zero(bitmask_traversed_out);
                if (tau_out.grid_coord != sigma_out.grid_coord) {
                    bit_idx_out += NDIM;
                }
                // add the output tau to faces_to_be_processed.
                tau_out.toBitset(bitset_tau_out);
                tau_to_be_processed.push_back(bitset_tau_out);
                // add the output sigma to cofaces_yet_to_be_processed.
                {
                    auto iter_out = cofaces_yet_to_be_processed.find(bitset_tau_out);
                    if (iter_out == cofaces_yet_to_be_processed.end()) {
                        const size_t fixed_coord = label_tau_out.to_ullong();
                        CofacesBitmask all_neighbors;
                        all_neighbors |= fixed_coord;
                        all_neighbors <<= NDIM;
                        all_neighbors |= fixed_coord;
                        iter_out = cofaces_yet_to_be_processed.emplace(bitset_tau_out, all_neighbors).first;
                    } else if (!iter_out->second[bit_idx_out]) {
                        throw std::runtime_error("output simplex cannot be already processed");
                    }
                    iter_out->second.reset(bit_idx_out);
                }
                // save edge: mfVert_in -> mfVert_out, coface (grid, label, vertices, the same way marching hypercubes).
                coordToBitset(sigma_out.grid_coord, bitset_sigma_coord);
                ++n_saved_cells;
                writeOutputCell(n_saved_cells,
                                bitset_sigma_coord,
                                label_sigma_out,
                                bitset_tau_in,
                                mf_vert_in,
                                bitset_tau_out,
                                mf_vert_out,
                                fout);
            }
        }

        if (yet_to_be_processed[unfix_coord_idx + NDIM]) {
            yet_to_be_processed.reset(unfix_coord_idx + NDIM);
            // checking if the coface is inside the domain.
            if (tau_in.grid_coord[unfix_coord_idx] > 0) {
                // computing coface_in.
                sigma_in.grid_coord = tau_in.grid_coord;
                sigma_in.grid_coord[unfix_coord_idx] -= 1;
                std::memcpy((sigma_in.e_idx.data() + 1), tau_in.e_idx.data(), KDIM * sizeof(UINT_EIDX));
                sigma_in.e_idx[0] = unfix_coord_idx;
                // call TraverseCoface and assign mfVert_out.
                size_t extra_vertex_idx = 0;
                sigma_out = sigma_in;
                sigma_out.toVertLabelsAndVertices(sigma_vert_labels, sigma_vertices);
                traverseHCoface(sigma_out.grid_coord,
                                sigma_vert_labels,
                                sigma_vertices,
                                extra_vertex_idx,
                                mf_vert_in,
                                mf_vert_out);
                // compute tau out and sigma out in the canonical notation.
                computeOutputTauAndSigma(sigma_in.grid_coord,
                                         sigma_vert_labels,
                                         extra_vertex_idx,
                                         tau_out,
                                         sigma_out);
                // compute their labels.
                tau_out.toCanonicalFaceLabel(label_tau_out);
                sigma_out.toCanonicalFaceLabel(label_sigma_out);
                // get the bitmask indicating the coface traversed from tau out.
                bitmask_traversed_out = label_tau_out.to_ullong();
                bitmask_traversed_out ^= label_sigma_out.to_ullong();
                bit_idx_out = std::countr_zero(bitmask_traversed_out);
                if (tau_out.grid_coord != sigma_out.grid_coord) {
                    bit_idx_out += NDIM;
                }
                // add the output tau to faces_to_be_processed.
                tau_out.toBitset(bitset_tau_out);
                tau_to_be_processed.push_back(bitset_tau_out);
                // add the output sigma to cofaces_yet_to_be_processed.
                {
                    auto iter_out = cofaces_yet_to_be_processed.find(bitset_tau_out);
                    if (iter_out == cofaces_yet_to_be_processed.end()) {
                        const size_t fixed_coord = label_tau_out.to_ullong();
                        CofacesBitmask all_neighbors;
                        all_neighbors |= fixed_coord;
                        all_neighbors <<= NDIM;
                        all_neighbors |= fixed_coord;
                        iter_out = cofaces_yet_to_be_processed.emplace(bitset_tau_out, all_neighbors).first;
                    } else if (!iter_out->second[bit_idx_out]) {
                        throw std::runtime_error("output simplex cannot be already processed");
                    }
                    iter_out->second.reset(bit_idx_out);
                }
                // save edge: mfVert_in -> mfVert_out, coface (grid, label, vertices, the same way marching hypercubes).
                coordToBitset(sigma_out.grid_coord, bitset_sigma_coord);
                ++n_saved_cells;
                writeOutputCell(n_saved_cells,
                                bitset_sigma_coord,
                                label_sigma_out,
                                bitset_tau_in,
                                mf_vert_in,
                                bitset_tau_out,
                                mf_vert_out,
                                fout);
            }
        }
    }
    if (yet_to_be_processed.any()) {
        throw std::runtime_error("after processing all cofaces, face_in.second must be zero");
    }
}

void addFirstPoint() {
    // run pta until it gets the first valid k-simplex. add this simplex to the queue.
    CanonicalSimplex<KDIM> found_simplex;
    if (!find_first_simplex(found_simplex)) {
        throw std::runtime_error("could not find first simplex");
    }
    BitsetSimplex<KDIM> bitset_simplex;
    found_simplex.toBitset(bitset_simplex);
    tau_to_be_processed.push_back(bitset_simplex);
    // adding to the dictionary.
    {
        CanonicalFaceLabel simplex_label;
        found_simplex.toCanonicalFaceLabel(simplex_label);
        const size_t fixed_coord = simplex_label.to_ullong();
        CofacesBitmask all_neighbors;
        all_neighbors |= fixed_coord;
        all_neighbors <<= NDIM;
        all_neighbors |= fixed_coord;
        cofaces_yet_to_be_processed.emplace(bitset_simplex, all_neighbors);
    }
}

void continuation_fch(const std::string &filename) {
    FILE *fout = fopen(filename.c_str(), "wb");
    if (fout == nullptr) {
        throw std::runtime_error("ERROR OPENING FILE");
    }

    BitsetSimplex<KDIM> bitset_tau;
    CanonicalSimplex<KDIM> tau;

    std::cout << tau_to_be_processed.size() << std::endl;

    while (!tau_to_be_processed.empty()) {
        bitset_tau = tau_to_be_processed.pop_front();

        // convert the label of the k-simplex to coordinates.
        tau.fromBitset(bitset_tau);

        // calculate the cofaces.
        // then, for each coface, calculate the faces.
        createCosimplicesFromSimplex(bitset_tau, tau, fout);
    }

    // close output file.
    {
        constexpr size_t EOF_FLAG = std::numeric_limits<size_t>::max();
        fwrite(&EOF_FLAG, sizeof(size_t), 1, fout);
    }
    fclose(fout);
    fout = nullptr;

    std::cout << n_saved_cells << " cells generated" << std::endl;
}

int main(int argc, char *argv[]) {
    // specifying the input arguments.
    argparse::ArgumentParser program("fch", "", argparse::default_arguments::help);
    program.add_description(
        "This program executes the FCH.");
    program.set_usage_max_line_width(80);
    program.add_usage_newline();
    program.add_argument("-o", "--output")
        .default_value("out_fch.bin")
        .required()
        .help("specify the output .bin file.");
    program.add_argument("--verbose")
        .help("prints more information in the terminal while the program runs.")
        .default_value(false)
        .implicit_value(true);

    // reading the input arguments.
    try {
        program.parse_args(argc, argv);
    } catch (const std::exception &err) {
        std::cerr << err.what() << std::endl;
        std::cerr << program << std::endl;
        return -1;
    }

    const std::string output_path = program.get<std::string>("--output");

    if (program["--verbose"] == true) {
        std::cout << "Verbosity enabled." << std::endl;
        std::cout << "Output file: " << output_path << std::endl;
    }

    addFirstPoint();

    const auto t0 = std::chrono::high_resolution_clock::now();
    continuation_fch(output_path);
    const auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "Time to generate bin (FCH): " << std::chrono::duration<double>(t1 - t0).count()
              << std::endl;

    return 0;
}