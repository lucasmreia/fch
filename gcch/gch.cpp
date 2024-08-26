///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#include "gch.h"

#include <iostream>
#include <set>
#include <unordered_map>
#include <unordered_set>

template <size_t DIM>
void vertLabelsToInductive(const std::array<size_t, DIM + 1> &vert_labels,
                           std::array<UINT_EIDX, DIM> &e_idx) {
    for (size_t i = 0; i < DIM; ++i) {
        e_idx[i] = std::countr_zero(vert_labels[i + 1] - vert_labels[i]);
    }
}

template <size_t DIM>
void inductiveToVertLabels(const size_t first_vert_label,
                           const std::array<UINT_EIDX, DIM> &e_idx,
                           std::array<size_t, DIM + 1> &vert_labels) {
    vert_labels[0] = first_vert_label;
    for (size_t i = 1; i <= DIM; ++i) {
        vert_labels[i] = vert_labels[i - 1] | (1ULL << e_idx[i - 1]);
    }
}

// based on the hypercube face label and the direction of the coface, inserts a vertex in the dim-simplex
// to make it a (dim+1)-simplex.
// returns the index, within the (dim+1)-simplex, of the new vertex added to the dim-simplex (either 0 or DIM+1).
// within the (k+1)-hypercube there is only one (k+1)-simplex that contains the k-simplex within the hypercube k-face.
template <size_t DIM>
size_t elevateDimension(const Label2N &in_hface_label, const std::array<UINT_EIDX, DIM> &in_e_idx,
                        const size_t unset_bit_idx,
                        Label2N &out_hface_label, std::array<UINT_EIDX, DIM + 1> &out_e_idx) {
    out_hface_label = in_hface_label;
    out_hface_label.reset(unset_bit_idx);
    const size_t unset_coord = unset_bit_idx % NDIM;
    if (unset_bit_idx < NDIM) {
        std::memcpy(out_e_idx.data(), in_e_idx.data(), DIM * sizeof(UINT_EIDX));
        out_e_idx[DIM] = unset_coord;
        return (DIM + 1);  // last vertex of the (k+1)-simplex.
        // NOTE: IT'S THE POSITION OF THE LAST VERTEX, NOT THE LAST e_idx, SO IT'S +1.
    } else {
        out_e_idx[0] = unset_coord;
        std::memcpy(&out_e_idx[1], in_e_idx.data(), DIM * sizeof(UINT_EIDX));
        return 0;  // first vertex of the (k+1)-simplex.
    }
}

// Get the index of the column that has to be replaced.
// I already pass lu_decomp and lambda to reuse those values.
template <size_t DIM>
int64_t getReplacedColumnIdx(const Eigen::PartialPivLU<Eigen::Matrix<double, DIM, DIM>> &lu_decomp,
                             const Eigen::Vector<double, DIM + 1> &lambda,
                             const Eigen::Vector<double, DIM> &b_mu) {
    // solve the system to calculate mu.
    Eigen::Vector<double, DIM + 1> mu = Eigen::Vector<double, DIM + 1>::Zero();
    mu.tail(DIM) = lu_decomp.solve(b_mu);
    mu(0) = 1. - mu.sum();

    // calcualte min{lambda/mu} and the index of the vertex that must be replaced.
    double min_lambda_mu = 0;
    int64_t idx_k = -1;
    for (int64_t ii = 0; ii < (DIM + 1); ++ii) {
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

// vertices, Vert_simplex_k_1 and extra_vertex_idx are progressively being modified as you pivot.
template <size_t DIM>
void traverseHCofaceHypercube(const std::vector<std::array<double, NDIM>> &vert_hcube,
                              std::array<size_t, DIM + 2> &vert_labels_hc,
                              size_t &extra_vertex_idx) {
    // keep solving the system and pivoting until you find the output simplex.
    // create the system.
    Eigen::Matrix<double, DIM, DIM> A{};
    Eigen::Matrix<double, NDIM, DIM + 1> V{};
    std::array<size_t, DIM + 1> columns{};  // label of the vertex assigned to each column.

    Eigen::Vector<double, DIM> b_lambda{};

    // calculating the positions of the vertices and evalVertices.
    // creating the matrices of the A*x=b system.
    for (int64_t j = 0, v_idx = 0; j < (DIM + 1); ++j, ++v_idx) {
        if (v_idx == extra_vertex_idx) {
            ++v_idx;
        }
        columns[j] = v_idx;
        std::memcpy(&V(0, j), vert_hcube[vert_labels_hc[v_idx]].data(), NDIM * sizeof(double));
        if (j == 0) {
            func(vert_hcube[vert_labels_hc[v_idx]].data(), b_lambda.data());
            b_lambda *= -1;
        } else {
            func(vert_hcube[vert_labels_hc[v_idx]].data(), &A(0, j - 1));
            A.col(j - 1) += b_lambda;
        }
    }

    // assigning b_mu.
    Eigen::Vector<double, DIM> b_mu{};
    func(vert_hcube[vert_labels_hc[extra_vertex_idx]].data(), b_mu.data());
    b_mu += b_lambda;

    // variables to solve the system.
    Eigen::Vector<double, DIM + 1> lambda{};
    Eigen::PartialPivLU<Eigen::Matrix<double, DIM, DIM>> lu_decomp;

    while (true) {
        // solve the system to calculate lambda.
        lu_decomp = A.partialPivLu();
        lambda(0) = 0;
        lambda.tail(KDIM) = lu_decomp.solve(b_lambda);
        lambda(0) = 1. - lambda.sum();

        if (!(lambda.minCoeff() >= LAMBDA_ZERO_THRESHOLD)) {
            throw std::runtime_error("door-in,door-out could not find output k-simplex");
        }

        const int64_t replaced_column_idx = getReplacedColumnIdx<DIM>(lu_decomp, lambda, b_mu);

        // replace the column with the extra_vertex_idx.
        std::swap(columns[replaced_column_idx], extra_vertex_idx);

        // pivoting the k+1 simplex on extra_vertex_idx to find the new extra vertex.
        // however, if it wants to pivot on the first or the last vertex of the (k+1)-simplex, if means the manifold
        // has left the (k+1)-hypercube, so the traversal ends.
        if ((extra_vertex_idx == 0) || (extra_vertex_idx == (DIM + 1))) {
            // end of the traversal inside the coface.
            // vert_labels_hc[columns] already has the vertex labels used by the output k-simplex.
            // only return.
            return;
        }

        // updating the column with the new vertex.
        std::memcpy(&V(0, replaced_column_idx), vert_hcube[vert_labels_hc[columns[replaced_column_idx]]].data(),
                    NDIM * sizeof(double));

        // updating the system of equations.
        if (replaced_column_idx > 0) {
            // if it is not the first vertex, just update the replaced column.
            func(vert_hcube[vert_labels_hc[columns[replaced_column_idx]]].data(), &A(0, replaced_column_idx - 1));
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

        // pivoting the k+1 simplex on extra_vertex_idx to find the new extra vertex.
        const size_t previous_v_idx = extra_vertex_idx - 1;
        const size_t next_v_idx = extra_vertex_idx + 1;
        const size_t dif_e = vert_labels_hc[next_v_idx] - vert_labels_hc[previous_v_idx];
        const size_t pivoted_v_label = vert_labels_hc[extra_vertex_idx] ^ dif_e;

        // assign a new pivoted vertex.
        vert_labels_hc[extra_vertex_idx] = pivoted_v_label;

        // updating b_mu with the new pivoted vertex.
        func(vert_hcube[vert_labels_hc[extra_vertex_idx]].data(), b_mu.data());
        b_mu += b_lambda;
    }
}

// loops through all k-simplices that make a hypercube k-face, and checks each simplex for intersection with the manifold.
void genVerticesFromFace(const std::array<bool, NDIM> &hface_fixed_coord,
                         const Label2N &hface_label,
                         const std::vector<std::array<double, NDIM>> &vert_hcube,
                         std::vector<std::array<double, NDIM>> &mf_vertices,
                         std::vector<std::array<size_t, KDIM + 1>> &mf_vert_labels,
                         std::vector<Label2N> &mf_hface_labels) {
    // generating a mapping with the positions of the free coordinates.
    // it is equivalent to the inductive representation they mention in the article.
    std::array<size_t, KDIM> free_coord_map{};
    size_t count = 0;
    for (size_t pos = 0; pos < NDIM; ++pos) {
        if (!hface_fixed_coord[pos]) {
            free_coord_map[count++] = pos;
        }
    }

#ifdef PRINT_DEBUG
    std::cout << "free coord: ";
    for (size_t j = 0; j < KDIM; ++j) {
        std::cout << ((int)free_coord_map[j]) << " ";
    }
    std::cout << std::endl;
#endif

    constexpr size_t mask_least_significant{ ((1ULL << NDIM) - 1ULL) };
    Label2N aux;

    // generating the vertices of the k-simplices that make up the hypercube k-face.
    // makes permutations of the inductive notation.
    // the generated values are the indices of the vertices already in the hypercube. now you can just access the hypercube directly.
    std::array<size_t, KDIM + 1> vertices_idx{};
    // taking the n most significant bits of the face label (equivalent to the index/label of the first vertex).
    aux = hface_label;
    aux >>= NDIM;
    vertices_idx[0] = aux.to_ullong();
    // taking the opposite of the n least significant bits of the face's label (equivalent to the index/label of the last vertex).
    aux = hface_label;
    aux &= mask_least_significant;
    aux ^= mask_least_significant;
    vertices_idx[KDIM] = aux.to_ullong();

    // initializing matrices.
    std::array<std::array<double, NDIM>, KDIM + 1> vert_simplex_k{};
    vert_simplex_k[0] = vert_hcube[vertices_idx[0]];
    vert_simplex_k[KDIM] = vert_hcube[vertices_idx[KDIM]];
    std::array<double, NDIM> mf_vertex{};

    // permute simplices.
    do {
        for (size_t j = 1; j < KDIM; ++j) {
            vertices_idx[j] = vertices_idx[j - 1] | (1ULL << free_coord_map[j - 1]);
            vert_simplex_k[j] = vert_hcube[vertices_idx[j]];
        }

        // indices of the generated vertices. now you can just use them.
        // the vertices can be acessed directly in the hypercube.

#ifdef PRINT_DEBUG
        std::cout << "vertices_idx: ";
        for (size_t j = 0; j < (KDIM + 1); ++j) {
            std::cout << vertices_idx[j] << " ";
        }
        std::cout << std::endl;
#endif

        const bool contains_manifold = solve(vert_simplex_k, mf_vertex);

        if (contains_manifold) {
            // storing the approximated point where the manifold intersects with the k-simplex.
            mf_vertices.emplace_back(mf_vertex);

            // storing the indices/labels of the vertices that make up the k-simplex.
            mf_vert_labels.emplace_back(vertices_idx);

            // storing the label of the hypercube k-face that generated the vertex.
            mf_hface_labels.emplace_back(hface_label);
        }

    } while (std::next_permutation(free_coord_map.begin(), free_coord_map.end()));
}

// loops through all k-faces of an n-hypercube.
void genManifoldVertices(const std::vector<std::array<double, NDIM>> &vert_hcube,
                         std::vector<std::array<double, NDIM>> &mf_vertices,
                         std::vector<std::array<size_t, KDIM + 1>> &mf_vert_labels,
                         std::vector<Label2N> &mf_hface_labels) {
    std::array<size_t, MDIM> idx{};

    for (size_t i = 0; i < MDIM; ++i) {
        idx[i] = i;
    }

    size_t count = 0;

    Label2N hface_label;
    std::array<bool, NDIM> hface_fixed_coord{};

    int64_t to_be_increased;
    bool reset;

    while (idx[0] <= (NDIM - MDIM)) {
        ++count;

        hface_fixed_coord.fill(false);
        for (size_t j = 0; j < MDIM; ++j) {
#ifdef PRINT_DEBUG
            std::cout << static_cast<int>(idx[j]) << " ";
#endif
            hface_fixed_coord[idx[j]] = true;
        }
#ifdef PRINT_DEBUG
        std::cout << std::endl;
        for (size_t j = 0; j < NDIM; ++j) {
            std::cout << static_cast<int>(hface_fixed_coord[j]) << " ";
        }
        std::cout << std::endl;
#endif

        for (size_t i = 0; i < (1 << MDIM); ++i) {
            hface_label.reset();
            for (size_t j = 0; j < MDIM; ++j) {
                hface_label.set(idx[j] + NDIM * ((i >> j) & 1));
            }
#ifdef PRINT_DEBUG
            for (size_t j = 0; j < (2 * NDIM); ++j) {
                std::cout << static_cast<int>(hface_label[j]) << " ";
            }
            std::cout << std::endl;
#endif
            genVerticesFromFace(hface_fixed_coord,
                                hface_label,
                                vert_hcube,
                                mf_vertices,
                                mf_vert_labels,
                                mf_hface_labels);
        }

        to_be_increased = MDIM - 1;
        reset = false;
        while ((to_be_increased >= 0) && (++idx[to_be_increased] > (NDIM - MDIM + to_be_increased))) {
            --to_be_increased;
            reset = true;
        }
        if (reset && (to_be_increased >= 0)) {
            for (size_t i = to_be_increased + 1; i < MDIM; ++i) {
                idx[i] = idx[i - 1] + 1;
            }
        }
    }

#ifdef PRINT_DEBUG
    std::cout << count << std::endl;
#endif
}

// connects structures of dimension dim to form structures of dimension dim+1.
// based on the adjacency rules of the combinatorial skeleton.
void genEdges(const size_t dim,
              const std::vector<Label2N> &hface_labels,
              std::vector<std::vector<size_t>> &edge_connections,
              std::vector<Label2N> &edge_hface_labels) {
    std::unordered_map<Label2N, std::set<size_t>> edge_map;
    // comparing one component to all the others.
    for (size_t i = 0; i < (hface_labels.size() - 1); ++i) {
        for (size_t j = (i + 1); j < hface_labels.size(); ++j) {
            // making the intersection between components.
            Label2N edge = hface_labels[i];
            edge &= hface_labels[j];
            if (edge.count() == (NDIM - dim)) {
                // if the intersection unlocks only 1 dimension, then I mark that
                // these 2 components belong to the same coface of 1 dimension higher.
                edge_map[edge].insert(i);
                edge_map[edge].insert(j);
            } else if (edge.count() == (NDIM - dim + 1)) {
                // if the intersection maintains the dimension, then there are two components on the
                // same face of same dimension.
                // I have to find all the cofaces of 1 dimension above that these two components can belong to.
                // I elevate the dimension of this face by 1 and mark the two components as belonging to
                // the resulting coface.
                // I have to do this in all the dimensions in which I can go up.
                // going through all the dimensions to see which one I can elevate.
                for (size_t m = 0; m < (2 * NDIM); ++m) {
                    if (edge[m]) {
                        Label2N edge_aux = edge;
                        edge_aux.reset(m);
                        edge_map[edge_aux].insert(i);
                        edge_map[edge_aux].insert(j);
                    }
                }
            }
        }
    }
    for (const auto &edge_set : edge_map) {
        edge_hface_labels.emplace_back(edge_set.first);
        std::vector<size_t> edges;
        edges.reserve(edge_set.second.size());
        for (const size_t edge : edge_set.second) {
            edges.emplace_back(edge);
        }
        edge_connections.emplace_back(edges);
    }
}

bool gch(const std::vector<std::array<double, NDIM>> &vert_hcube,
         HypercubeApprox &approx, Label2N &neighbor_cells) {
    approx.reset();
    neighbor_cells.reset();

    // Calculating the vertices on the faces of the hypercube.
    std::vector<Label2N> mf_hface_labels_aux;
    genManifoldVertices(vert_hcube,
                        approx.vertices,
                        approx.vertex_labels,
                        mf_hface_labels_aux);

    // Checking if any vertices have been found.
    if (approx.vertex_labels.empty()) {
        return false;
    }

    // Setting the neighbors of the current cube.
    for (size_t i = 0; i < mf_hface_labels_aux.size(); ++i) {
        neighbor_cells |= mf_hface_labels_aux[i];
    }

    // Generating the edges that connect the vertices.
    std::vector<std::vector<size_t>> edge_connections_aux;
    std::vector<Label2N> edge_hface_labels_aux;
    genEdges((KDIM + 1), mf_hface_labels_aux, edge_connections_aux, edge_hface_labels_aux);

    // Checking if any edge found.
    if (edge_hface_labels_aux.empty()) {
        return false;
    }

    // check for ambiguity.
    for (size_t i = 0; i < edge_connections_aux.size(); ++i) {
        if ((edge_connections_aux[i].size() % 2) || (edge_connections_aux[i].size() < 2)) {
            throw std::runtime_error("Hface of dim k+1 must contain positive even number of hfaces of dim k");
        } else if (edge_connections_aux[i].size() == 2) {
            // Add edge without ambiguity in approx.
            approx.edge_connections.emplace_back(std::array<size_t, 2>{ edge_connections_aux[i][0],
                                                                        edge_connections_aux[i][1] });
            approx.edge_hface_labels.emplace_back(edge_hface_labels_aux[i]);
        } else {
            // Edge with ambiguity.
            // taking the vertices and the hypercube coface that present ambiguity.
            std::set<size_t> connected_vertices{ edge_connections_aux[i].begin(), edge_connections_aux[i].end() };
            const Label2N &hcoface_label{ edge_hface_labels_aux[i] };

            std::vector<std::array<size_t, 2>> new_edge_connections;
            std::array<UINT_EIDX, KDIM> hface_e_idx_in{};
            std::array<UINT_EIDX, KDIM + 1> hcoface_e_idx{};
            std::array<size_t, KDIM + 2> hcoface_vertex_labels{};
            std::array<size_t, KDIM + 1> hface_vertex_labels_out{};
            Label2N label_aux;

            // traversing through the hypercube coface to connect all vertices within this coface.
            // making connections as long as there are unconnected vertices.
            while (!connected_vertices.empty()) {
                // taking the vertex to start the traversal.
                const size_t v_idx_in = *connected_vertices.begin();
                connected_vertices.erase(connected_vertices.begin());

                // creating the representation of the k-simplex that generated this vertex, and its hypercube k-face.
                const std::array<size_t, KDIM + 1> &hface_vertex_labels_in = approx.vertex_labels[v_idx_in];
                const Label2N &hface_label_in = mf_hface_labels_aux[v_idx_in];
                vertLabelsToInductive(hface_vertex_labels_in, hface_e_idx_in);

                // identifying the extra dimension of the hypercube coface.
                label_aux = hcoface_label;
                label_aux ^= hface_label_in;
                const size_t unset_bit_idx = std::countr_zero(label_aux.to_ullong());

                // elevating the dimension of the k-simplex to generate the simplex coface.
                // identifying the extra vertex that forms the coface.
                size_t extra_vertex = elevateDimension(hface_label_in, hface_e_idx_in, unset_bit_idx, label_aux,
                                                       hcoface_e_idx);
                if (hcoface_label != label_aux) {
                    throw std::runtime_error("hcoface_label != hcoface_label_aux");
                }
                inductiveToVertLabels(hcoface_label.getMsb(), hcoface_e_idx, hcoface_vertex_labels);

                // traversing inside the hypercube coface.
                traverseHCofaceHypercube<KDIM>(vert_hcube, hcoface_vertex_labels, extra_vertex);

                // creating the output simplex representation.
                for (size_t j0 = 0, j1 = 0; j0 < (KDIM + 1); ++j0, ++j1) {
                    if (j0 == extra_vertex) {
                        ++j1;
                    }
                    hface_vertex_labels_out[j0] = hcoface_vertex_labels[j1];
                }

                // identifying the output simplex in connected_vertices.
                bool found = false;
                for (auto iter = connected_vertices.begin(); iter != connected_vertices.end(); ++iter) {
                    if (hface_vertex_labels_out == approx.vertex_labels[*iter]) {
                        // inserting new edge in new_edge_connections and removing index from connected_vertices.
                        new_edge_connections.emplace_back(std::array<size_t, 2>{ v_idx_in, *iter });
                        connected_vertices.erase(iter);
                        found = true;
                        break;
                    }
                }
                if (!found) {
                    throw std::runtime_error("out simplex not found");
                }
            }

            // Add new edges without ambiguity in approx.
            for (size_t j = 0; j < new_edge_connections.size(); ++j) {
                approx.edge_connections.emplace_back(std::array<size_t, 2>{ new_edge_connections[j][0],
                                                                            new_edge_connections[j][1] });
                approx.edge_hface_labels.emplace_back(hcoface_label);
            }
        }
    }

    // returning a valid result.
    return true;
}

void writeOutputHypercube(const size_t g,
                          const BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH> &bitset_coord,
                          const HypercubeApprox &approx, FILE *fout) {
    // saving the hypercube number.
    fwrite(&g, sizeof(size_t), 1, fout);
    // saving the hypercube position.
    bitset_coord.write(fout);
    // saving the number of vertices found in the hypercube k-faces.
    {
        const size_t nv = approx.vertex_labels.size();
        fwrite(&nv, sizeof(size_t), 1, fout);
    }
    // saving the vertices of the approximation.
    for (size_t i = 0; i < approx.vertex_labels.size(); ++i) {
        // index/label of the vertices of the k-simplex that generated the approximation vertex.
        fwrite(approx.vertex_labels[i].data(), sizeof(size_t), (KDIM + 1), fout);
        // approximation vertex position.
        fwrite(approx.vertices[i].data(), sizeof(double), NDIM, fout);
    }
    // saving the number of hypercube (k+1)-faces (the cofaces).
    {
        const size_t nf = approx.edge_hface_labels.size();
        fwrite(&nf, sizeof(size_t), 1, fout);
    }
    // saving the cofaces and the connections inside them (edges).
    for (size_t i = 0; i < approx.edge_hface_labels.size(); ++i) {
        approx.edge_hface_labels[i].write(fout);
        fwrite(approx.edge_connections[i].data(), sizeof(size_t), 2, fout);
    }
}

void readOutputHypercube(FILE *fin,
                         BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH> &bitset_coord,
                         HypercubeApprox &approx) {
    approx.reset();
    // reading the coord of the hypercube.
    bitset_coord.read(fin);
    // reading the number of vertices found in the hypercube k-faces.
    size_t nv;
    fread(&nv, sizeof(size_t), 1, fin);
    // reading the vertices of the approximation.
    approx.vertex_labels.resize(nv);
    approx.vertices.resize(nv);
    for (size_t i = 0; i < nv; ++i) {
        // index/label of the vertices of the k-simplex that generated the approximation vertex.
        fread(approx.vertex_labels[i].data(), sizeof(size_t), (KDIM + 1), fin);
        // approximation vertex position.
        fread(approx.vertices[i].data(), sizeof(double), NDIM, fin);
    }
    // reading the number of hypercube (k+1)-faces (the cofaces).
    size_t nf;
    fread(&nf, sizeof(size_t), 1, fin);
    // reading the cofaces and the connections inside them (edges).
    approx.edge_hface_labels.resize(nf);
    approx.edge_connections.resize(nf);
    for (size_t i = 0; i < nf; ++i) {
        approx.edge_hface_labels[i].read(fin);
        fread(approx.edge_connections[i].data(), sizeof(size_t), 2, fin);
    }
}
