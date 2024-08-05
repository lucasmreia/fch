///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#include <algorithm>
#include <argparse/argparse.hpp>
#include <cstdio>
#include <iostream>
#include <unordered_map>
#include <unordered_set>

#include "../common/bmi2.h"
#include "permutahedron.h"

// constexpr bool checkVertexLabelsEnumSize() {
//     // if size_t can't store the entire size of the integer, it throws an exception.
//     size_t aux = std::numeric_limits<size_t>::max();
//     aux /= (1ULL << NDIM);
//     for (size_t i = 1; i <= NDIM; ++i) {
//         if (aux < (1ULL << NDIM)) {
//             return false;
//         }
//         aux /= (1ULL << NDIM);
//     }
//     return true;
// }
//
// static_assert(checkVertexLabelsEnumSize(), "vertex label enum integer is too big");

// a simplex can have up to (NDIM + 1) vertices, and each vertex has a label of NDIM bits.
typedef std::bitset<NDIM *(NDIM + 1)> BitsetVertLabels;

struct Simplex {
    std::array<size_t, NDIM> grid_coord_canonical{};
    std::array<size_t, (KDIM + 1)> vert_labels{};
    std::array<double, NDIM> mf_vertex{};
};

void vertLabelsToBitset(const std::vector<size_t> &vert_labels, BitsetVertLabels &b) {
    if (vert_labels.size() > (NDIM + 1)) {
        throw std::runtime_error("too many vertices in a simplex");
    }
    b.reset();
    for (size_t vl : vert_labels) {
        b <<= NDIM;
        b |= vl;
    }
}

typedef std::unordered_map<BitsetVertLabels, std::unordered_set<size_t>> EdgeMap;
typedef std::unordered_map<BitsetVertLabels, std::vector<size_t>> LabelMap;

std::vector<Simplex> simplices;
std::unordered_map<size_t, std::unordered_set<size_t>> hcubes;

void createCofacesFromFace(const SFace &face, const size_t idx_simplex) {
    // print_face(face);
    SCoface coface, coface_canonical;
    size_t label_grid;
    coface.grid_coord = face.grid_coord;
    size_t pcount, idx = 0, part_a, part_b;
    P_MASK pdep_mask;
    for (size_t i = 0; i < (KDIM + 1); ++i) {  // looping through the face vectors to see which one can be split.
        pcount = std::popcount(face.e[i]);     // popcount gives the number of vectors added together.
        if (pcount <= 1) {                     // if it can't be divided, just assign it.
            coface.e[idx] = face.e[i];
            ++idx;
        } else {
            pdep_mask = pmask(face.e[i]);  // the mask for pdep is the result vector itself.
            // looping from 000...0001 to 111...1110 to split the groups of vectors.
            for (size_t j = 1; j <= ((1ULL << pcount) - 2); ++j) {
                part_a = j;  // separating in two groups, group a and group b.
                part_b = part_a ^ ((1ULL << pcount) - 1);
                coface.e[idx] = pdep(part_a, pdep_mask);
                coface.e[idx + 1] = pdep(part_b, pdep_mask);
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
                    label_grid = coordToEnum(coface_canonical.grid_coord, DOMAIN_GRID_BASIS1);
                    auto iter = hcubes.find(label_grid);
                    if (iter != hcubes.end()) {
                        iter->second.insert(idx_simplex);
                    } else {
                        hcubes[label_grid] = { idx_simplex };
                    }
                }
            }
            coface.e[idx] = face.e[i];  // now just assign this one to split the next ones.
            ++idx;
        }
    }
}

void connectEdges(const std::array<size_t, NDIM> &grid_coord,
                  const std::vector<Simplex> &mf_vertices,
                  const std::vector<size_t> &v_idx,
                  EdgeMap &edges, LabelMap &labels) {
    edges.clear();
    labels.clear();
    size_t pi, pj, count;
    size_t grid_offset_i, grid_offset_j;
    std::array<size_t, KDIM + 1> vlabels_simplex_k_i{};
    std::array<size_t, KDIM + 1> vlabels_simplex_k_j{};
    std::vector<size_t> edge_label;
    BitsetVertLabels bitset_edge_label;
    for (size_t i = 0; i < (v_idx.size() - 1); ++i) {
        // calculating the grid offset of the canonical representation for the current position.
        grid_offset_i = 0;
        for (size_t ii = 0; ii < NDIM; ++ii) {
            if (mf_vertices[v_idx[i]].grid_coord_canonical[ii] > grid_coord[ii]) {
                grid_offset_i += (1ULL << ii);
            } else if (mf_vertices[v_idx[i]].grid_coord_canonical[ii] < grid_coord[ii]) {
                grid_offset_i -= (1ULL << ii);
            }
        }
        // calculating the labels of vertices of the k-simplex at the current grid position.
        for (size_t ii = 0; ii < (KDIM + 1); ++ii) {
            vlabels_simplex_k_i[ii] = mf_vertices[v_idx[i]].vert_labels[ii] + grid_offset_i;
        }

        for (size_t j = (i + 1); j < v_idx.size(); ++j) {
            // calculating the grid offset of the canonical representation for the current position.
            grid_offset_j = 0;
            for (size_t ii = 0; ii < NDIM; ++ii) {
                if (mf_vertices[v_idx[j]].grid_coord_canonical[ii] > grid_coord[ii]) {
                    grid_offset_j += (1ULL << ii);
                } else if (mf_vertices[v_idx[j]].grid_coord_canonical[ii] < grid_coord[ii]) {
                    grid_offset_j -= (1ULL << ii);
                }
            }
            // calculating the labels of vertices of the k-simplex at the current grid position.
            for (size_t ii = 0; ii < (KDIM + 1); ++ii) {
                vlabels_simplex_k_j[ii] = mf_vertices[v_idx[j]].vert_labels[ii] + grid_offset_j;
            }

            edge_label.clear();
            pi = 0;
            pj = 0;
            count = 0;
            while ((pi < (KDIM + 1)) && (pj < (KDIM + 1))) {
                if (vlabels_simplex_k_i[pi] == vlabels_simplex_k_j[pj]) {
                    edge_label.push_back(vlabels_simplex_k_i[pi]);
                    ++pi;
                    ++pj;
                    ++count;
                } else if (vlabels_simplex_k_i[pi] > vlabels_simplex_k_j[pj]) {
                    edge_label.push_back(vlabels_simplex_k_j[pj]);
                    ++pj;
                } else {
                    edge_label.push_back(vlabels_simplex_k_i[pi]);
                    ++pi;
                }
            }
            if (count == KDIM) {
                // pi or pj will be at the end. at least one of the two will have reached the end.
                // now assign the leftovers to the label.
                while (pi < (KDIM + 1)) {
                    edge_label.push_back(vlabels_simplex_k_i[pi]);
                    ++pi;
                }
                while (pj < (KDIM + 1)) {
                    edge_label.push_back(vlabels_simplex_k_j[pj]);
                    ++pj;
                }
                if (edge_label.size() != (KDIM + 2)) {
                    throw std::runtime_error("connect_edges: wrong resulting dimension");
                }
                // add to the list of edges.
                vertLabelsToBitset(edge_label, bitset_edge_label);
                auto iter = edges.find(bitset_edge_label);
                if (iter != edges.end()) {
                    iter->second.insert(i);
                    iter->second.insert(j);
                } else {
                    edges[bitset_edge_label] = { i, j };
                    labels[bitset_edge_label] = edge_label;
                }
            }
        }
    }
}

void skeleton(const size_t k_vert,
              const std::vector<std::vector<size_t>> &label_vert,
              EdgeMap &edges, LabelMap &labels) {
    edges.clear();
    labels.clear();
    size_t pi, pj, count;
    std::vector<size_t> edge_label;
    BitsetVertLabels bitset_edge_label;
    for (size_t i = 0; i < (label_vert.size() - 1); ++i) {
        for (size_t j = (i + 1); j < label_vert.size(); ++j) {
            edge_label.clear();
            pi = 0;
            pj = 0;
            count = 0;
            while ((pi < (k_vert + 1)) && (pj < (k_vert + 1))) {
                if (label_vert[i][pi] == label_vert[j][pj]) {
                    edge_label.push_back(label_vert[i][pi]);
                    ++pi;
                    ++pj;
                    ++count;
                } else if (label_vert[i][pi] > label_vert[j][pj]) {
                    edge_label.push_back(label_vert[j][pj]);
                    ++pj;
                } else {
                    edge_label.push_back(label_vert[i][pi]);
                    ++pi;
                }
            }
            if (count == k_vert) {
                // pi or pj will be at the end. at least one of the two will have reached the end.
                // now assign the leftovers to the label.
                while (pi < (k_vert + 1)) {
                    edge_label.push_back(label_vert[i][pi]);
                    ++pi;
                }
                while (pj < (k_vert + 1)) {
                    edge_label.push_back(label_vert[j][pj]);
                    ++pj;
                }
                if (edge_label.size() != (k_vert + 2)) {
                    throw std::runtime_error("connect_edges: wrong resulting dimension");
                }
                // add to the list of edges.
                vertLabelsToBitset(edge_label, bitset_edge_label);
                auto iter = edges.find(bitset_edge_label);
                if (iter != edges.end()) {
                    iter->second.insert(i);
                    iter->second.insert(j);
                } else {
                    edges[bitset_edge_label] = { i, j };
                    labels[bitset_edge_label] = edge_label;
                }
            }
        }
    }
}

void connectSimplicesOfGridCell(const size_t g_idx, const size_t &label_grid,
                                const std::unordered_set<size_t> &s_idx,
                                const std::string &float_format,
                                FILE *fout) {
    // position of the cube in the grid.
    std::array<size_t, NDIM> grid{};
    enumToCoord(label_grid, DOMAIN_GRID_BASIS1, grid);

    // listing the vertices of the manifold.
    std::vector<Simplex> mfVertices;
    mfVertices.reserve(s_idx.size());
    for (size_t idx : s_idx) {
        mfVertices.push_back(simplices[idx]);
    }
    std::vector<bool> mfVertices_used(mfVertices.size(), false);

    // canonical basis in R^n.
    std::array<size_t, NDIM> e_n{};
    for (size_t i = 0; i < NDIM; ++i) {
        e_n[i] = 1ULL << i;
    }

    // break the n-hypercube in into n-simplices.
    // n-simplex.
    std::array<size_t, NDIM + 1> simplex_n{};

    size_t pn, pk, count;
    std::vector<std::vector<size_t>> v_idx_all;

    // variables to adjust the labels of the vertices in relation to the grid.
    size_t grid_offset;
    std::array<size_t, KDIM + 1> vlabels_simplex_k{};

    // looping through the n-simplices.
    simplex_n[0] = 0;
    simplex_n[NDIM] = (1ULL << NDIM) - 1;
    do {
        for (size_t i = 0; i < (NDIM - 1); ++i) {
            simplex_n[i + 1] = simplex_n[i] + e_n[i];
        }

        // see which k-simplices this n-simplex contains.
        // looping through the k-simplices.
        std::vector<size_t> v_idx;
        for (size_t idx = 0; idx < mfVertices.size(); ++idx) {
            // calculating the grid offset of the canonical representation for the current position.
            grid_offset = 0;
            for (size_t i = 0; i < NDIM; ++i) {
                if (mfVertices[idx].grid_coord_canonical[i] > grid[i]) {
                    grid_offset += (1ULL << i);
                } else if (mfVertices[idx].grid_coord_canonical[i] < grid[i]) {
                    grid_offset -= (1ULL << i);
                }
            }

            // calculating the labels of the vertices that make up the k-simplex at the current grid position.
            for (size_t i = 0; i < (KDIM + 1); ++i) {
                vlabels_simplex_k[i] = mfVertices[idx].vert_labels[i] + grid_offset;
            }

            // checking if the k-simplex belongs to the n-simplex.
            pn = 0;
            pk = 0;
            count = 0;
            while ((pn < (NDIM + 1)) && (pk < (KDIM + 1))) {
                if (simplex_n[pn] == vlabels_simplex_k[pk]) {
                    ++count;
                    ++pn;
                    ++pk;
                } else if (simplex_n[pn] > vlabels_simplex_k[pk]) {
                    ++pk;
                } else {
                    ++pn;
                }
            }
            if (count == (KDIM + 1)) {
                v_idx.push_back(idx);
            }
        }
        // the manifold has dimension n-k. at least n-k+1 vertices are required.
        if (v_idx.size() >= (NDIM - KDIM + 1)) {
            // marking the vertices that were used.
            for (const size_t idx : v_idx) {
                mfVertices_used[idx] = true;
            }
            // adding to the list of vertices of each simplex.
            v_idx_all.push_back(v_idx);
        }

    } while (std::next_permutation(e_n.begin(), e_n.end()));

    // remapping the vertex indices to discount the unused vertices.
    std::vector<size_t> remap_idx(mfVertices.size());
    size_t idx_new = 0;
    for (size_t i = 0; i < mfVertices.size(); ++i) {
        remap_idx[i] = -1;
        if (mfVertices_used[i]) {
            remap_idx[i] = idx_new;
            ++idx_new;
        }
    }

    // remove unused vertices from mfVertices.
    for (int64_t i = (static_cast<int64_t>(mfVertices_used.size()) - 1); i >= 0; --i) {
        if (!mfVertices_used[i]) {
            mfVertices.erase(mfVertices.begin() + i);
        }
    }

    // update the vertex indices in v_idx_all
    for (size_t i = 0; i < v_idx_all.size(); ++i) {
        for (size_t j = 0; j < v_idx_all[i].size(); ++j) {
            v_idx_all[i][j] = remap_idx[v_idx_all[i][j]];
        }
    }

    if (mfVertices.empty()) {
        return;
    }

    // saving only the region of the cut.
    // if((grid[5] < (DOMAIN_DIV[5] / 2)) || (grid[5] > (1 + DOMAIN_DIV[5] / 2))) {
    // // if(grid[5] != (DOMAIN_DIV[5] / 2)) {
    //     return;
    // }

    // save the header of the cube.
    fprintf(fout, "%3lld ", static_cast<uint64_t>(g_idx));
    for (int64_t i = 0; i < NDIM; ++i) {
        fprintf(fout, "%3lld ", static_cast<uint64_t>(grid[i]));
    }
    fprintf(fout, "\n");

    // save the number of connected components in the file (number of simplices).
    {
        const uint64_t nc = v_idx_all.size();
        fprintf(fout, "%3lld ", nc);
        fprintf(fout, "\n");
    }

    // knowing which n-simplex each k-simplex belongs to, connect the vertices into edges.
    EdgeMap edges;
    LabelMap labels;
    std::vector<std::vector<size_t>> label_vert;
    std::vector<std::vector<size_t>> conn;
    size_t k_dim;
    for (size_t i = 0; i < v_idx_all.size(); ++i) {
        // save the number of manifold vertices in the simplex in the file.
        {
            const uint64_t nv = v_idx_all[i].size();
            fprintf(fout, "%3lld ", nv);
            fprintf(fout, "\n");
        }
        // save the vertices of the simplex in the file.
        for (size_t ii = 0; ii < v_idx_all[i].size(); ++ii) {
            // saving the label of the vertices of the face that contains the manifold vertex.
            // calculating the offset in the grid of the canonical representation for the current position.
            grid_offset = 0;
            for (size_t j = 0; j < NDIM; ++j) {
                if (mfVertices[v_idx_all[i][ii]].grid_coord_canonical[j] > grid[j]) {
                    grid_offset += (1ULL << j);
                } else if (mfVertices[v_idx_all[i][ii]].grid_coord_canonical[j] < grid[j]) {
                    grid_offset -= (1ULL << j);
                }
            }
            // calculating the vertex labels of the k-simplex at the current grid position.
            for (size_t j = 0; j < (KDIM + 1); ++j) {
                const uint64_t l = mfVertices[v_idx_all[i][ii]].vert_labels[j] + grid_offset;
                fprintf(fout, "%3lld ", l);
            }
            // saving the position of the manifold vertex.
            for (size_t j = 0; j < NDIM; ++j) {
                const double v = mfVertices[v_idx_all[i][ii]].mf_vertex[j];
                fprintf(fout, float_format.c_str(), v);
            }
            fprintf(fout, "\n");
        }

        connectEdges(grid, mfVertices, v_idx_all[i], edges, labels);

        // checking if the number of vertices per edge is correct.
        // assigning the labels and connections of the edges.
        label_vert.clear();
        conn.clear();
        for (const auto &e : edges) {
            if (e.second.size() != 2) {
                throw std::runtime_error("edge has more than 2 vertices");
            }

            label_vert.push_back(labels[e.first]);
            conn.emplace_back();
            for (const auto idx : e.second) {
                conn.back().push_back(idx);
            }
        }

        // save the number of edges in the simplex in the file.
        {
            const uint64_t ne = conn.size();
            fprintf(fout, "%3lld ", ne);
            fprintf(fout, "\n");
        }

        // save the edges in the simplex in the file.
        for (size_t ii = 0; ii < conn.size(); ++ii) {
            for (size_t j = 0; j < conn[ii].size(); ++j) {
                const uint64_t e = conn[ii][j] + 1;
                fprintf(fout, "%3lld ", e);
            }
            fprintf(fout, "\n");
        }

        // knowing the edges, run the complete combinatorial skeleton.
        k_dim = KDIM + 1;
        while (k_dim < NDIM) {
            if (label_vert.size() < 2) {
                throw std::runtime_error("label_vert.size() must be >= 2");
            }

            // if(label_vert.size() > 1) {
            // running the current step of the combinatorial skeleton.
            skeleton(k_dim, label_vert, edges, labels);

            // assigning labels and connections.
            label_vert.clear();
            conn.clear();
            for (const auto &e : edges) {
                label_vert.push_back(labels[e.first]);
                conn.emplace_back();
                for (const auto idx : e.second) {
                    conn.back().push_back(idx);
                }
            }

            // saving the current step of the combinatorial skeleton in the output file.
            // saving the current number of connections of the simplex in the file.
            {
                const uint64_t ne = conn.size();
                fprintf(fout, "%3lld ", ne);
                fprintf(fout, "\n");
            }

            // save the current connections of the simplex in the file.
            for (size_t ii = 0; ii < conn.size(); ++ii) {
                for (size_t j = 0; j < conn[ii].size(); ++j) {
                    const uint64_t e = conn[ii][j] + 1;
                    fprintf(fout, "%3lld ", e);
                }
                fprintf(fout, "\n");
            }
            // } else {
            //     // saving the current step of the combinatorial skeleton in the output file.
            //     // saving the current number of connections of the simplex in the file.
            //     {
            //         const uint64_t ne = 1;
            //         fprintf(fout, "%3lld ", ne);
            //         fprintf(fout, "\n");
            //
            //         const uint64_t e = 1;
            //         fprintf(fout, "%3lld ", e);
            //         fprintf(fout, "\n");
            //     }
            // }

            ++k_dim;
        }
    }

    fprintf(fout, "\n");
}

int main(int argc, char *argv[]) {
    // specifying the input arguments.
    argparse::ArgumentParser program("skeleton_pta", "", argparse::default_arguments::help);
    program.add_description(
        "This program executes the Combinatorial Skeleton for the binary output of PTA.");
    program.set_usage_max_line_width(80);
    program.add_usage_newline();
    program.add_argument("-i", "--input")
        .default_value("out_pta.bin")
        .required()
        .help("specify the input .bin file.");
    program.add_argument("-o", "--output")
        .default_value("out_pta.pol")
        .required()
        .help("specify the output .pol file.");
    program.add_argument("-ff", "--float_format")
        .default_value("15.8")
        .help("specify format of the float when printing to the .pol file.");
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

    const std::string input_path = program.get<std::string>("--input");
    const std::string output_path = program.get<std::string>("--output");
    const std::string float_format = "%" + program.get<std::string>("--float_format") + "f ";

    if (program["--verbose"] == true) {
        std::cout << "Verbosity enabled." << std::endl;
        std::cout << "Input file: " << input_path << std::endl;
        std::cout << "Output file: " << output_path << std::endl;
        std::cout << "Float format: " << float_format << std::endl;
    }

    FILE *fin = fopen(input_path.c_str(), "rb");
    if (fin == nullptr) {
        std::cout << "input file not found" << std::endl;
        return -1;
    }

    FILE *fout = fopen(output_path.c_str(), "w");

    size_t count;
    BitsetPermutahedronSFace bitset_face;
    SFace face;
    Simplex s;
    size_t simplex_idx, g_idx;

    uint64_t n_saved_simplices = 0;

    fprintf(fout, "%3lu %3lu\n", NDIM, KDIM);
    for (int64_t i = 0; i < NDIM; ++i) {
        fprintf(fout, "%3lu ", DOMAIN_DIV[i]);
    }
    fprintf(fout, "\n\n");

    while (true) {
        // reading the cell number.
        fread(&count, sizeof(size_t), 1, fin);

        // checking if it is eof.
        if (count == std::numeric_limits<uint64_t>::max()) {
            break;
        }

        ++n_saved_simplices;

        // reading the label.
        bitset_face.reset();
        bitset_face.read(fin);

        face.fromBitset(bitset_face);

        // setting the grid of the simplex (in canonical notation).
        for (size_t j = 0; j < NDIM; ++j) {
            s.grid_coord_canonical[j] = face.grid_coord[j];
        }

        // determining the vertices of the simplex.
        s.vert_labels[0] = 0;
        for (size_t i = 0; i < KDIM; ++i) {  // notação canonica, referencia é o vértice 0.
            s.vert_labels[i + 1] = s.vert_labels[i] + face.e[i];
        }

        // reading the vertex of the manifold.
        fread(s.mf_vertex.data(), sizeof(double), NDIM, fin);

        // adding simplex to the list.
        simplex_idx = simplices.size();
        simplices.push_back(s);

        // adding simplex to the grid positions.
        createCofacesFromFace(face, simplex_idx);
    }

    // run the skeleton for each cube now.
    g_idx = 1;
    for (const auto &cube : hcubes) {
        connectSimplicesOfGridCell(g_idx, cube.first, cube.second, float_format, fout);
        ++g_idx;
    }

    fprintf(fout, "-1\n");

    fclose(fin);
    fclose(fout);

    std::cout << n_saved_simplices << " cells generated" << std::endl;

    return 0;
}
