///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#include <argparse/argparse.hpp>
#include <cstdio>
#include <deque>
#include <iostream>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>

#include "../common/bmi2.h"
#include "fch.h"

struct Edge {
    std::array<size_t, NDIM> hcoface_grid_coord{};
    Label2N hcoface_label;
    std::array<size_t, NDIM> hface_grid_coord_in{};
    std::array<size_t, KDIM + 1> simplex_vert_label_in{};
    std::array<double, NDIM> mf_vert_in{};
    std::array<size_t, NDIM> hface_grid_coord_out{};
    std::array<size_t, KDIM + 1> simplex_vert_label_out{};
    std::array<double, NDIM> mf_vert_out{};
};

std::vector<Edge> edges;                                 // list of all calculated simplices.
std::unordered_map<size_t, std::vector<size_t>> hcubes;  // marking the indices of the simplices that belong to each cube.

int main(int argc, char *argv[]) {
    // specifying the input arguments.
    argparse::ArgumentParser program("fch_to_hypercubes", "", argparse::default_arguments::help);
    program.add_description(
        "This program converts the output of the FCH to the output hypercube format of the GCCH.");
    program.set_usage_max_line_width(80);
    program.add_usage_newline();
    program.add_argument("-i", "--input")
        .default_value("out_fch.bin")
        .required()
        .help("specify the input .bin file.");
    program.add_argument("-o", "--output")
        .default_value("out_fch_hypercubes.bin")
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

    const std::string input_path = program.get<std::string>("--input");
    const std::string output_path = program.get<std::string>("--output");

    if (program["--verbose"] == true) {
        std::cout << "Verbosity enabled." << std::endl;
        std::cout << "Input file: " << input_path << std::endl;
        std::cout << "Output file: " << output_path << std::endl;
    }

    FILE *fin = fopen(input_path.c_str(), "rb");

    if (fin == nullptr) {
        std::cout << "input file not found" << std::endl;
        return -1;
    }

    FILE *fout = fopen(output_path.c_str(), "wb");

    constexpr size_t EOF_FLAG = std::numeric_limits<size_t>::max();

    size_t n_read_cells = 0;

    while (true) {
        // reading the hypercube position.
        size_t hcoface_g, hface_g_in, hface_g_out;
        fread(&hcoface_g, sizeof(uint64_t), 1, fin);

        // checking if it is eof.
        if (hcoface_g == EOF_FLAG) {
            break;
        }

        ++n_read_cells;

        Edge edge;
        enumToCoord(hcoface_g, DOMAIN_GRID_BASIS, edge.hcoface_grid_coord);

        // reading the other parameters of the edge.
        readOutputCell(fin,
                       edge.hcoface_label,
                       hface_g_in,
                       edge.simplex_vert_label_in,
                       edge.mf_vert_in,
                       hface_g_out,
                       edge.simplex_vert_label_out,
                       edge.mf_vert_out);

        enumToCoord(hface_g_in, DOMAIN_GRID_BASIS, edge.hface_grid_coord_in);
        enumToCoord(hface_g_out, DOMAIN_GRID_BASIS, edge.hface_grid_coord_out);

        // adding to the vector of cofaces and taking the index.
        const size_t idx = edges.size();
        edges.emplace_back(edge);

        // assigning idx to all the hypercubes to which the coface belongs.
        size_t fixed_coord, free_coord;
        edge.hcoface_label.getFixedAndFreeCoord(fixed_coord, free_coord);
        // {
        //     std::array<size_t, NDIM> grid_coord_aux = edge.hcoface_grid_coord;
        //     size_t grid_coord_enm = coordToEnum(grid_coord_aux, DOMAIN_GRID_BASIS);
        //     auto iter = hcubes.find(grid_coord_enm);
        //     if (iter == hcubes.end()) {
        //         hcubes[grid_coord_enm] = {idx};
        //     } else {
        //         iter->second.push_back(idx);
        //     }
        // }
        const size_t n_fixed_coord = std::popcount(fixed_coord);
        const P_MASK pdep_mask = pmask(fixed_coord);
        for (size_t i = 0; i < (1ULL << n_fixed_coord); ++i) {
            const size_t neighbor = pdep(i, pdep_mask);
            std::array<size_t, NDIM> grid_coord_aux = edge.hcoface_grid_coord;
            bool skip = false;
            for (size_t j = 0; j < NDIM; ++j) {
                if (neighbor & (1ULL << j)) {
                    if (edge.hcoface_label[j]) {
                        if (grid_coord_aux[j] <= 0) {
                            skip = true;
                            break;
                        }
                        grid_coord_aux[j] -= 1;
                    } else if (edge.hcoface_label[j + NDIM]) {
                        if (grid_coord_aux[j] >= DOMAIN_DIV[j] - 1) {
                            skip = true;
                            break;
                        }
                        grid_coord_aux[j] += 1;
                    } else {
                        throw std::runtime_error("fixed coord is not fixed");
                    }
                }
            }
            if (skip) {
                continue;
            }
            size_t grid_coord_enm = coordToEnum(grid_coord_aux, DOMAIN_GRID_BASIS);
            auto iter = hcubes.find(grid_coord_enm);
            if (iter == hcubes.end()) {
                hcubes[grid_coord_enm] = { idx };
            } else {
                iter->second.push_back(idx);
            }
        }
    }

    size_t n_saved_hypercubes = 0;

    // loop through all the hypercubes, assemble the HyperCubeApprox for each cube and save it in another binary file.
    // then I run the combinatorial skeleton in this other binary file to assemble the pol.
    HypercubeApprox approx;
    for (const auto &hcube : hcubes) {
        const size_t g = hcube.first;
        const std::vector<size_t> &edges_idx = hcube.second;
        std::array<size_t, NDIM> grid_coord{};
        enumToCoord(g, DOMAIN_GRID_BASIS, grid_coord);
        approx.reset();
        // vertices on the faces of the hypercube.
        std::array<size_t, KDIM + 1> vert_label{};
        std::map<std::array<size_t, KDIM + 1>, size_t> vertex_map;
        size_t vertex_counter = 0, v_in, v_out;
        for (const size_t idx : edges_idx) {
            // adding the vert in.
            vert_label = edges[idx].simplex_vert_label_in;
            // checking if the cosimplex was specified in the same hypercube.
            for (size_t i = 0; i < NDIM; ++i) {
                if (edges[idx].hface_grid_coord_in[i] != grid_coord[i]) {
                    // adjusting the simplex's vert_labels.
                    for (size_t j = 0; j < (KDIM + 1); ++j) {
                        vert_label[j] ^= (1ULL << i);
                    }
                }
            }
            auto iter_in = vertex_map.find(vert_label);
            if (iter_in == vertex_map.end()) {
                v_in = vertex_counter;
                vertex_map.insert({ vert_label, vertex_counter });
                ++vertex_counter;
                approx.vertex_labels.push_back(vert_label);
                approx.vertices.push_back(edges[idx].mf_vert_in);
            } else {
                v_in = iter_in->second;
            }

            // adding the vert out.
            vert_label = edges[idx].simplex_vert_label_out;
            // checking if the cosimplex was specified in the same hypercube.
            for (size_t i = 0; i < NDIM; ++i) {
                if (edges[idx].hface_grid_coord_out[i] != grid_coord[i]) {
                    // adjusting the simplex's vert_labels.
                    for (size_t j = 0; j < (KDIM + 1); ++j) {
                        vert_label[j] ^= (1ULL << i);
                    }
                }
            }
            auto iter_out = vertex_map.find(vert_label);
            if (iter_out == vertex_map.end()) {
                v_out = vertex_counter;
                vertex_map.insert({ vert_label, vertex_counter });
                ++vertex_counter;
                approx.vertex_labels.push_back(vert_label);
                approx.vertices.push_back(edges[idx].mf_vert_out);
            } else {
                v_out = iter_out->second;
            }

            // determining the label of the hcoface.
            Label2N hcoface_label = edges[idx].hcoface_label;
            // checking if the cosimplex was specified in the same hypercube.
            for (size_t i = 0; i < NDIM; ++i) {
                if (edges[idx].hcoface_grid_coord[i] != grid_coord[i]) {
                    // error check.
                    if (!hcoface_label[i] && !hcoface_label[i + NDIM]) {
                        throw std::runtime_error("hcoface_label[i] or hcoface_label[i+NDIM] must be true");
                    }
                    // adjusting the simplex's vert_labels.
                    hcoface_label.flip(i);
                    hcoface_label.flip(i + NDIM);
                }
            }

            // adding the edge.
            approx.edge_hface_labels.push_back(hcoface_label);
            approx.edge_connections.push_back({ v_in, v_out });
        }

        // saving to the output file.
        if (approx.vertex_labels.size() >= (MDIM + 1)) {
            ++n_saved_hypercubes;
            writeOutputHypercube(g, approx, fout);
        }
    }

    // close output file.
    fwrite(&EOF_FLAG, sizeof(size_t), 1, fout);
    fclose(fin);
    fclose(fout);

    std::cout << n_read_cells << " cells read" << std::endl;
    std::cout << n_saved_hypercubes << " hypercubes generated" << std::endl;

    return 0;
}
