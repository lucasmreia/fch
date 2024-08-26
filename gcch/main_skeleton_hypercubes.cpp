///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#include <argparse/argparse.hpp>
#include <cstdio>
#include <deque>
#include <filesystem>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../common/definitions.h"
#include "gch.h"

std::vector<HypercubeApprox> separateConnectedComponents(const HypercubeApprox &approx) {
    std::vector<HypercubeApprox> connected_components;
    // creating a list of links.
    std::vector<std::vector<size_t>> links(approx.vertices.size());
    for (size_t i = 0; i < approx.edge_connections.size(); ++i) {
        links[approx.edge_connections[i][0]].push_back(approx.edge_connections[i][1]);
        links[approx.edge_connections[i][1]].push_back(approx.edge_connections[i][0]);
    }
    // creating markers for visited nodes.
    std::vector<bool> visited(approx.vertices.size(), false);
    while (true) {
        // finding the first one not visited.
        size_t first_idx;
        for (first_idx = 0; first_idx < approx.vertices.size(); ++first_idx) {
            if (!visited[first_idx]) {
                break;
            }
        }
        if (first_idx == approx.vertices.size()) {
            // all have been visited.
            break;
        }
        // executing the search and retrieving the visited vertices.
        std::deque<size_t> to_be_visited;
        to_be_visited.push_back(first_idx);
        std::vector<size_t> visited_vertices;
        while (!to_be_visited.empty()) {
            const size_t idx = to_be_visited.front();
            to_be_visited.pop_front();
            if (!visited[idx]) {
                visited[idx] = true;
                visited_vertices.push_back(idx);
                for (size_t i = 0; i < links[idx].size(); ++i) {
                    to_be_visited.push_back(links[idx][i]);
                }
            }
        }
        // inserting in the connected components.
        if (visited_vertices.size() >= (MDIM + 1)) {
            std::sort(visited_vertices.begin(), visited_vertices.end());
            connected_components.emplace_back();
            // inserting the vertices.
            std::unordered_map<size_t, size_t> remap_vertex;
            for (size_t i = 0; i < visited_vertices.size(); ++i) {
                remap_vertex[visited_vertices[i]] = i;
                connected_components.back().vertex_labels.push_back(approx.vertex_labels[visited_vertices[i]]);
                connected_components.back().vertices.push_back(approx.vertices[visited_vertices[i]]);
            }
            // inserting the remapped edges.
            for (size_t i = 0; i < approx.edge_connections.size(); ++i) {
                std::vector<std::array<size_t, 2>> this_edge_conn;
                if (remap_vertex.contains(approx.edge_connections[i][0]) &&
                    remap_vertex.contains(approx.edge_connections[i][1])) {
                    this_edge_conn.push_back(
                        { remap_vertex[approx.edge_connections[i][0]], remap_vertex[approx.edge_connections[i][1]] });
                }
                if (!this_edge_conn.empty()) {
                    connected_components.back().edge_hface_labels.push_back(approx.edge_hface_labels[i]);
                    for (size_t j = 0; j < this_edge_conn.size(); ++j) {
                        connected_components.back().edge_connections.push_back(this_edge_conn[j]);
                    }
                }
            }
        }
    }
    return connected_components;
}

int main(int argc, char *argv[]) {
    // specifying the input arguments.
    argparse::ArgumentParser program("skeleton_hypercubes", "", argparse::default_arguments::help);
    program.add_description(
        "This program executes the Combinatorial Skeleton for the binary output of GCCH in hypercube format.");
    program.set_usage_max_line_width(80);
    program.add_usage_newline();
    program.add_argument("-i", "--input")
        .default_value("out_gcch.bin")
        .required()
        .help("specify the input .bin file.");
    program.add_argument("-o", "--output")
        .default_value("out_gcch.pol")
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

    constexpr size_t EOF_FLAG = std::numeric_limits<size_t>::max();

    size_t g;
    std::array<UINT_COORD, NDIM> grid_coord{};
    BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH> bitset_coord{};
    HypercubeApprox approx;
    size_t n_saved_hypercubes = 0;
    size_t n_read_hypercubes = 0;

    // writing header.
    fprintf(fout, "%3lu %3lu\n", NDIM, KDIM);
    for (int64_t i = 0; i < NDIM; ++i) {
        fprintf(fout, "%3lu ", DOMAIN_DIV[i]);
    }
    fprintf(fout, "\n\n");

    while (true) {
        // reading the hypercube number.
        fread(&g, sizeof(size_t), 1, fin);

        // checking if it is eof.
        if (g == EOF_FLAG) {
            break;
        }

        ++n_read_hypercubes;

        // reading the parameters of the hypercube.
        readOutputHypercube(fin, bitset_coord, approx);
        bitsetToCoord(bitset_coord, grid_coord);

        // separating the connected components.
        const std::vector<HypercubeApprox> connected_components = separateConnectedComponents(approx);

        // if there are no connected components (because the number of vertices didn't match with the manifold dimension), skip.
        if (connected_components.empty()) {
            continue;
        }

        // saving only the region of the cut.
        // if((grid_coord[5] < (DOMAIN_DIV[5] / 2)) || (grid_coord[5] > (1 + DOMAIN_DIV[5] / 2))) {
        // // if(grid_coord[5] != (DOMAIN_DIV[5] / 2)) {
        //     continue;
        // }

        // writing the index and position of the hypercube.
        ++n_saved_hypercubes;
        fprintf(fout, "%3lu ", n_saved_hypercubes);
        for (int64_t i = 0; i < NDIM; ++i) {
            fprintf(fout, "%3lu ", static_cast<size_t>(grid_coord[i]));
        }
        fprintf(fout, "\n");

        // saving the number of connected components.
        fprintf(fout, "%3lu ", connected_components.size());
        fprintf(fout, "\n");

        // looping through all connected components.
        for (const auto &c : connected_components) {
            // saving the number of vertices.
            fprintf(fout, "%3lu ", c.vertices.size());
            fprintf(fout, "\n");

            // saving the vertices.
            for (size_t i = 0; i < c.vertices.size(); ++i) {
                for (size_t j = 0; j < KDIM + 1; ++j) {
                    fprintf(fout, "%3lu ", c.vertex_labels[i][j]);
                }
                for (size_t j = 0; j < NDIM; ++j) {
                    fprintf(fout, float_format.c_str(), c.vertices[i][j]);
                }
                fprintf(fout, "\n");
            }

            // checking the faces and connections.
            if (c.edge_hface_labels.empty()) {
                throw std::runtime_error("number of faces should not be zero");
            }

            // saving the edges.
            fprintf(fout, "%3lu ", c.edge_connections.size());
            fprintf(fout, "\n");
            for (size_t i = 0; i < c.edge_connections.size(); ++i) {
                fprintf(fout, "%3lu ", (c.edge_connections[i][0] + 1));
                fprintf(fout, "%3lu ", (c.edge_connections[i][1] + 1));
                fprintf(fout, "\n");
            }

            // creating auxiliary variables to run the combinatorial skeleton.
            std::vector<Label2N> current_hfaces{ c.edge_hface_labels };
            std::vector<std::vector<size_t>> hface_connections;
            std::vector<Label2N> next_hfaces;

            // running the combinatorial skeleton.
            for (size_t dim = KDIM + 2; dim <= NDIM; ++dim) {
                genEdges(dim,
                         current_hfaces,
                         hface_connections,
                         next_hfaces);
                fprintf(fout, "%3lu ", hface_connections.size());
                fprintf(fout, "\n");
                for (size_t i = 0; i < hface_connections.size(); ++i) {
                    for (size_t j = 0; j < hface_connections[i].size(); ++j) {
                        fprintf(fout, "%3lu ", hface_connections[i][j] + 1);
                    }
                    fprintf(fout, "\n");
                }
                current_hfaces = next_hfaces;
                hface_connections.clear();
                next_hfaces.clear();
            }
        }
        fprintf(fout, "\n");
    }

    fprintf(fout, "-1\n");

    fclose(fin);
    fclose(fout);

    std::cout << "n_read_hypercubes: " << n_read_hypercubes << std::endl;
    std::cout << "n_saved_hypercubes: " << n_saved_hypercubes << std::endl;

    return 0;
}
