///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#include <argparse/argparse.hpp>
#include <cstdio>
#include <map>
#include <string>
#include <unordered_map>

#include "../common/definitions.h"
#include "gch.h"

std::map<BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH>, HypercubeApprox> all_approx;

void readApproxOutput(const std::string &filename_in) {
    FILE *fin = fopen(filename_in.c_str(), "rb");

    if (fin == nullptr) {
        std::cout << "input file not found" << std::endl;
        return;
    }

    constexpr size_t EOF_FLAG = std::numeric_limits<size_t>::max();

    size_t g;
    BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH> bitset_coord{};
    HypercubeApprox approx;
    size_t n_saved_hypercubes = 0;

    while (true) {
        // reading the hypercube number.
        fread(&g, sizeof(size_t), 1, fin);

        // checking if it is eof.
        if (g == EOF_FLAG) {
            break;
        }

        ++n_saved_hypercubes;

        // reading the parameters of the hypercube.
        readOutputHypercube(fin, bitset_coord, approx);

        // inserting in the map.
        all_approx.insert({ bitset_coord, approx });
    }

    fclose(fin);
}

void saveApproxOutput(const std::string &filename_out) {
    FILE *fout = fopen(filename_out.c_str(), "wb");

    // writing all approx.
    size_t count = 1;
    for (const auto &approx : all_approx) {
        writeOutputHypercube(count++, approx.first, approx.second, fout);
    }

    // writing eof.
    {
        constexpr size_t EOF_FLAG = std::numeric_limits<size_t>::max();
        fwrite(&EOF_FLAG, sizeof(size_t), 1, fout);
    }

    fclose(fout);
}

void reorderApprox() {
    for (auto &approx : all_approx) {
        // sorting the vertices by label.
        std::map<std::array<size_t, KDIM + 1>, size_t> vert_label_map;
        for (size_t i = 0; i < approx.second.vertex_labels.size(); ++i) {
            vert_label_map.insert({ approx.second.vertex_labels[i], i });
        }
        // taking the index correspondence before sorting -> after sorting.
        std::vector<size_t> reordered_vert_idx;
        std::unordered_map<size_t, size_t> reordered_vert_idx_map;
        for (const auto &v : vert_label_map) {
            reordered_vert_idx_map.insert({ v.second, reordered_vert_idx.size() });
            reordered_vert_idx.push_back(v.second);
        }
        // creating the ordered vertices.
        std::vector<std::array<double, NDIM> > vertices;
        std::vector<std::array<size_t, KDIM + 1> > vertex_labels;
        for (size_t i = 0; i < reordered_vert_idx.size(); ++i) {
            vertices.push_back(approx.second.vertices[reordered_vert_idx[i]]);
            vertex_labels.push_back(approx.second.vertex_labels[reordered_vert_idx[i]]);
        }
        // replacing in approx.
        approx.second.vertices = vertices;
        approx.second.vertex_labels = vertex_labels;
        // updating the vertex indices in the edges.
        for (size_t i = 0; i < approx.second.edge_connections.size(); ++i) {
            approx.second.edge_connections[i][0] = reordered_vert_idx_map[approx.second.edge_connections[i][0]];
            approx.second.edge_connections[i][1] = reordered_vert_idx_map[approx.second.edge_connections[i][1]];
            // reordering the vertices of this edge.
            std::sort(approx.second.edge_connections[i].begin(), approx.second.edge_connections[i].end());
        }
        // sorting the edges by the indices of their vertices.
        std::map<std::array<size_t, 2>, size_t> edge_vert_map;
        for (size_t i = 0; i < approx.second.edge_connections.size(); ++i) {
            edge_vert_map.insert({ approx.second.edge_connections[i], i });
        }
        // taking the index correspondence before sorting -> after sorting.
        std::vector<size_t> reordered_edge_idx;
        std::unordered_map<size_t, size_t> reordered_edge_idx_map;
        for (const auto &v : edge_vert_map) {
            reordered_edge_idx_map.insert({ v.second, reordered_edge_idx.size() });
            reordered_edge_idx.push_back(v.second);
        }
        // creating the ordered edges.
        std::vector<Label2N> edge_hface_labels;
        std::vector<std::array<size_t, 2> > edge_connections;
        for (size_t i = 0; i < reordered_edge_idx.size(); ++i) {
            edge_hface_labels.push_back(approx.second.edge_hface_labels[reordered_edge_idx[i]]);
            edge_connections.push_back(approx.second.edge_connections[reordered_edge_idx[i]]);
        }
        // replacing in approx.
        approx.second.edge_hface_labels = edge_hface_labels;
        approx.second.edge_connections = edge_connections;
    }
}

int main(int argc, char *argv[]) {
    // specifying the input arguments.
    argparse::ArgumentParser program("reorder_hypercubes", "", argparse::default_arguments::help);
    program.add_description(
        "This program sorts, in ascending order, the binary output of GCCH in hypercube format. The execution of this program is not required; it is used only when you want to compare between the outputs of GCCH and FCH.");
    program.set_usage_max_line_width(80);
    program.add_usage_newline();
    program.add_argument("-i", "--input")
        .default_value("out_gcch.bin")
        .required()
        .help("specify the input .bin file.");
    program.add_argument("-o", "--output")
        .default_value("out_gcch_reordered.bin")
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

    readApproxOutput(input_path);
    reorderApprox();
    saveApproxOutput(output_path);

    return 0;
}
