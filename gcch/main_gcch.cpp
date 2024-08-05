///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#include <argparse/argparse.hpp>
#include <chrono>
#include <iostream>
#include <string>
#include <unordered_set>

#include "../common/utility.h"
#include "gch.h"

void continuation(const std::string &filename) {
    // calculate which n-hypercube the starting point belongs to.
    std::array<size_t, NDIM> grid_coord = getFirstPointHypercubeCoord();
    size_t g = coordToEnum(grid_coord, DOMAIN_GRID_BASIS);

    // use two unordered_sets, one for those that have yet to be processed and one for those that have already been processed.
    std::unordered_set<size_t> to_be_processed;
    std::unordered_set<size_t> already_processed;
    to_be_processed.insert(g);

    // create output file and save header.
    size_t n_saved_hypercubes = 0;

    FILE *fout = fopen(filename.c_str(), "wb");
    if (fout == nullptr) {
        throw std::runtime_error("ERROR OPENING FILE");
    }

    std::array<std::array<double, NDIM>, (1ULL << NDIM)> vert_hypercube{};
    HypercubeApprox approx;
    Label2N neighbor_cells;
    bool valid;

    // while to_be_processed not empty, take a cube and process it.
    while (!to_be_processed.empty()) {
        g = *to_be_processed.begin();
        to_be_processed.erase(to_be_processed.begin());

        // add the current cube to those that have already been traversed.
        already_processed.insert(g);

        // convert the cube index to coordinates.
        enumToCoord(g, DOMAIN_GRID_BASIS, grid_coord);
        for (size_t j = 0; j < (1ULL << NDIM); ++j) {
            for (size_t i = 0; i < NDIM; ++i) {
                vert_hypercube[j][i] = DOMAIN_MIN[i] + DOMAIN_RANGE[i] * static_cast<double>(grid_coord[i]);
                if (j & (1ULL << i)) {
                    vert_hypercube[j][i] += DOMAIN_RANGE[i];
                }
            }
            addPerturbationVertex(grid_coord, j, vert_hypercube[j]);
        }

        // call gch.
        valid = gch(vert_hypercube, approx, neighbor_cells);

        // save result of the gcmh.
        if (valid) {
            ++n_saved_hypercubes;
            writeOutputHypercube(g, approx, fout);
        }

        // for each neighbor marked, add it to the list of cubes that need to be traversed.
        for (size_t i = 0; i < (2 * NDIM); ++i) {
            if (neighbor_cells[i]) {
                if ((i < NDIM) && (g >= DOMAIN_GRID_BASIS[i])) {
                    const size_t ng = g - DOMAIN_GRID_BASIS[i];
                    if (!already_processed.contains(ng)) {
                        to_be_processed.insert(ng);
                    }
                } else if ((i >= NDIM) && ((g + DOMAIN_GRID_BASIS[i % NDIM]) < DOMAIN_GRID_BASIS[NDIM])) {
                    const size_t ng = g + DOMAIN_GRID_BASIS[i % NDIM];
                    if (!already_processed.contains(ng)) {
                        to_be_processed.insert(ng);
                    }
                }
            }
        }
    }

    // close the output file.
    {
        constexpr size_t EOF_FLAG = std::numeric_limits<size_t>::max();
        fwrite(&EOF_FLAG, sizeof(size_t), 1, fout);
    }
    fclose(fout);
    fout = nullptr;

    std::cout << n_saved_hypercubes << " cells generated" << std::endl;
}

int main(int argc, char *argv[]) {
    // specifying the input arguments.
    argparse::ArgumentParser program("gcch", "", argparse::default_arguments::help);
    program.add_description(
        "This program executes the GCCH.");
    program.set_usage_max_line_width(80);
    program.add_usage_newline();
    program.add_argument("-o", "--output")
        .default_value("out_gcch.bin")
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

    const auto t0 = std::chrono::high_resolution_clock::now();
    continuation(output_path);
    const auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "Time to generate bin (GCCH): " << std::chrono::duration<double>(t1 - t0).count()
              << std::endl;

    return 0;
}