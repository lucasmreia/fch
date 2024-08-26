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

// use two unordered_sets, one for those that have yet to be processed and one for those that have already been processed.
std::unordered_set<BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH>> to_be_processed;
std::unordered_set<BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH>> already_processed;

void addFirstPoint() {
    // calculate which n-hypercube the starting point belongs to.
    std::array<UINT_COORD, NDIM> grid_coord = getFirstPointHypercubeCoord();
    BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH> coord_bitset{};
    coordToBitset(grid_coord, coord_bitset);

    to_be_processed.insert(coord_bitset);
}

void continuation(const std::string &filename) {
    std::array<UINT_COORD, NDIM> grid_coord{};
    BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH> coord_bitset{};

    // create output file and save header.
    size_t n_saved_hypercubes = 0;

    FILE *fout = fopen(filename.c_str(), "wb");
    if (fout == nullptr) {
        throw std::runtime_error("ERROR OPENING FILE");
    }

    std::vector<std::array<double, NDIM>> vert_hypercube((1ULL << NDIM));
    HypercubeApprox approx;
    Label2N neighbor_cells;
    bool valid;

    // while to_be_processed not empty, take a cube and process it.
    while (!to_be_processed.empty()) {
        coord_bitset = *to_be_processed.begin();
        to_be_processed.erase(to_be_processed.begin());

        // add the current cube to those that have already been traversed.
        already_processed.insert(coord_bitset);

        // convert to coordinates.
        bitsetToCoord(coord_bitset, grid_coord);
        for (size_t j = 0; j < (1ULL << NDIM); ++j) {
            for (size_t i = 0; i < NDIM; ++i) {
                vert_hypercube[j][i] = DOMAIN_MIN[i] + DOMAIN_STEP[i] * static_cast<double>(grid_coord[i]);
                if (j & (1ULL << i)) {
                    vert_hypercube[j][i] += DOMAIN_STEP[i];
                }
            }
            addPerturbationVertex(grid_coord, j, vert_hypercube[j]);
        }

        // call gch.
        valid = gch(vert_hypercube, approx, neighbor_cells);

        // save result of the gcmh.
        if (valid) {
            ++n_saved_hypercubes;
            writeOutputHypercube(n_saved_hypercubes, coord_bitset, approx, fout);
        }

        // for each neighbor marked, add it to the list of cubes that need to be traversed.
        for (size_t i = 0; i < (2 * NDIM); ++i) {
            if (neighbor_cells[i]) {
                if ((i < NDIM) && (grid_coord[i] > 0)) {
                    std::array<UINT_COORD, NDIM> ngrid_coord = grid_coord;
                    ngrid_coord[i] -= 1;
                    coordToBitset(ngrid_coord, coord_bitset);
                    if (!already_processed.contains(coord_bitset)) {
                        to_be_processed.insert(coord_bitset);
                    }
                } else if ((i >= NDIM) && (grid_coord[i % NDIM] < (DOMAIN_DIV[i % NDIM] - 1))) {
                    std::array<UINT_COORD, NDIM> ngrid_coord = grid_coord;
                    ngrid_coord[i % NDIM] += 1;
                    coordToBitset(ngrid_coord, coord_bitset);
                    if (!already_processed.contains(coord_bitset)) {
                        to_be_processed.insert(coord_bitset);
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

    addFirstPoint();

    const auto t0 = std::chrono::high_resolution_clock::now();
    continuation(output_path);
    const auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "Time to generate bin (GCCH): " << std::chrono::duration<double>(t1 - t0).count()
              << std::endl;

    return 0;
}