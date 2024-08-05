///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#include <algorithm>
#include <argparse/argparse.hpp>
#include <array>
#include <bit>
#include <chrono>
#include <iostream>
#include <unordered_set>

#include "../common/bmi2.h"
#include "permutahedron.h"

std::unordered_set<BitsetPermutahedronSFace> faces_to_be_processed;
std::unordered_set<BitsetPermutahedronSFace> faces_already_processed;

size_t n_saved_simplices = 0;

void createFacesFromCoface(const SCoface &coface, FILE *fout) {
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
                    // save vertex to file.
                    ++n_saved_simplices;
                    fwrite(&n_saved_simplices, sizeof(size_t), 1, fout);
                    bitset_face.write(fout);
                    fwrite(mf_vertex.data(), sizeof(double), NDIM, fout);
                    //
                    // printing the point found if it is aligned with the grid.
                    // size_t w= 0;
                    // for(int idx = 0; idx < KDIM; ++idx) {
                    //     w += face.e[idx];
                    // }
                    // if(std::popcount(w) == KDIM) {
                    //     for(int idx = 0; idx < (NDIM - 1); ++idx) {
                    //         printf("%15.8f, ", mf_vertex[idx]);
                    //     }
                    //     printf("%15.8f", mf_vertex[NDIM - 1]);
                    //     std::cout << std::endl;
                    // }
                    //
                }
            }
        }
    }
    // combining the last vector with the first.
    // I have to change the reference vertex: advance 1 vector.
    face.grid_coord = coface.grid_coord;
    for (size_t j = 0; j < NDIM; ++j) {   // coface.grid += face.e[0], only considers canonical notation.
        if (coface.e[0] & (1ULL << j)) {  // e0 will never have the negative (e_n).
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
                // save vertex to file.
                ++n_saved_simplices;
                fwrite(&n_saved_simplices, sizeof(size_t), 1, fout);
                bitset_face.write(fout);
                fwrite(mf_vertex.data(), sizeof(double), NDIM, fout);
            }
        }
    }
}

void createCofacesFromFace(const SFace &face, FILE *fout) {
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
                    createFacesFromCoface(coface_canonical, fout);
                }
            }
            coface.e[idx] = face.e[i];  // now just assign this one to split the next ones.
            ++idx;
        }
    }
}

void addFirstPoint() {
    // calculate which n-hypercube the first point belongs to.
    constexpr std::array<size_t, NDIM> GRID = getFirstPointHypercubeCoord();

    // determine, in the permutahedral representation, which n-simplex the point belongs to.
    std::array<double, NDIM + 1> aux{};
    std::array<size_t, NDIM + 1> e_idx{};
    for (size_t i = 0; i < NDIM; ++i) {
        // taking the percentage of each coordinate of the point inside the unit hypercube.
        aux[i] = (FIRST_POINT[i] - (DOMAIN_MIN[i] + static_cast<double>(GRID[i]) * DOMAIN_RANGE[i])) / DOMAIN_RANGE[i];
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
                         // sort indices according to corresponding array element
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
    std::array<size_t, NDIM + 1> simplex_n{};
    simplex_n[0] = 0;
    for (size_t i = 1; i < (NDIM + 1); ++i) {
        simplex_n[i] = simplex_n[i - 1] | e_idx[i - 1];
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
            face_canonical.toBitset(bitset_face);
            faces_to_be_processed.insert(bitset_face);
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

void continuation_pta(const std::string &filename) {
    FILE *fout = fopen(filename.c_str(), "wb");
    if (fout == nullptr) {
        throw std::runtime_error("ERROR OPENING FILE");
    }

    BitsetPermutahedronSFace bitset_face;
    SFace face;

    addFirstPoint();

    std::cout << faces_to_be_processed.size() << std::endl;

    while (!faces_to_be_processed.empty()) {
        bitset_face = *faces_to_be_processed.begin();
        faces_to_be_processed.erase(faces_to_be_processed.begin());

        // convert the label of the k-simplex to coordinates.
        face.fromBitset(bitset_face);

        // calculate the cofaces.
        // then, for each coface, calculate the faces.
        createCofacesFromFace(face, fout);
    }

    // close output file.
    {
        constexpr size_t EOF_FLAG = std::numeric_limits<size_t>::max();
        fwrite(&EOF_FLAG, sizeof(size_t), 1, fout);
    }
    fclose(fout);
    fout = nullptr;

    std::cout << n_saved_simplices << " simplices generated" << std::endl;
}

int main(int argc, char *argv[]) {
    // specifying the input arguments.
    argparse::ArgumentParser program("pta", "", argparse::default_arguments::help);
    program.add_description(
        "This program executes the manifold tracing algorithm based on the permutahedral representation.");
    program.set_usage_max_line_width(80);
    program.add_usage_newline();
    program.add_argument("-o", "--output")
        .default_value("out_pta.bin")
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
    continuation_pta(output_path);
    const auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "Time to generate bin (PTA): " << std::chrono::duration<double>(t1 - t0).count()
              << std::endl;

    return 0;
}
