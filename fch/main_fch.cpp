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

std::unordered_map<BitsetSimplex<KDIM>, Label2N> cosimplices_already_processed;
DequeUniqueElements simplices_to_be_processed;

size_t n_saved_cells = 0;

void computeCosimplexAndSimplex(const std::array<size_t, KDIM + 2> &cosimplex_vert_labels,
                                const size_t extra_vertex_idx,
                                Simplex<KDIM> &simplex,
                                Simplex<KDIM + 1> &cosimplex) {
    // --- error check.
    if (!std::is_sorted(cosimplex_vert_labels.begin(), cosimplex_vert_labels.end())) {
        std::cout << "coface_vertices must be sorted" << std::endl;
        throw std::runtime_error("coface_vertices not sorted");
    }
    // ---

    std::array<size_t, KDIM + 1> simplex_vert_labels{};
    for (size_t i = 0, v_idx = 0; i < (KDIM + 1); ++i, ++v_idx) {
        if (v_idx == extra_vertex_idx) {
            ++v_idx;
        }
        simplex_vert_labels[i] = cosimplex_vert_labels[v_idx];
    }

    computeSimplexLabelAndInductive(cosimplex_vert_labels,
                                    cosimplex.hface.label,
                                    cosimplex.e_idx);
    computeSimplexLabelAndInductive(simplex_vert_labels,
                                    simplex.hface.label,
                                    simplex.e_idx);
}

size_t getCosimplexBitIdx(const Simplex<KDIM> &simplex_canonical,
                          const Simplex<KDIM + 1> &cosimplex_canonical) {
    Label2N aux = simplex_canonical.hface.label;
    aux ^= cosimplex_canonical.hface.label;
    if (aux.count() != 1) {
        throw std::runtime_error("getCosimplexBitIdx: not valid cosimplex of the given simplex");
    }
    size_t msb, lsb;
    aux.getMsbAndLsb(msb, lsb);
    const size_t idx_flip = msb ? (std::countr_zero(msb) + NDIM) : std::countr_zero(lsb);
    if ((idx_flip >= NDIM) || (simplex_canonical.hface.grid_coord[idx_flip] == cosimplex_canonical.hface.grid_coord[idx_flip])) {
        return idx_flip;
    } else {
        return ((idx_flip + NDIM) % (2 * NDIM));
    }
}

void createCosimplicesFromSimplex(const BitsetSimplex<KDIM> &bitset_simplex_in,
                                  const Simplex<KDIM> &simplex_in_canonical,
                                  const std::array<size_t, KDIM + 1> &simplex_vert_labels,
                                  const std::array<double, NDIM> &mf_vert_in,
                                  FILE *fout) {
    Simplex<KDIM + 1> cosimplex{}, cosimplex_canonical{}, cosimplex_out{}, cosimplex_out_canonical{};
    BitsetSimplex<KDIM + 1> bitset_cosimplex, bitset_cosimplex_out;
    std::array<size_t, KDIM + 2> cosimplex_vert_labels{};
    std::array<std::array<double, NDIM>, KDIM + 2> cosimplex_vertices{};

    Simplex<KDIM> simplex_out{}, simplex_out_canonical{};
    std::array<size_t, KDIM + 1> simplex_vert_labels_out{};
    std::array<double, NDIM> mf_vert_out{};
    BitsetSimplex<KDIM> bitset_simplex_out;

    Simplex<KDIM> mirrored_simplex{};
    mirrored_simplex.e_idx = simplex_in_canonical.e_idx;

    const P_MASK mask = pmask(simplex_in_canonical.hface.label.to_ullong());

    std::pair<BitsetSimplex<KDIM>, Label2N> face_in;
    {
        auto iter_in = cosimplices_already_processed.find(bitset_simplex_in);
        face_in = *iter_in;
        cosimplices_already_processed.erase(iter_in);
    }

    for (size_t i = 0; i < (NDIM - KDIM); ++i) {
        const size_t unset_bit = pdep(1ULL << i, mask);
        const size_t bit_idx = std::countr_zero(unset_bit);
        if (face_in.second[bit_idx]) {
            cosimplex.hface.grid_coord = simplex_in_canonical.hface.grid_coord;
            // define face, and coface vectors.
            size_t extra_vertex_idx = elevateDimension(simplex_in_canonical.hface.label, simplex_in_canonical.e_idx,
                                                       bit_idx,
                                                       cosimplex.hface.label, cosimplex.e_idx);
            // coface to canonical.
            cosimplex.toCanonical(cosimplex_canonical);
            // check if coface is not in cofaces_already_processed.
            cosimplex_canonical.toBitset(bitset_cosimplex);

            // add coface to cofaces_already_processed.
            face_in.second.reset(bit_idx);
            // get the vertices labels and positions.
            inductiveToVertLabels(cosimplex_canonical.hface.label.getMsb(), cosimplex_canonical.e_idx,
                                  cosimplex_vert_labels);
            cosimplex_canonical.toVertices(cosimplex_vertices);
            // call TraverseCoface and assign mfVert_out.
            traverseHCoface<KDIM>(cosimplex_canonical.hface.grid_coord,
                                  cosimplex_vertices, cosimplex_vert_labels,
                                  extra_vertex_idx,
                                  mf_vert_out);
            // compute the Simplex<KDIM> of the k-face where the manifold leaves the (k+1)-face.
            simplex_out.hface.grid_coord = cosimplex_canonical.hface.grid_coord;
            cosimplex_out.hface.grid_coord = cosimplex_canonical.hface.grid_coord;
            computeCosimplexAndSimplex(cosimplex_vert_labels, extra_vertex_idx,
                                       simplex_out, cosimplex_out);
            // convert the output simplex to canonical notation.
            simplex_out.toCanonical(simplex_out_canonical);
            cosimplex_out.toCanonical(cosimplex_out_canonical);
            const size_t bit_idx_out = getCosimplexBitIdx(simplex_out_canonical, cosimplex_out_canonical);
            // calculate the vert labels of the output simplex.
            inductiveToVertLabels(simplex_out_canonical.hface.label.getMsb(), simplex_out_canonical.e_idx,
                                  simplex_vert_labels_out);
            // add the output simplex to faces_to_be_processed.
            simplex_out_canonical.toBitset(bitset_simplex_out);
            simplices_to_be_processed.push_back(bitset_simplex_out);
            // add the output coface to cofaces_already_processed.
            cosimplex_out_canonical.toBitset(bitset_cosimplex_out);
            {
                auto iter_out = cosimplices_already_processed.find(bitset_simplex_out);
                if (iter_out == cosimplices_already_processed.end()) {
                    size_t fixed_coord, free_coord;
                    simplex_out_canonical.hface.label.getFixedAndFreeCoord(fixed_coord, free_coord);
                    Label2N all_neighbors;
                    all_neighbors |= fixed_coord;
                    all_neighbors <<= NDIM;
                    all_neighbors |= fixed_coord;
                    iter_out = cosimplices_already_processed.emplace(bitset_simplex_out, all_neighbors).first;
                } else if (!iter_out->second[bit_idx_out]) {
                    throw std::runtime_error("output simplex cannot be already processed");
                }
                iter_out->second.reset(bit_idx_out);
            }
            // save edge: mfVert_in -> mfVert_out, coface (grid, label, vertices, the same way marching hypercubes).
            ++n_saved_cells;
            writeOutputCell(coordToEnum(cosimplex_canonical.hface.grid_coord, DOMAIN_GRID_BASIS),
                            cosimplex_out_canonical.hface.label,
                            coordToEnum(simplex_in_canonical.hface.grid_coord, DOMAIN_GRID_BASIS),
                            simplex_vert_labels,
                            mf_vert_in,
                            coordToEnum(simplex_out_canonical.hface.grid_coord, DOMAIN_GRID_BASIS),
                            simplex_vert_labels_out,
                            mf_vert_out,
                            fout);
        }

        const size_t mirrored_bit_idx = (bit_idx + NDIM) % (2 * NDIM);
        if (!face_in.second[mirrored_bit_idx]) {
            continue;
        } else if (mirrorHFace(simplex_in_canonical.hface.grid_coord, simplex_in_canonical.hface.label,
                               bit_idx,
                               mirrored_simplex.hface.grid_coord, mirrored_simplex.hface.label)) {
            cosimplex.hface.grid_coord = mirrored_simplex.hface.grid_coord;
            // define face, and coface vectors.
            size_t mirrored_extra_vertex_idx = elevateDimension<KDIM>(mirrored_simplex.hface.label,
                                                                      mirrored_simplex.e_idx,
                                                                      (bit_idx + NDIM) % (2 * NDIM),
                                                                      cosimplex.hface.label, cosimplex.e_idx);
            // coface to canonical.
            cosimplex.toCanonical(cosimplex_canonical);
            // check if coface is not in cofaces_already_processed.
            cosimplex_canonical.toBitset(bitset_cosimplex);

            // add coface to cofaces_already_processed.
            face_in.second.reset(mirrored_bit_idx);
            // get the vertices labels and positions.
            inductiveToVertLabels(cosimplex_canonical.hface.label.getMsb(), cosimplex_canonical.e_idx,
                                  cosimplex_vert_labels);
            cosimplex_canonical.toVertices(cosimplex_vertices);
            // call TraverseCoface and assign mfVert_out.
            traverseHCoface<KDIM>(cosimplex_canonical.hface.grid_coord,
                                  cosimplex_vertices, cosimplex_vert_labels,
                                  mirrored_extra_vertex_idx,
                                  mf_vert_out);
            // compute the Simplex<KDIM> of the k-face where the manifold leaves the (k+1)-face.
            simplex_out.hface.grid_coord = cosimplex_canonical.hface.grid_coord;
            cosimplex_out.hface.grid_coord = cosimplex_canonical.hface.grid_coord;
            computeCosimplexAndSimplex(cosimplex_vert_labels, mirrored_extra_vertex_idx,
                                       simplex_out, cosimplex_out);
            // convert the output simplex to canonical notation.
            simplex_out.toCanonical(simplex_out_canonical);
            cosimplex_out.toCanonical(cosimplex_out_canonical);
            const size_t bit_idx_out = getCosimplexBitIdx(simplex_out_canonical, cosimplex_out_canonical);
            // calculate the vert labels of the output simplex.
            inductiveToVertLabels(simplex_out_canonical.hface.label.getMsb(), simplex_out_canonical.e_idx,
                                  simplex_vert_labels_out);
            // add the output simplex to faces_to_be_processed.
            simplex_out_canonical.toBitset(bitset_simplex_out);
            simplices_to_be_processed.push_back(bitset_simplex_out);
            // add the output coface to cofaces_already_processed.
            cosimplex_out_canonical.toBitset(bitset_cosimplex_out);
            {
                auto iter_out = cosimplices_already_processed.find(bitset_simplex_out);
                if (iter_out == cosimplices_already_processed.end()) {
                    size_t fixed_coord, free_coord;
                    simplex_out_canonical.hface.label.getFixedAndFreeCoord(fixed_coord, free_coord);
                    Label2N all_neighbors;
                    all_neighbors |= fixed_coord;
                    all_neighbors <<= NDIM;
                    all_neighbors |= fixed_coord;
                    iter_out = cosimplices_already_processed.emplace(bitset_simplex_out, all_neighbors).first;
                } else if (!iter_out->second[bit_idx_out]) {
                    throw std::runtime_error("output simplex cannot be already processed");
                }
                iter_out->second.reset(bit_idx_out);
            }
            // save edge: mfVert_in -> mfVert_out, coface (grid, label, vertices, the same way marching hypercubes).
            ++n_saved_cells;
            writeOutputCell(coordToEnum(cosimplex_canonical.hface.grid_coord, DOMAIN_GRID_BASIS),
                            cosimplex_out_canonical.hface.label,
                            coordToEnum(simplex_in_canonical.hface.grid_coord, DOMAIN_GRID_BASIS),
                            simplex_vert_labels,
                            mf_vert_in,
                            coordToEnum(simplex_out_canonical.hface.grid_coord, DOMAIN_GRID_BASIS),
                            simplex_vert_labels_out,
                            mf_vert_out,
                            fout);
        } else {
            face_in.second.reset(mirrored_bit_idx);
        }
    }
    if (face_in.second.any()) {
        throw std::runtime_error("after processing all cofaces, face_in.second must be zero");
    }
}

void addFirstPoint() {
    // calculate which n-hypercube the first point belongs to.
    constexpr std::array<size_t, NDIM> GRID_COORD = getFirstPointHypercubeCoord();

    // determine, in permutahedron notation, which n-simplex the point belongs to.
    std::array<double, NDIM> grid_pct{};  // percentage of the position of the point inside the unit hypercube.
    std::array<size_t, NDIM> e_idx{};     // vector indices.
    std::array<size_t, NDIM> e{};         // binary representation of the vectors of the whole n-simplex.
    std::array<size_t, KDIM> e_face{};    // binary representation of the vectors of the face's simplex.
    for (size_t i = 0; i < NDIM; ++i) {
        // taking the percentage of each coordinate of the point inside the unit hypercube.
        grid_pct[i] = (FIRST_POINT[i] - (DOMAIN_MIN[i] + static_cast<double>(GRID_COORD[i]) * DOMAIN_RANGE[i])) /
                      DOMAIN_RANGE[i];
        e_idx[i] = i;
    }
    // printing aux to help with debugging if no point is found.
    std::cout << "FIRST_POINT coords inside the unitary hypercube: (";
    for (size_t i = 0; i < NDIM; ++i) {
        std::cout << grid_pct[i];
        if (i < (NDIM - 1)) {
            std::cout << ", ";
        }
    }
    std::cout << ")" << std::endl;
    std::cout << "if two coords are equal, might not find the manifold (may get the wrong n-simplex)" << std::endl;
    // argsort: rearranges e_idx by sorting aux in decreasing order.
    std::stable_sort(e_idx.begin(), e_idx.end(),
                     [&grid_pct](int left, int right) -> bool {
                         // sort indices according to corresponding array element
                         return grid_pct[left] > grid_pct[right];
                     });
    // now e_idx represents the indices of the ordered vectors of the n-simplex that contains the first point.
    // converting it to the binary representation.
    for (size_t i = 0; i < NDIM; ++i) {
        e[i] = (1ULL << e_idx[i]);
    }

    // now decompose the n-simplex in k-simplices and mark the faces to be processed.
    // n-simplex has n+1 vectors. k-simplex has k+1 vectors. I need to add (n+1)-(k+1)=n-k vectors to
    // reduce the dimension.
    // if I represent them as indices of the vertices of the cube, I can do nchoosek(n+1, k+1), and then go back to the
    // permutahedron notation.
    // calculating the vertices of the n-simplex.
    std::array<size_t, NDIM + 1> n_simplex_vert_labels{};
    n_simplex_vert_labels[0] = 0;
    for (size_t i = 1; i < (NDIM + 1); ++i) {
        n_simplex_vert_labels[i] = n_simplex_vert_labels[i - 1] | e[i - 1];
    }

    // now I need to do nchoosek(n+1, k+1).
    // making n choose k (actually (n+1) choose (k+1)) and taking the face.
    Simplex<KDIM> simplex, simplex_canonical;
    BitsetSimplex<KDIM> bitset_simplex;
    std::array<size_t, KDIM + 1> idx{};
    for (size_t i = 0; i < (KDIM + 1); ++i) {
        idx[i] = i;
    }
    while (idx[0] <= ((NDIM + 1) - (KDIM + 1))) {
        // assigning the grid and part of the face of simplex k.
        simplex.hface.grid_coord = GRID_COORD;
        simplex.hface.label.reset();
        for (size_t j = 0; j < NDIM; ++j) {
            if (n_simplex_vert_labels[idx[0]] & (1ULL << j)) {
                simplex.hface.label.set(j + NDIM);
            }
        }

        // assigning the vectors of the k-simplex.
        for (size_t j = 0; j < KDIM; ++j) {
            e_face[j] = n_simplex_vert_labels[idx[j + 1]] - n_simplex_vert_labels[idx[j]];
        }

        const size_t free_coord = n_simplex_vert_labels[idx[KDIM]] - n_simplex_vert_labels[idx[0]];
        const size_t fixed_coord = ((1ULL << NDIM) - 1) ^ free_coord;

        // popcount(free_coord) must be equals k.
        if (std::popcount(free_coord) == KDIM) {
            // assigning the vectors of the face.
            for (size_t j = 0; j < KDIM; ++j) {
                if (std::popcount(e_face[j]) != 1) {
                    throw std::runtime_error("wrong number of vectors");
                }
                simplex.e_idx[j] = std::countr_zero(e_face[j]);
            }

            // assigning the rest of the face label.
            for (size_t j = 0; j < NDIM; ++j) {
                if ((fixed_coord & (1ULL << j)) && !simplex.hface.label[j + NDIM]) {
                    simplex.hface.label.set(j);
                }
            }

            // converting to canonical and adding to faces_to_be_processed.
            simplex.toCanonical(simplex_canonical);
            simplex_canonical.toBitset(bitset_simplex);
            simplices_to_be_processed.push_back(bitset_simplex);

            // initializing the already processed.
            Label2N all_neighbors;
            all_neighbors |= fixed_coord;
            all_neighbors <<= NDIM;
            all_neighbors |= fixed_coord;
            cosimplices_already_processed.emplace(bitset_simplex, all_neighbors);
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

void continuation_fch(const std::string &filename) {
    FILE *fout = fopen(filename.c_str(), "wb");
    if (fout == nullptr) {
        throw std::runtime_error("ERROR OPENING FILE");
    }

    BitsetSimplex<KDIM> bitset_simplex;
    Simplex<KDIM> simplex;
    std::array<std::array<double, NDIM>, KDIM + 1> vertices{};
    std::array<size_t, KDIM + 1> vert_labels{};
    std::array<double, NDIM> mf_vert{};

    addFirstPoint();

    std::cout << simplices_to_be_processed.size() << std::endl;

    while (!simplices_to_be_processed.empty()) {
        bitset_simplex = simplices_to_be_processed.pop_front();

        // convert the label of the k-simplex to coordinates.
        simplex.fromBitset(bitset_simplex);

        // calculate simplex vertices.
        simplex.toVertices(vertices);
        inductiveToVertLabels(simplex.hface.label.getMsb(), simplex.e_idx, vert_labels);

        // calculate mfVert.
        if (solve(vertices, mf_vert)) {
            // calculate the cofaces.
            // then, for each coface, calculate the faces.
            createCosimplicesFromSimplex(bitset_simplex, simplex, vert_labels, mf_vert, fout);
        } else {
            cosimplices_already_processed.erase(bitset_simplex);
        }
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

    const auto t0 = std::chrono::high_resolution_clock::now();
    continuation_fch(output_path);
    const auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << "Time to generate bin (FCH): " << std::chrono::duration<double>(t1 - t0).count()
              << std::endl;

    return 0;
}