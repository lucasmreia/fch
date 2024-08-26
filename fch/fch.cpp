///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#include "fch.h"

void writeOutputCell(const size_t g,
                     const BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH> &bitset_coord,
                     const CanonicalFaceLabel &hcoface_label,
                     const BitsetSimplex<KDIM> &bitset_simplex_in,
                     const std::array<double, NDIM> &mf_vert_in,
                     const BitsetSimplex<KDIM> &bitset_simplex_out,
                     const std::array<double, NDIM> &mf_vert_out,
                     FILE *fout) {
    fwrite(&g, sizeof(size_t), 1, fout);
    bitset_coord.write(fout);
    hcoface_label.write(fout);
    bitset_simplex_in.write(fout);
    fwrite(mf_vert_in.data(), sizeof(double), NDIM, fout);
    bitset_simplex_out.write(fout);
    fwrite(mf_vert_out.data(), sizeof(double), NDIM, fout);
}

void readOutputCell(FILE *fin,
                    BitLabel<DOMAIN_DIV_TOTAL_BIT_WIDTH> &bitset_coord,
                    CanonicalFaceLabel &hcoface_label,
                    BitsetSimplex<KDIM> &bitset_simplex_in,
                    std::array<double, NDIM> &mf_vert_in,
                    BitsetSimplex<KDIM> &bitset_simplex_out,
                    std::array<double, NDIM> &mf_vert_out) {
    bitset_coord.read(fin);
    hcoface_label.read(fin);
    bitset_simplex_in.read(fin);
    fread(mf_vert_in.data(), sizeof(double), NDIM, fin);
    bitset_simplex_out.read(fin);
    fread(mf_vert_out.data(), sizeof(double), NDIM, fin);
}
