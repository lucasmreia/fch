///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#include "fch.h"

bool mirrorHFace(const std::array<size_t, NDIM> &in_grid, const Label2N &in_label,
                 const size_t bit_idx,
                 std::array<size_t, NDIM> &out_grid, Label2N &out_label) {
    const size_t dimIdx = bit_idx % NDIM;
    out_grid = in_grid;
    if (bit_idx < NDIM) {
        if (out_grid[dimIdx] <= 0) {
            return false;
        }
        out_grid[dimIdx] -= 1;
    } else {
        if (out_grid[dimIdx] >= (DOMAIN_DIV[dimIdx] - 1)) {
            return false;
        }
        out_grid[dimIdx] += 1;
    }
    out_label = in_label;
    out_label.flip(dimIdx);
    out_label.flip(dimIdx + NDIM);
    return true;
}

void writeOutputCell(const size_t hcoface_g,
                     const Label2N &hcoface_label,
                     const size_t hface_g_in,
                     const std::array<size_t, KDIM + 1> &simplex_vert_label_in,
                     const std::array<double, NDIM> &mf_vert_in,
                     const size_t hface_g_out,
                     const std::array<size_t, KDIM + 1> &simplex_vert_label_out,
                     const std::array<double, NDIM> &mf_vert_out,
                     FILE *fout) {
    fwrite(&hcoface_g, sizeof(size_t), 1, fout);
    hcoface_label.write(fout);
    fwrite(&hface_g_in, sizeof(size_t), 1, fout);
    fwrite(simplex_vert_label_in.data(), sizeof(size_t), KDIM + 1, fout);
    fwrite(mf_vert_in.data(), sizeof(double), NDIM, fout);
    fwrite(&hface_g_out, sizeof(size_t), 1, fout);
    fwrite(simplex_vert_label_out.data(), sizeof(size_t), KDIM + 1, fout);
    fwrite(mf_vert_out.data(), sizeof(double), NDIM, fout);
}

void readOutputCell(FILE *fin,
                    Label2N &hcoface_label,
                    size_t &hface_g_in,
                    std::array<size_t, KDIM + 1> &simplex_vert_label_in,
                    std::array<double, NDIM> &mf_vert_in,
                    size_t &hface_g_out,
                    std::array<size_t, KDIM + 1> &simplex_vert_label_out,
                    std::array<double, NDIM> &mf_vert_out) {
    hcoface_label.read(fin);
    fread(&hface_g_in, sizeof(size_t), 1, fin);
    fread(simplex_vert_label_in.data(), sizeof(size_t), KDIM + 1, fin);
    fread(mf_vert_in.data(), sizeof(double), NDIM, fin);
    fread(&hface_g_out, sizeof(size_t), 1, fin);
    fread(simplex_vert_label_out.data(), sizeof(size_t), KDIM + 1, fin);
    fread(mf_vert_out.data(), sizeof(double), NDIM, fin);
}
