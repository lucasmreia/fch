///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#include "permutahedron.h"

bool isSFaceInsideGrid(const SFace &face) {
    size_t acc_e = 0;
    for (size_t i = 0; i < KDIM; ++i) {  // I assume it's always in canonical notation, so I don't use the last one.
        acc_e |= face.e[i];
    }
    // since in the canonical notation the vertices always accumulate positively, it's enough to find out if the coordinate of the last
    // vertex is inside or on the edge of the grid.
    for (size_t j = 0; j < NDIM; ++j) {
        // I assume it's always in canonical notation, so all you have to do is add it up.
        if ((face.grid_coord[j] + static_cast<int>((acc_e & (1ULL << j)) > 0)) > DOMAIN_DIV[j]) {
            return false;  // returns indicating outside the domain grid.
        }
    }
    return true;
}

bool isSCofaceInsideGrid(const SCoface &coface) {
    // in canonical notation, the first vertex of the coface must be inside the grid.
    // as in canonical notation the vertices always accumulate positively, incrementing by 1 at each coordinate, if
    // the first one is inside the grid, it's enough to guarantee that the entire simplex is inside.
    for (size_t j = 0; j < NDIM; ++j) {
        if (coface.grid_coord[j] >= DOMAIN_DIV[j]) {
            return false;
        }
    }
    return true;
}
