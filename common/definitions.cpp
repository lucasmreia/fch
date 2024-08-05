///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#include "definitions.h"

void func(const double *x, double *f){
    return TestCaseImpl<TestCase::SELECTED_CASE>::func(x, f);
}
