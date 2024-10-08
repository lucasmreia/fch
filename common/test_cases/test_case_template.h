///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#pragma once

enum class TestCase {
    KLEIN_5TO3,
    HYPERSHPERE_5TO1,
    CIRCLE_9TO8,
    HYPERSHPERE_N_TO_K,
    S1_K_TIMES,
};

template <TestCase CASE>
struct TestCaseImpl;
