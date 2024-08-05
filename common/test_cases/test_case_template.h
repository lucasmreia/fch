///////////////////////////////////
// Author: Lucas M. Reia
// Original source: https://github.com/lucasmreia/fch
///////////////////////////////////

#pragma once

enum class TestCase {
    KLEIN_5TO3,
    HYPERSHPERE_5TO1,
    CIRCLE_9TO8,
};

template <TestCase CASE>
struct TestCaseImpl;
