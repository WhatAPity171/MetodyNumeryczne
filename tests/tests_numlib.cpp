#include "catch_amalgamated.hpp"
#include "numlib.h"
#include <fstream>

TEST_CASE("Numlib - printVector/Matrix", "[io]") {
    std::vector<double> v = { 1.0, 2.0 };
    std::vector<std::vector<double>> m = { {1.0, 0.0}, {0.0, 1.0} };
    printVector(v, "v");
    printMatrix(m, "m");
    REQUIRE(true); // rêczne sprawdzenie wydruku
}

TEST_CASE("Numlib - loadMatrixFile (not found)", "[load]") {
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    int N;
    REQUIRE_FALSE(loadMatrixFile("nonexistent.txt", A, b, N));
}
