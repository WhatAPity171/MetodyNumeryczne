#include "catch_amalgamated.hpp"
#include "linear.h"
#include <vector>

TEST_CASE("Linear - Gaussian Elimination", "[gauss]") {
    std::vector<std::vector<double>> A = { {2, 1}, {1, 3} };
    std::vector<double> b = { 3, 4 };
    std::vector<double> x;
    linear::gaussianElimination(A, b, x);
    REQUIRE(x[0] == Catch::Approx(1.0).epsilon(1e-6));
    REQUIRE(x[1] == Catch::Approx(1.0).epsilon(1e-6));
}

TEST_CASE("Linear - Permutations", "[permute]") {
    std::vector<int> p = { 1, 0 };
    std::vector<double> in = { 10, 20 }, out;
    linear::applyPermutation(p, in, out);
    REQUIRE(out[0] == Catch::Approx(20.0));
    REQUIRE(out[1] == Catch::Approx(10.0));
}
