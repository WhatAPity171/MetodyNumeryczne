#include "catch_amalgamated.hpp"
#include "interp.h"

TEST_CASE("Interp - Lagrange", "[lagrange]") {
    double X[4] = { 0, 1, 2, 3 };
    double Y[4] = { 0, 1, 4, 9 };
    REQUIRE(interp::lagrange(1, 4, X, Y) == Catch::Approx(1.0).epsilon(1e-4));
}

TEST_CASE("Interp - Lagrange (edge case)", "[lagrange]") {
    double X[4] = { 0, 1, 2, 3 };
    double Y[4] = { 1, 1, 1, 1 };
    REQUIRE(interp::lagrange(2, 4, X, Y) == Catch::Approx(1.0).epsilon(1e-4));
}
