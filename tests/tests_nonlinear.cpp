#include "catch_amalgamated.hpp"
#include "nonlinear.h"
#include <cmath>

TEST_CASE("Nonlinear - Newton", "[newton]") {
    auto f = [](double x) { return x * x - 2; };
    auto df = [](double x) { return 2 * x; };
    REQUIRE(nonlinear::newton(f, df, 1.0) == Catch::Approx(std::sqrt(2)).epsilon(1e-6));
}

TEST_CASE("Nonlinear - Secant", "[secant]") {
    auto f = [](double x) { return x * x - 2; };
    REQUIRE(nonlinear::secant(f, 1.0, 2.0) == Catch::Approx(std::sqrt(2)).epsilon(1e-6));
}

TEST_CASE("Nonlinear - Bisection", "[bisection]") {
    auto f = [](double x) { return x - 3; };
    REQUIRE(nonlinear::bisection(f, 2.0, 4.0) == Catch::Approx(3.0).epsilon(1e-8));
}

TEST_CASE("Nonlinear - False Position", "[false_position]") {
    auto f = [](double x) { return x - 5; };
    REQUIRE(nonlinear::false_position(f, 4.0, 6.0) == Catch::Approx(5.0).epsilon(1e-6));
}
