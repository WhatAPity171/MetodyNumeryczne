#include "catch_amalgamated.hpp"
#include "integral.h"
#include <cmath>

TEST_CASE("Integral - Simpson", "[simpson]") {
    auto f = [](double x) { return x * x; };
    REQUIRE(integral::simpson(f, 0, 2, 10) == Catch::Approx(8.0 / 3.0).epsilon(1e-4));
}

TEST_CASE("Integral - Simpson throws on odd n", "[simpson]") {
    auto f = [](double x) { return x; };
    REQUIRE_THROWS_AS(integral::simpson(f, 0, 1, 3), std::invalid_argument);
}

TEST_CASE("Integral - Trapezoid", "[trapezoid]") {
    auto f = [](double x) { return x; };
    REQUIRE(integral::trapezoid(f, 0, 1, 10) == Catch::Approx(0.5).epsilon(1e-4));
}

TEST_CASE("Integral - Rectangle", "[rectangle]") {
    auto f = [](double x) { return 1.0; };
    REQUIRE(integral::rectangle(f, 0, 1, 100) == Catch::Approx(1.0).epsilon(1e-4));
}
