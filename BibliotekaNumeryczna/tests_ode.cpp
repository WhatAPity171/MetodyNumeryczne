#include "catch_amalgamated.hpp"
#include "ode.h"

TEST_CASE("ODE - Heun", "[heun]") {
    auto f = [](double t, double T) { return -T; };
    auto result = ode::heun(f, 1.0, 1.0, 0.1);
    REQUIRE_FALSE(result.empty());
    REQUIRE(result.front() == Catch::Approx(1.0).epsilon(1e-12));
}

TEST_CASE("ODE - RK4", "[rk4]") {
    auto f = [](double t, double T) { return -T; };
    auto result = ode::rk4(f, 1.0, 1.0, 0.1);
    REQUIRE_FALSE(result.empty());
    REQUIRE(result.front() == Catch::Approx(1.0).epsilon(1e-12));
}
