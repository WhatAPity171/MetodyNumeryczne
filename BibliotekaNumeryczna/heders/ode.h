#pragma once
#include <vector>
#include <functional>

namespace ode {
    std::vector<double> heun(std::function<double(double, double)> f, double T0, double t_max, double dt);
    std::vector<double> midpoint(std::function<double(double, double)> f, double T0, double t_max, double dt);
    std::vector<double> rk4(std::function<double(double, double)> f, double T0, double t_max, double dt);
}