#pragma once
#include <vector>
#include <functional>

namespace integral {
    double horner(const std::vector<double>& coeffs, double x);
    double simpson(std::function<double(double)> f, double a, double b, int n);
    double trapezoid(std::function<double(double)> f, double a, double b, int n);
    double rectangle(std::function<double(double)> f, double a, double b, int n);
}
