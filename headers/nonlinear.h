#pragma once
#include <functional>

namespace nonlinear {
    double newton(std::function<double(double)> f, std::function<double(double)> df, double x0, int max_iter = 100, double tol = 1e-8);
    double secant(std::function<double(double)> f, double x0, double x1, int max_iter = 100, double tol = 1e-8);
    double bisection(std::function<double(double)> f, double a, double b, int max_iter = 100, double tol = 1e-8);
    double false_position(std::function<double(double)> f, double a, double b, int max_iter = 100, double tol = 1e-8);
}
