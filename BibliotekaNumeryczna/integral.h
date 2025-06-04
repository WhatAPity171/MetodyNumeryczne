#pragma once
#include "numlib.h"
namespace integral {
    inline double horner(const std::vector<double>& coeffs, double x) {
        double result = coeffs[0];
        for (size_t i = 1; i < coeffs.size(); i++) result = result * x + coeffs[i];
        return result;
    }

    inline double simpson(std::function<double(double)> f, double a, double b, int n) {
        if (n % 2 != 0) throw std::invalid_argument("n musi byæ parzyste dla Simpsona");
        double h = (b - a) / n;
        double sum = f(a) + f(b);
        for (int i = 1; i < n; ++i) {
            double x = a + i * h;
            sum += (i % 2 == 0 ? 2 : 4) * f(x);
        }
        return sum * h / 3.0;
    }

    inline double trapezoid(std::function<double(double)> f, double a, double b, int n) {
        double h = (b - a) / n;
        double sum = 0.5 * (f(a) + f(b));
        for (int i = 1; i < n; ++i) sum += f(a + i * h);
        return sum * h;
    }

    inline double rectangle(std::function<double(double)> f, double a, double b, int n) {
        double h = (b - a) / n;
        double sum = 0.0;
        for (int i = 0; i < n; ++i) sum += f(a + i * h);
        return sum * h;
    }
}