#include "nonlinear.h"
#include <cmath>

constexpr double EPS = 1e-12;

namespace nonlinear {
    double newton(std::function<double(double)> f, std::function<double(double)> df, double x0, int max_iter, double tol) {
        for (int i = 0; i < max_iter; ++i) {
            double fx = f(x0);
            double dfx = df(x0);
            if (fabs(dfx) < EPS) return NAN;
            double x1 = x0 - fx / dfx;
            if (fabs(x1 - x0) < tol) return x1;
            x0 = x1;
        }
        return NAN;
    }

    double secant(std::function<double(double)> f, double x0, double x1, int max_iter, double tol) {
        for (int i = 0; i < max_iter; ++i) {
            double f0 = f(x0), f1 = f(x1);
            if (fabs(f1 - f0) < EPS) return NAN;
            double x2 = x1 - f1 * (x1 - x0) / (f1 - f0);
            if (fabs(x2 - x1) < tol) return x2;
            x0 = x1;
            x1 = x2;
        }
        return NAN;
    }

    double bisection(std::function<double(double)> f, double a, double b, int max_iter, double tol) {
        if (f(a) * f(b) >= 0) return NAN;
        for (int i = 0; i < max_iter; ++i) {
            double c = (a + b) / 2;
            double fc = f(c);
            if (fabs(fc) < tol || (b - a) / 2 < tol) return c;
            if (f(a) * fc < 0) b = c; else a = c;
        }
        return NAN;
    }

    double false_position(std::function<double(double)> f, double a, double b, int max_iter, double tol) {
        double fa = f(a), fb = f(b);
        if (fa * fb >= 0) return NAN;
        for (int i = 0; i < max_iter; ++i) {
            double c = (a * fb - b * fa) / (fb - fa);
            double fc = f(c);
            if (fabs(fc) < tol) return c;
            if (fa * fc < 0) { b = c; fb = fc; }
            else { a = c; fa = fc; }
        }
        return NAN;
    }
}
