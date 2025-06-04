#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <queue>
#include <sstream>
#include <string>
#include <algorithm>
#include <functional>
#include <numeric>

namespace numlib {
//=================== Układy liniowe ===================//
namespace linear {
    inline void gaussianElimination(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& x) {
        int N = A.size();
        for (int i = 0; i < N; i++) {
            int maxRow = i;
            for (int k = i + 1; k < N; k++) {
                if (fabs(A[k][i]) > fabs(A[maxRow][i])) maxRow = k;
            }
            std::swap(A[i], A[maxRow]);
            std::swap(b[i], b[maxRow]);

            for (int k = i + 1; k < N; k++) {
                double factor = A[k][i] / A[i][i];
                for (int j = i; j < N; j++) A[k][j] -= factor * A[i][j];
                b[k] -= factor * b[i];
            }
        }

        x.assign(N, 0);
        for (int i = N - 1; i >= 0; i--) {
            x[i] = b[i];
            for (int j = i + 1; j < N; j++) x[i] -= A[i][j] * x[j];
            x[i] /= A[i][i];
        }
    }

    inline void LUDecomposition(const std::vector<std::vector<double>>& input, std::vector<std::vector<double>>& L,
                                 std::vector<std::vector<double>>& U, std::vector<int>& rowPerm, std::vector<int>& colPerm) {
        int N = input.size();
        auto A = input;
        L.assign(N, std::vector<double>(N, 0));
        U.assign(N, std::vector<double>(N, 0));
        rowPerm.resize(N);
        colPerm.resize(N);
        std::iota(rowPerm.begin(), rowPerm.end(), 0);
        std::iota(colPerm.begin(), colPerm.end(), 0);

        for (int k = 0; k < N; k++) {
            double maxVal = 0;
            int maxRow = k, maxCol = k;
            for (int i = k; i < N; i++) {
                for (int j = k; j < N; j++) {
                    if (fabs(A[i][j]) > fabs(maxVal)) {
                        maxVal = A[i][j];
                        maxRow = i;
                        maxCol = j;
                    }
                }
            }
            if (fabs(maxVal) < EPS) throw std::runtime_error("Macierz osobliwa");

            std::swap(A[k], A[maxRow]);
            std::swap(rowPerm[k], rowPerm[maxRow]);
            for (int i = 0; i < N; i++) std::swap(A[i][k], A[i][maxCol]);
            std::swap(colPerm[k], colPerm[maxCol]);

            for (int i = k + 1; i < N; i++) {
                A[i][k] /= A[k][k];
                for (int j = k + 1; j < N; j++) A[i][j] -= A[i][k] * A[k][j];
            }
        }

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i > j) L[i][j] = A[i][j];
                else if (i == j) {
                    L[i][j] = 1.0;
                    U[i][j] = A[i][j];
                } else U[i][j] = A[i][j];
            }
        }
    }

    inline void forwardSubstitution(const std::vector<std::vector<double>>& L, const std::vector<double>& b, std::vector<double>& z) {
        int N = L.size();
        z.resize(N);
        for (int i = 0; i < N; i++) {
            z[i] = b[i];
            for (int j = 0; j < i; j++) z[i] -= L[i][j] * z[j];
        }
    }

    inline void backwardSubstitution(const std::vector<std::vector<double>>& U, const std::vector<double>& z, std::vector<double>& x) {
        int N = U.size();
        x.resize(N);
        for (int i = N - 1; i >= 0; i--) {
            x[i] = z[i];
            for (int j = i + 1; j < N; j++) x[i] -= U[i][j] * x[j];
            x[i] /= U[i][i];
        }
    }

    inline void applyPermutation(const std::vector<int>& perm, const std::vector<double>& in, std::vector<double>& out) {
        int N = perm.size();
        out.resize(N);
        for (int i = 0; i < N; i++) out[i] = in[perm[i]];
    }

    inline void undoPermutation(const std::vector<int>& perm, const std::vector<double>& in, std::vector<double>& out) {
        int N = perm.size();
        out.resize(N);
        for (int i = 0; i < N; i++) out[perm[i]] = in[i];
    }
}

//=================== Całkowanie ===================//
namespace integral {
    inline double horner(const std::vector<double>& coeffs, double x) {
        double result = coeffs[0];
        for (size_t i = 1; i < coeffs.size(); i++) result = result * x + coeffs[i];
        return result;
    }

    inline double simpson(std::function<double(double)> f, double a, double b, int n) {
        if (n % 2 != 0) throw std::invalid_argument("n musi być parzyste dla Simpsona");
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
//=================== Równania nieliniowe ===================//
namespace nonlinear {
    inline double newton(std::function<double(double)> f, std::function<double(double)> df, double x0, int max_iter = 100, double tol = 1e-8) {
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

    inline double secant(std::function<double(double)> f, double x0, double x1, int max_iter = 100, double tol = 1e-8) {
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

    inline double bisection(std::function<double(double)> f, double a, double b, int max_iter = 100, double tol = 1e-8) {
        if (f(a) * f(b) >= 0) return NAN;
        for (int i = 0; i < max_iter; ++i) {
            double c = (a + b) / 2;
            double fc = f(c);
            if (fabs(fc) < tol || (b - a) / 2 < tol) return c;
            if (f(a) * fc < 0) b = c; else a = c;
        }
        return NAN;
    }

    inline double false_position(std::function<double(double)> f, double a, double b, int max_iter = 100, double tol = 1e-8) {
        double fa = f(a), fb = f(b);
        if (fa * fb >= 0) return NAN;
        for (int i = 0; i < max_iter; ++i) {
            double c = (a * fb - b * fa) / (fb - fa);
            double fc = f(c);
            if (fabs(fc) < tol) return c;
            if (fa * fc < 0) { b = c; fb = fc; } else { a = c; fa = fc; }
        }
        return NAN;
    }
}
namespace ode {
    // f: funkcja pochodnej, T0: wartość początkowa, t_max: koniec przedziału, dt: krok czasu
    inline std::vector<double> heun(std::function<double(double, double)> f, double T0, double t_max, double dt) {
        std::vector<double> temps;
        double T = T0;
        for (double t = 0.0; t <= t_max + 1e-9; t += dt) {
            temps.push_back(T);
            double k1 = f(t, T);
            double k2 = f(t + dt, T + dt * k1);
            T += 0.5 * dt * (k1 + k2);
        }
        return temps;
    }

    inline std::vector<double> midpoint(std::function<double(double, double)> f, double T0, double t_max, double dt) {
        std::vector<double> temps;
        double T = T0;
        for (double t = 0.0; t <= t_max + 1e-9; t += dt) {
            temps.push_back(T);
            double k1 = f(t, T);
            double k2 = f(t + dt / 2.0, T + dt * k1 / 2.0);
            T += dt * k2;
        }
        return temps;
    }

    inline std::vector<double> rk4(std::function<double(double, double)> f, double T0, double t_max, double dt) {
        std::vector<double> temps;
        double T = T0;
        for (double t = 0.0; t <= t_max + 1e-9; t += dt) {
            temps.push_back(T);
            double k1 = f(t, T);
            double k2 = f(t + dt / 2.0, T + dt * k1 / 2.0);
            double k3 = f(t + dt / 2.0, T + dt * k2 / 2.0);
            double k4 = f(t + dt, T + dt * k3);
            T += dt / 6.0 * (k1 + 2 * k2 + 2 * k3 + k4);
        }
        return temps;
    }
}

} // namespace numlib


