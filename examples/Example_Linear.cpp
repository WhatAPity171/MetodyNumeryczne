#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include "libnum.h"

constexpr double pi = 3.14159265358979323846;

int main() {
    // === linear: rozwiązanie Ax = b ===
    std::vector<std::vector<double>> A = {
        {2, -1, 1},
        {3, 3, 9},
        {3, 3, 5}
    };
    std::vector<double> b = { 2, -1, 4 };

    // 1. Eliminacja Gaussa
    std::vector<double> x_gauss;
    auto A1 = A;
    auto b1 = b;
    linear::gaussianElimination(A1, b1, x_gauss);
    std::cout << "Rozwiazanie (Gauss): ";
    for (double xi : x_gauss) std::cout << xi << " ";
    std::cout << "\n";

    // 2. Rozkład LU z permutacjami
    std::vector<std::vector<double>> L, U;
    std::vector<int> rowPerm, colPerm;
    linear::LUDecomposition(A, L, U, rowPerm, colPerm);

    std::vector<double> b_perm;
    linear::applyPermutation(rowPerm, b, b_perm);

    std::vector<double> z, y_perm, y;
    linear::forwardSubstitution(L, b_perm, z);
    linear::backwardSubstitution(U, z, y_perm);
    linear::undoPermutation(colPerm, y_perm, y);

    std::cout << "Rozwiazanie (LU):    ";
    for (double yi : y) std::cout << yi << " ";
    std::cout << "\n";

    // === integral: obliczenia całek ===
    auto f = [](double x) { return std::sin(x); };
    double a = 0.0, b_int = pi;
    int n = 100;

    double rect = integral::rectangle(f, a, b_int, n);
    double trap = integral::trapezoid(f, a, b_int, n);
    double simp = integral::simpson(f, a, b_int, n);

    std::cout << "Rectangle integral: " << rect << "\n";
    std::cout << "Trapezoid integral: " << trap << "\n";
    std::cout << "Simpson integral:   " << simp << "\n";

    // Horner: f(x) = 2x^3 - 3x^2 + x - 5
    std::vector<double> coeffs = { 2, -3, 1, -5 };
    double val = integral::horner(coeffs, 2.0);
    std::cout << "Horner (x=2):        " << val << "\n";

    const int N = 9;
    double arrX[N] = { 0, 1, 2, 3, 4, 5, 6, 7, 8 };
    double arrY[N] = { 1, 2, 4, 8, 16, 32, 64, 128, 256 };  // f(x) = 2^x

    int pointIndex = 4; // interpolujemy dla x = arrX[4] = 4
    double interpolated = interp::lagrange(pointIndex, N, arrX, arrY);

    std::cout << "Interpolacja Lagrange'a dla x = " << arrX[pointIndex]
        << " daje wynik: " << interpolated << "\n";

    // === ODE: rozwiązywanie równań różniczkowych zwyczajnych ===
    std::cout << "\n=== Rozwiazywanie rownan rozniczkowych ===\n";

    // Przykład 1: Wzrost populacji dT/dt = k*T, gdzie k=0.1
    auto population_growth = [](double t, double T) { return 0.1 * T; };
    double T0 = 100.0;
    double t_max = 10.0;
    double dt = 0.5;

    std::vector<double> heun_result = ode::heun(population_growth, T0, t_max, dt);
    std::vector<double> midpoint_result = ode::midpoint(population_growth, T0, t_max, dt);
    std::vector<double> rk4_result = ode::rk4(population_growth, T0, t_max, dt);

    std::cout << "Wzrost populacji dT/dt = 0.1*T, T(0) = " << T0 << "\n";
    std::cout << "Czas\tHeun\tMidpoint\tRK4\tAnalityczne\n";

    for (size_t i = 0; i < heun_result.size(); ++i) {
        double t = i * dt;
        double analytical = T0 * std::exp(0.1 * t);
        std::cout << t << "\t"
            << heun_result[i] << "\t"
            << midpoint_result[i] << "\t"
            << rk4_result[i] << "\t"
            << analytical << "\n";
    }

    // Przykład 2: Schładzanie według prawa Newtona
    std::cout << "\nSchladzanie wedlug prawa Newtona dT/dt = -0.2*(T-20)\n";
    auto cooling = [](double t, double T) { return -0.2 * (T - 20.0); };
    double T0_cooling = 80.0;
    double t_max_cooling = 15.0;
    double dt_cooling = 1.0;

    std::vector<double> cooling_heun = ode::heun(cooling, T0_cooling, t_max_cooling, dt_cooling);
    std::vector<double> cooling_rk4 = ode::rk4(cooling, T0_cooling, t_max_cooling, dt_cooling);

    std::cout << "Czas\tHeun\tRK4\tAnalityczne\n";
    for (size_t i = 0; i < cooling_heun.size(); ++i) {
        double t = i * dt_cooling;
        double analytical = 20.0 + (T0_cooling - 20.0) * std::exp(-0.2 * t);
        std::cout << t << "\t"
            << cooling_heun[i] << "\t"
            << cooling_rk4[i] << "\t"
            << analytical << "\n";
    }

    // === NONLINEAR: rozwiązywanie równań nieliniowych ===
    std::cout << "\n=== Rozwiazywanie rownan nieliniowych ===\n";

    int max_iter = 100;
    double tol = 1e-10;

    // Przykład 1: f(x) = x^3 - 2x - 5 = 0
    // Pierwiastek około x = 2.094...
    std::cout << "Rownanie: x^3 - 2x - 5 = 0\n";

    auto f1 = [](double x) { return x * x * x - 2 * x - 5; };
    auto df1 = [](double x) { return 3 * x * x - 2; };  // pochodna

    // Metoda Newtona
    double newton_result1 = nonlinear::newton(f1, df1, 2.0, max_iter, tol);
    std::cout << "Newton (x0=2.0):     " << newton_result1 << ", f(x) = " << f1(newton_result1) << "\n";

    // Metoda siecznych
    double secant_result1 = nonlinear::secant(f1, 1.5, 2.5, max_iter, tol);
    std::cout << "Secant (x0=1.5, x1=2.5): " << secant_result1 << ", f(x) = " << f1(secant_result1) << "\n";

    // Metoda bisekcji
    double bisection_result1 = nonlinear::bisection(f1, 2.0, 3.0, max_iter, tol);
    std::cout << "Bisection [2,3]:     " << bisection_result1 << ", f(x) = " << f1(bisection_result1) << "\n";

    // Metoda falsi
    double falsi_result1 = nonlinear::false_position(f1, 2.0, 3.0, max_iter, tol);
    std::cout << "False position [2,3]: " << falsi_result1 << ", f(x) = " << f1(falsi_result1) << "\n";

    // Przykład 2: f(x) = cos(x) - x = 0
    // Pierwiastek około x = 0.739...
    std::cout << "\nRownanie: cos(x) - x = 0\n";

    auto f2 = [](double x) { return std::cos(x) - x; };
    auto df2 = [](double x) { return -std::sin(x) - 1; };  // pochodna

    double newton_result2 = nonlinear::newton(f2, df2, 0.5, max_iter, tol);
    std::cout << "Newton (x0=0.5):     " << newton_result2 << ", f(x) = " << f2(newton_result2) << "\n";

    double secant_result2 = nonlinear::secant(f2, 0.0, 1.0, max_iter, tol);
    std::cout << "Secant (x0=0.0, x1=1.0): " << secant_result2 << ", f(x) = " << f2(secant_result2) << "\n";

    double bisection_result2 = nonlinear::bisection(f2, 0.0, 1.0, max_iter, tol);
    std::cout << "Bisection [0,1]:     " << bisection_result2 << ", f(x) = " << f2(bisection_result2) << "\n";

    // Przykład 3: f(x) = x^2 - 2 = 0 (pierwiastek z 2)
    // Pierwiastki: x = ±√2 ≈ ±1.414...
    std::cout << "\nRownanie: x^2 - 2 = 0 (pierwiastek z 2)\n";

    auto f3 = [](double x) { return x * x - 2; };
    auto df3 = [](double x) { return 2 * x; };

    // Znajdowanie dodatniego pierwiastka
    double newton_result3 = nonlinear::newton(f3, df3, 1.5, max_iter, tol);
    std::cout << "Newton (x0=1.5):     " << newton_result3 << ", sqrt(2) = " << std::sqrt(2) << "\n";

    double secant_result3 = nonlinear::secant(f3, 1.0, 2.0, max_iter, tol);
    std::cout << "Secant (x0=1.0, x1=2.0): " << secant_result3 << "\n";

    // Znajdowanie ujemnego pierwiastka
    double newton_result3_neg = nonlinear::newton(f3, df3, -1.5, max_iter, tol);
    std::cout << "Newton ujemny (x0=-1.5): " << newton_result3_neg << ", -sqrt(2) = " << -std::sqrt(2) << "\n";

    // Przykład 4: f(x) = e^x - 3x = 0
    // Ma dwa pierwiastki
    std::cout << "\nRownanie: e^x - 3x = 0\n";

    auto f4 = [](double x) { return std::exp(x) - 3 * x; };
    auto df4 = [](double x) { return std::exp(x) - 3; };

    // Pierwszy pierwiastek (około x = 0.619)
    double newton_result4a = nonlinear::newton(f4, df4, 0.5, max_iter, tol);
    std::cout << "Newton (x0=0.5):     " << newton_result4a << ", f(x) = " << f4(newton_result4a) << "\n";

    // Drugi pierwiastek (około x = 1.512)
    double newton_result4b = nonlinear::newton(f4, df4, 1.8, max_iter, tol);
    std::cout << "Newton (x0=1.8):     " << newton_result4b << ", f(x) = " << f4(newton_result4b) << "\n";

    // Używamy bisekcji dla pierwszego pierwiastka
    double bisection_result4 = nonlinear::bisection(f4, 0.0, 1.0, max_iter, tol);
    std::cout << "Bisection [0,1]:     " << bisection_result4 << ", f(x) = " << f4(bisection_result4) << "\n";

    // Przykład 5: Równanie styczne - f(x) = x^3 - x = 0
    // Pierwiastki: x = -1, 0, 1
    std::cout << "\nRownanie: x^3 - x = 0 (pierwiastki: -1, 0, 1)\n";

    auto f5 = [](double x) { return x * x * x - x; };
    auto df5 = [](double x) { return 3 * x * x - 1; };

    // Znajdowanie wszystkich pierwiastków
    double root_neg = nonlinear::newton(f5, df5, -0.8, max_iter, tol);
    double root_zero = nonlinear::bisection(f5, -0.5, 0.5, max_iter, tol);
    double root_pos = nonlinear::newton(f5, df5, 0.8, max_iter, tol);

    std::cout << "Pierwiastek ujemny:   " << root_neg << "\n";
    std::cout << "Pierwiastek zero:     " << root_zero << "\n";
    std::cout << "Pierwiastek dodatni:  " << root_pos << "\n";

    return 0;
}
