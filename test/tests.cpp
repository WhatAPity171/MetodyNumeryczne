#include "tests.h"
#include "integral.h"
#include "linear.h"
#include "nonlinear.h"
#include "interp.h"
#include "ode.h"
#include "numlib.h"

#include <iostream>
#include <cmath>
#include <vector>
#include <stdexcept>

namespace tests {

// Prosty helper do testów: porównanie double z tolerancją
bool approxEqual(double a, double b, double tol = 1e-6) {
    return std::fabs(a - b) < tol;
}

void testIntegral() {
    using namespace integral;
    std::cout << "Testing integral functions...\n";

    // Simpson - test 1
    auto f1 = [](double x) { return x * x; };
    double res1 = simpson(f1, 0, 2, 10);
    if (!approxEqual(res1, 8.0/3.0)) std::cout << "[FAIL] simpson test1\n";

    // Simpson - test 2 (throw on odd n)
    bool caught = false;
    try {
        simpson(f1, 0, 1, 3);
    } catch (const std::invalid_argument&) {
        caught = true;
    }
    if (!caught) std::cout << "[FAIL] simpson test2 (no exception)\n";

    // Trapezoid - test 1
    auto f2 = [](double x) { return x; };
    double res2 = trapezoid(f2, 0, 1, 10);
    if (!approxEqual(res2, 0.5)) std::cout << "[FAIL] trapezoid test1\n";

    // Trapezoid - test 2 (integral of constant)
    auto f3 = [](double) { return 3.0; };
    double res3 = trapezoid(f3, 0, 1, 10);
    if (!approxEqual(res3, 3.0)) std::cout << "[FAIL] trapezoid test2\n";

    // Rectangle - test 1
    auto f4 = [](double) { return 1.0; };
    double res4 = rectangle(f4, 0, 1, 100);
    if (!approxEqual(res4, 1.0)) std::cout << "[FAIL] rectangle test1\n";

    // Rectangle - test 2 (linear)
    double res5 = rectangle(f2, 0, 1, 1000);
    if (!approxEqual(res5, 0.5, 1e-3)) std::cout << "[FAIL] rectangle test2\n";

    std::cout << "Integral tests done.\n";
}

void testLinear() {
    using namespace linear;
    std::cout << "Testing linear algebra functions...\n";

    // Gaussian elimination test 1: solve x + y = 2, x - y = 0 => x=1,y=1
    {
        std::vector<std::vector<double>> A = {{1,1},{1,-1}};
        std::vector<double> b = {2,0}, x;
        gaussianElimination(A, b, x);
        if (!(approxEqual(x[0], 1.0) && approxEqual(x[1], 1.0))) 
            std::cout << "[FAIL] gaussianElimination test1\n";
    }

    // Gaussian elimination test 2: singular matrix check (should handle without crash)
    {
        std::vector<std::vector<double>> A = {{1,2},{2,4}};
        std::vector<double> b = {3,6}, x;
        try {
            gaussianElimination(A, b, x);
            // Solution might exist or not, no exception here, so just check result
        } catch (...) {
            std::cout << "[FAIL] gaussianElimination test2 threw exception\n";
        }
    }

    // LU decomposition test 1
    {
        std::vector<std::vector<double>> A = {{4,3},{6,3}};
        std::vector<std::vector<double>> L,U;
        std::vector<int> rowPerm, colPerm;
        LUDecomposition(A,L,U,rowPerm,colPerm);
        // Check L diagonal = 1
        bool diagOK = true;
        for(size_t i=0;i<L.size();++i)
            if(!approxEqual(L[i][i],1.0)) diagOK=false;
        if(!diagOK) std::cout << "[FAIL] LUDecomposition test1 (L diagonal)\n";
    }

    // LU decomposition test 2: singular matrix should throw
    {
        std::vector<std::vector<double>> A = {{0,0},{0,0}};
        std::vector<std::vector<double>> L,U;
        std::vector<int> rowPerm, colPerm;
        bool caught = false;
        try {
            LUDecomposition(A,L,U,rowPerm,colPerm);
        } catch (const std::runtime_error&) {
            caught = true;
        }
        if(!caught) std::cout << "[FAIL] LUDecomposition test2 (no exception for singular)\n";
    }

    std::cout << "Linear algebra tests done.\n";
}

void testNonlinear() {
    using namespace nonlinear;
    std::cout << "Testing nonlinear solvers...\n";

    auto f = [](double x) { return x*x - 2; };
    auto df = [](double x) { return 2*x; };

    // Newton test 1 (root sqrt(2))
    double root1 = newton(f, df, 1.0);
    if (!approxEqual(root1, std::sqrt(2), 1e-6)) std::cout << "[FAIL] newton test1\n";

    // Newton test 2 (derivative zero case returns NAN)
    auto dfzero = [](double) { return 0.0; };
    double res = newton(f, dfzero, 1.0);
    if (!std::isnan(res)) std::cout << "[FAIL] newton test2 (derivative zero case)\n";

    // Secant test 1
    double root2 = secant(f, 1.0, 2.0);
    if (!approxEqual(root2, std::sqrt(2), 1e-6)) std::cout << "[FAIL] secant test1\n";

    // Secant test 2 (fail on flat function)
    auto flat = [](double) { return 1.0; };
    double fail = secant(flat, 0.0, 1.0);
    if (!std::isnan(fail)) std::cout << "[FAIL] secant test2 (flat function)\n";

    // Bisection test 1
    double root3 = bisection(f, 0, 2);
    if (!approxEqual(root3, std::sqrt(2), 1e-6)) std::cout << "[FAIL] bisection test1\n";

    // Bisection test 2 (invalid interval returns NAN)
    double fail2 = bisection(f, 2, 3);
    if (!std::isnan(fail2)) std::cout << "[FAIL] bisection test2 (invalid interval)\n";

    // False position test 1
    double root4 = false_position(f, 0, 2);
    if (!approxEqual(root4, std::sqrt(2), 1e-6)) std::cout << "[FAIL] false_position test1\n";

    // False position test 2 (invalid interval returns NAN)
    double fail3 = false_position(f, 2, 3);
    if (!std::isnan(fail3)) std::cout << "[FAIL] false_position test2 (invalid interval)\n";

    std::cout << "Nonlinear solver tests done.\n";
}

void testInterp() {
    using namespace interp;
    std::cout << "Testing interpolation...\n";

    // Setup sample points
    const int N = 5;
    double arrX[N] = {0,1,2,3,4};
    double arrY[N] = {0,1,4,9,16};

    // Lagrange test 1 (interpolate at x=2, should be 4)
    double val1 = lagrange(2, N, arrX, arrY);
    if (!approxEqual(val1, 4.0)) std::cout << "[FAIL] lagrange test1\n";

    // Lagrange test 2 (interpolate at x=3, should be 9)
    double val2 = lagrange(3, N, arrX, arrY);
    if (!approxEqual(val2, 9.0)) std::cout << "[FAIL] lagrange test2\n";

    std::cout << "Interpolation tests done.\n";
}

void testOde() {
    using namespace ode;
    std::cout << "Testing ODE solvers...\n";

    auto f = [](double t, double y) { return y; }; // dy/dt = y, solution: y = e^t

    // Heun test 1: y(0)=1, t=1, dt=0.1
    auto sol1 = heun(f, 1.0, 1.0, 0.1);
    if (!approxEqual(sol1.back(), std::exp(1.0), 1e-2)) std::cout << "[FAIL] heun test1\n";

    // Heun test 2: check length of output vector
    if (sol1.size() != 11) std::cout << "[FAIL] heun test2 (length)\n";

    // Midpoint test 1
    auto sol2 = midpoint(f, 1.0, 1.0, 0.1);
    if (!approxEqual(sol2.back(), std::exp(1.0), 1e-2)) std::cout << "[FAIL] midpoint test1\n";

    // Midpoint test 2 (length check)
    if (sol2.size() != 11) std::cout << "[FAIL] midpoint test2 (length)\n";

    // RK4 test 1
    auto sol3 = rk4(f, 1.0, 1.0, 0.1);
    if (!approxEqual(sol3.back(), std::exp(1.0), 1e-4)) std::cout << "[FAIL] rk4 test1\n";

    // RK4 test 2 (length check)
    if (sol3.size() != 11) std::cout << "[FAIL] rk4 test2 (length)\n";

    std::cout << "ODE solver tests done.\n";
}

void testNumlib() {
    std::cout << "Testing numlib utils...\n";
    // printMatrix and printVector just print, no return, so minimal test:
    std::vector<std::vector<double>> mat = {{1,2},{3,4}};
    std::vector<double> vec = {5,6};

    printMatrix(mat, "Matrix");
    printVector(vec, "Vector");

    int N;
    std::vector<std::vector<double>> A;
    std::vector<double> b;

    bool loaded = loadMatrixFile("non_existent_file.txt", A, b, N);
    if (loaded) std::cout << "[FAIL] loadMatrixFile test1 (should fail for missing file)\n";

    std::cout << "numlib tests done.\n";
}

void runAll() {
    testIntegral();
    testLinear();
    testNonlinear();
    testInterp();
    testOde();
    testNumlib();
    std::cout << "All tests finished.\n";
}

} // namespace tests

int main() {
    tests::runAll();
    return 0;
}
