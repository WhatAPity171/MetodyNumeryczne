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



//=================== Wspólne narzêdzia ===================//
constexpr double EPS = 1e-12;

inline void printMatrix(const std::vector<std::vector<double>>& matrix, const std::string& name = "") {
    if (!name.empty()) std::cout << name << ":\n";
    for (const auto& row : matrix) {
        for (double val : row) std::cout << std::setw(10) << std::fixed << std::setprecision(4) << val << " ";
        std::cout << "\n";
    }
    std::cout << "\n";
}

inline void printVector(const std::vector<double>& v, const std::string& name = "") {
    if (!name.empty()) std::cout << name << " = [ ";
    for (double val : v) std::cout << std::fixed << std::setprecision(4) << val << " ";
    std::cout << "]\n";
}

inline bool loadMatrixFile(const std::string& filename, std::vector<std::vector<double>>& A, std::vector<double>& b, int& N) {
    std::ifstream file(filename);
    if (!file) return false;

    std::string temp;
    file >> temp >> N;
    b.resize(N);
    A.assign(N, std::vector<double>(N));

    file >> temp;
    for (int i = 0; i < N; i++) file >> b[i];

    file >> temp;
    for (int i = 0; i < N; i++)
        for (int j = 0; j < N; j++)
            file >> A[i][j];

    return true;
}

