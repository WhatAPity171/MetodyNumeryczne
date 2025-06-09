#pragma once
#include <vector>
#include <string>

constexpr double EPS = 1e-12;

void printMatrix(const std::vector<std::vector<double>>& matrix, const std::string& name = "");
void printVector(const std::vector<double>& v, const std::string& name = "");
bool loadMatrixFile(const std::string& filename, std::vector<std::vector<double>>& A, std::vector<double>& b, int& N);
