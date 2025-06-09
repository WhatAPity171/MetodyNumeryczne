#pragma once
#include <vector>

namespace linear {
    void gaussianElimination(std::vector<std::vector<double>>& A, std::vector<double>& b, std::vector<double>& x);

    void LUDecomposition(const std::vector<std::vector<double>>& input,
        std::vector<std::vector<double>>& L,
        std::vector<std::vector<double>>& U,
        std::vector<int>& rowPerm,
        std::vector<int>& colPerm);

    void forwardSubstitution(const std::vector<std::vector<double>>& L, const std::vector<double>& b, std::vector<double>& z);
    void backwardSubstitution(const std::vector<std::vector<double>>& U, const std::vector<double>& z, std::vector<double>& x);

    void applyPermutation(const std::vector<int>& perm, const std::vector<double>& in, std::vector<double>& out);
    void undoPermutation(const std::vector<int>& perm, const std::vector<double>& in, std::vector<double>& out);
}
