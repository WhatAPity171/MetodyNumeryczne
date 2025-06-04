#include "numlib.h"

//=================== Uk³ady liniowe ===================//
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
                }
                else U[i][j] = A[i][j];
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