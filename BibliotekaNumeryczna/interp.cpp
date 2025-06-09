#include "interp.h"

namespace interp {
    constexpr int SKIP = 3;

    double lagrange(int x, int N, const double* arrX, const double* arrY) {
        double answer = 0;
        for (int i = 0; i < N - 1; i++) {
            if (i % SKIP != 0) continue;
            double Li = 1;
            for (int j = 0; j < N - 1; j++) {
                if (j % SKIP != 0 || j == i) continue;
                Li *= (arrX[x] - arrX[j]) / (arrX[i] - arrX[j]);
            }
            answer += arrY[i] * Li;
        }
        return answer;
    }
}
