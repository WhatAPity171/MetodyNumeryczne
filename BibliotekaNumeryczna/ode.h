#pragma once
#include "numlib.h"
namespace ode {
    // f: funkcja pochodnej, T0: wartoœæ pocz¹tkowa, t_max: koniec przedzia³u, dt: krok czasu
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