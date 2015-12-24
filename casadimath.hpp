/* 
 * File:   casadiri.hpp
 * Author: Abuenameh
 *
 * Created on November 29, 2015, 11:11 PM
 */

#ifndef CASADIMATH_HPP
#define	CASADIMATH_HPP

#include <casadi/casadi.hpp>
#include <casadi/solvers/rk_integrator.hpp>
#include <casadi/solvers/collocation_integrator.hpp>
#include <casadi/interfaces/sundials/cvodes_interface.hpp>
#include <casadi/core/function/custom_function.hpp>

using namespace casadi;

#include "gutzwiller.hpp"

namespace casadi {

    inline bool isnan(SX& sx) {
        return sx.at(0).isNan();
    }

    inline bool isinf(SX sx) {
        return sx.at(0).isInf();
    }
}

inline double eps(vector<double>& U, int i, int j, int n, int m) {
	return n * U[i] - (m - 1) * U[j];
}

inline double JW(double W) {
    return alpha * (W * W) / (Ng * Ng + W * W);
}

inline double JWij(double Wi, double Wj) {
    return alpha * (Wi * Wj) / (sqrt(Ng * Ng + Wi * Wi) * sqrt(Ng * Ng + Wj * Wj));
}

inline double UW(double W) {
    return -2 * (g24 * g24) / Delta * (Ng * Ng * W * W) / ((Ng * Ng + W * W) * (Ng * Ng + W * W));
}

inline SX JW(SX W) {
    return alpha * (W * W) / (Ng * Ng + W * W);
}

inline SX JWij(SX Wi, SX Wj) {
    return alpha * (Wi * Wj) / (sqrt(Ng * Ng + Wi * Wi) * sqrt(Ng * Ng + Wj * Wj));
}

inline SX UW(SX W) {
    return -2 * (g24 * g24) / Delta * (Ng * Ng * W * W) / ((Ng * Ng + W * W) * (Ng * Ng + W * W));
}

complex<SX> energy(SX& fin, SX& J, SX& U0, SX& dU, double mu);
complex<SX> energy(SX& fin, SX& J, SX& U0, SX& dU, SX& mu);
complex<SX> energy(int i, SX& fin, SX& J, SX& U0, SX& dU, SX& mu);
complex<SX> energy(int i, int n, SX& fin, SX& J, SX& U0, SX& dU, SX& mu);
complex<SX> canonical(SX& fin, SX& J, SX& U0, SX& dU, SX mu);
complex<SX> canonical(int i, SX& fin, SX& J, SX& U0, SX& dU, SX mu);
complex<SX> canonical(int i, int n, SX& fin, SX& J, SX& U0, SX& dU, SX mu);

#include "casadimath.hincl"

#endif	/* CASADIRI_HPP */

