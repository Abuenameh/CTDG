/* 
 * File:   casadiri.hpp
 * Author: Abuenameh
 *
 * Created on November 29, 2015, 11:11 PM
 */

#ifndef CASADIMATH_HPP
#define	CASADIMATH_HPP

#include "gutzwiller.hpp"

inline double eps(vector<double>& U, int i, int j, int n, int m) {
	return n * U[i] - (m - 1) * U[j];
}

SX energy(SX& fin, SX& J, SX& U0, SX& dU, double mu);
SX energy(SX& fin, SX& J, SX& U0, SX& dU, SX& mu);
SX canonical(SX& fin, SX& J, SX& U0, SX& dU, SX mu);

#include "casadimath.hincl"

#endif	/* CASADIRI_HPP */

