#pragma once
#include "DifferentialEvolution.h"
#include "CPPDifferentialEvolution.h"

extern "C"{
    double* solve(bool verbose, const int _peaks, const int maxiter, const double beta, const double atol, const double tol, const int popsize, double* _bounds, double* _x_data, double* _y_data, size_t _xy_size);
}