#pragma once
#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include "DifferentialEvolution.h"


std::vector<double> subtract(const std::vector<double>& v1, const std::vector<double>& v2);
double sum_of_squares(const std::vector<double>& v);


class Model : public de::IOptimizable
{
public:
    Model(const int _peaks, const double _beta, const double _k_J_per_K, const double _N, const std::array<double, 6> _bounds, const std::vector<double>& _x_data, const std::vector<double>& _y_data);
    double EvaluteCost(std::vector<double> inputs) const override;
    unsigned int NumberOfParameters() const override {return peaks * 3;}
    std::vector<Constraints> GetConstraints() const override;
    double antiderivative(const double T_prime, const double E_J) const {return k_J_per_K / E_J * (T_prime * T_prime) * exp(-E_J / (k_J_per_K * T_prime));}
    std::vector<double> evaluate(const std::vector<double>& inputs, const std::vector<double>& x) const;
private:
    int peaks = 0;
    double k_J_per_K = 0.0;
    double beta = 0.0;
    double N = 0.0;
    std::array<double, 6> bounds;
    std::vector<double> x_data;
    std::vector<double> y_data;
};
