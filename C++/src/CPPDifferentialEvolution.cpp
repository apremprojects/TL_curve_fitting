#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include "../include/CPPDifferentialEvolution.h"
#include "../include/DifferentialEvolution.h"
std::vector<double> subtract(const std::vector<double>& v1, const std::vector<double>& v2) {
    assert(v1.size() == v2.size());
    std::vector<double> v3(v1.size(), 0.0);
    for (int i = 0; i < v3.size(); i++) {
        v3[i] = v1[i] - v2[i];
    }
    return v3;
}

double sum_of_squares(const std::vector<double>& v) {
    double sum = 0;
    for (const double i : v) {
        sum += i * i;
    }
    return sum;
}

Model::Model(const int _peaks, const double _beta, const double _k_J_per_K, const double _N, const std::array<double, 6> _bounds, const std::vector<double>& _x_data, const std::vector<double>& _y_data) : peaks(_peaks), beta(_beta), k_J_per_K(_k_J_per_K), N(_N), bounds(_bounds), x_data(_x_data), y_data(_y_data) {
    //std::cout << "Model::Model()\n";
}

double Model::EvaluteCost(std::vector<double> inputs) const {
    //residual
    assert(inputs.size() == NumberOfParameters());
    std::vector<double> res = evaluate(inputs, x_data);
    return sum_of_squares(subtract(y_data, res));
}

std::vector<de::IOptimizable::Constraints> Model::GetConstraints() const {
    std::vector<Constraints> constr(NumberOfParameters());
    for (int i = 0; i < NumberOfParameters(); i += 3)
    {
        constr[i] = Constraints(bounds[0], bounds[1], true);
        constr[i + 1] = Constraints(bounds[2], bounds[3], true);
        constr[i + 2] = Constraints(bounds[4], bounds[5], true);
    }
    return constr;
}

std::vector<double> Model::evaluate(const std::vector<double>& inputs, const std::vector<double>& x) const {
    //std::cout << "Input -> " << inputs.size() << ", " << x.size() << "\n";
    assert(inputs.size() == NumberOfParameters());
    std::vector<double> res(x.size(), 0.0);
    for (int i = 0; i < inputs.size(); i += 3) {
        double E_ev = inputs[i];
        double S = inputs[i + 1];
        double n0 = inputs[i + 2];
        double E_J = E_ev * (1.602e-19);
        double T0_antiderivative = antiderivative(x[0], E_J);
        for (int j = 0; j < x.size(); j++) {
            double integral_result = antiderivative(x[j], E_J) - T0_antiderivative;
            double numerator = (n0 * n0) * S * exp(-E_J / (k_J_per_K * x[j]));
            double denominator = N * (1 + (n0 * S) / (beta * N) * integral_result) * (1 + (n0 * S) / (beta * N) * integral_result);
            res[j] += numerator / denominator;
        }
    }
    return res;
}