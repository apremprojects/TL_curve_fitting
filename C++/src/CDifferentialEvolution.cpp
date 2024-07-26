#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include "../include/CPPDifferentialEvolution.h"
#include "../include/DifferentialEvolution.h"
#include "../include/progressbar.hpp"


void callback(const de::DifferentialEvolution& c) {
}


extern "C" {
    double* solve(bool verbose, const int _peaks, const int maxiter, const double atol, const double tol, const int popsize, double* _bounds, double* _x_data, double* _y_data, size_t _xy_size) {
        std::vector<double> x_data(_xy_size);
        for (int i = 0; i < _xy_size; i++) {
            x_data[i] = _x_data[i];
        }
        std::vector<double> y_data(_xy_size);
        for (int i = 0; i < _xy_size; i++) {
            y_data[i] = _y_data[i];
        }
        std::vector<double> bounds(8 * _peaks);
        for (int i = 0; i < 8 * _peaks; i++) {
            bounds[i] = _bounds[i];
        }
        //std::cout << "Here\n";
        Model model(_peaks, 2.0, 1.381e-23, bounds, x_data, y_data);
        //std::cout << "Here2\n";
        progressbar pbar(100);
        pbar.show_bar(!verbose);
        int count = 0;
        std::function<bool(const de::DifferentialEvolution&)> termcond = [atol, tol,&pbar, &count, maxiter, verbose](const de::DifferentialEvolution& c) {
            std::vector<std::pair<std::vector<double>, double>> pop = c.GetPopulationWithCosts();
            double total_pop_energy = 0.0;
            for (const auto& i : pop) {
                total_pop_energy += i.second;
            }
            double mean = total_pop_energy / double(pop.size());
            double std_dev = 0.0;
            for (const auto& i : pop) {
                std_dev += (i.second - mean) * (i.second - mean);
            }
            double ctol = sqrt(std_dev / double(pop.size()));
            double etol = atol + tol * abs(total_pop_energy / pop.size());
            if (count % (maxiter / 3000) == 0) {
                //update progress bar or display update
                if (verbose) {
                    std::cout << "Current tolerance: " << etol / ctol * 100.0 << "% - Current best cost: " << c.GetBestCost() << " - Current best solution: ";
                    for (const double i : c.GetBestAgent()) {
                        std::cout << i << " ";
                    }
                    std::cout << "\n";
                }
                else {
                }
            }
            pbar.set_progress(etol * 100.0 / ctol);
            pbar.refresh_two();
            count++;
            return ctol <= etol;
        };
        std::cout << "Here3\n";
        de::DifferentialEvolution de(model, popsize * (_peaks * 4), 3, true, callback, termcond);
        std::cout << "Here4\n";
        de.Optimize((maxiter + 1) * popsize * (_peaks * 4), false);
        std::cout << "Here5\n";
        //pbar.show_bar(false);
        std::cout << "\n";
        double* res = (double*)malloc((1 + _peaks * 4) * sizeof(double));
        res[0] = de.GetBestCost();
        std::vector<double> agent = de.GetBestAgent();
        for (int i = 1; i < (1 + _peaks * 4); i++) {
            res[i] = agent[i - 1];
        }
        return res;
    }
}