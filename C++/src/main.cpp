#include <iostream>
#include <cmath>
#include <vector>
#include <array>
#include "../include/DifferentialEvolution.h"
#include "../include/csv.h"
#include "../include/CPPDifferentialEvolution.h"
#include "../include/CDifferentialEvolution.h"

int main(){
    io::CSVReader<2> in("../../real_data.csv");
    std::cout << "here\n";
    in.read_header(io::ignore_no_column, "temp_c", "intensity_cps");
    double x, y;
    std::vector<double> x_data;
    std::vector<double> y_data;
    while (in.read_row(x, y)) {
        // do stuff with the data
        x_data.push_back(x + 273.15);
        y_data.push_back(y);
        std::cout << x << ", " << y << "\n";
    }
    std::array<double, 6> bounds = {0.5, 2.0, 1.0e6, 1.0e11, 1.0E5, 10.0E5};
    Model model(5, 2.0, 1.381e-23, 1.0e5, bounds, x_data, y_data);
    de::DifferentialEvolution de(model, 100, std::time(nullptr));

    //de.Optimize(1000, false);
    //std::cout << "ANSWER -> " << de.GetBestCost() << "\n";
    std::cout << "ANS2 -> " << solve(true, 8, 3000, 0.0, 0.1, 15, bounds.data(), x_data.data(), y_data.data(), x_data.size())[0] << "\n";
    //test(bounds.data(), 0);
}