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
    std::vector<double> bounds = {
        0.8, 0.9, 1.0e+11, 1.0e+12, 1.0e+03, 1.0e+07, 1.0e+05, 1.0e+07,
        0.8, 0.9, 1.0e+10, 1.0e+11, 1.0e+03, 1.0e+07, 1.0e+05, 1.0e+07,
        0.9, 0.94, 1.0e+10, 1.0e+11, 1.0e+03, 1.0e+07, 1.0e+05, 1.0e+07,
        1.08, 1.12, 1.0e+11, 1.0e+12, 1.0e+03, 1.0e+07, 1.0e+05, 1.0e+07,
        1.18, 1.22, 1.0e+10, 1.0e+11, 1.0e+03, 1.0e+07, 1.0e+05, 1.0e+07,
        1.08, 1.12, 1.0e+08, 1.0e+9, 1.0e+03, 1.0e+07, 1.0e+05, 1.0e+07,
        1.18, 1.22, 1.0e+08, 1.0e+9, 1.0e+03, 1.0e+07, 1.0e+05, 1.0e+07,
        1.38, 1.42, 1.0e+08, 1.0e+9, 1.0e+03, 1.0e+07, 1.0e+05, 1.0e+07
    };
    //Model model(5, 2.0, 1.381e-23, 1.0e5, bounds, x_data, y_data);
    //de::DifferentialEvolution de(model, 100, std::time(nullptr));

    //de.Optimize(1000, false);
    //std::cout << "ANSWER -> " << de.GetBestCost() << "\n";
    std::cout << "ANS2 -> " << solve(true, 8, 3000, 0.0, 1.0, 15, bounds.data(), x_data.data(), y_data.data(), x_data.size())[0] << "\n";
    //test(bounds.data(), 0);
}