#pragma once
#include <fstream>
#include <vector>
#include <string>
void inputGrid(std::vector<std::vector<std::vector<double>>>& grid, std::string str, std::vector<int>& gridSizes) {
    double x, y;
    int rows, cols;

    std::ifstream ifs(str);
    ifs >> rows >> cols;
    gridSizes.push_back(rows);
    gridSizes.push_back(cols);
    for (int i = 0; i < rows; ++i) {
        std::vector<std::vector<double>> row;
        for (int j = 0; j < cols; ++j) {
            std::vector<double> dots;
            ifs >> x >> y;
            dots.push_back(x);
            dots.push_back(y);
            row.push_back(dots);
        }
        grid.push_back(row);
    }
    ifs.close();
}
