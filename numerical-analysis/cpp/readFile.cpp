#include <vector>
#include <iostream>
#include <fstream>
#include <bits/stdc++.h>
#include "readFile.hpp"

using namespace std;

vector<vector<double>> fill_matrix(string filename)
{
    std::ifstream file(filename);
    int num_row, num_col, num_lines;

    // Ignore comments headers
    while (file.peek() == '%')
        file.ignore(2048, '\n');

    // Read number of rows and columns
    file >> num_row >> num_col >> num_lines;

    vector<vector<double>> matrix(num_row, vector<double>(num_col, 0.));

    for (int l = 0; l < num_lines; l++)
    {
        double data;
        int row, col;
        file >> row >> col >> data;
        matrix[row - 1][col - 1] = data;
    }

    file.close();
    return matrix;
}

// void recordTimes()
// {
//     for (int i = 0; i < 1024; i++)
//     {
//     }
// }

// }

// int main()
// {
//     vector<vector<double>> matrix = fill_matrix("../matrices/bcspwr04.mtx");

//     for (const auto &rows : matrix)
//     {
//         for (const auto &el : rows)
//         {
//             cout << el << " ";
//         }
//         cout << '\n';
//     }
// }