#include <cstdio>
#include <chrono>
#include <vector>
#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include "readFile.hpp"

using namespace std;
using namespace std::chrono;

void naive_multiply(vector<vector<double>> &A, vector<vector<double>> &B, vector<vector<double>> R)
{
    for (int i = 0; i < A.size(); i++)
    {
        for (int j = 0; j < B.size(); j++)
        {
            for (int k = 0; k < R.size(); k++)
            {
                R[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void recordTimes(string filename)
{
    string directory = "../matrices/" + filename;
    std::ofstream myfile;
    myfile.open("../benchmarks.csv", std::ios_base::app); // append instead of overwrite
    vector<vector<double>> matrix = fill_matrix(directory);
    vector<vector<double>> res(matrix.size(), vector<double>(matrix[0].size(), 0.));

    for (int i = 0; i < 5; i++)
    {
        auto start = high_resolution_clock::now();
        naive_multiply(matrix, matrix, res);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        myfile << "cpp,naive," << filename << "," << res.size() << "," << duration.count() << "\n";
    }
    myfile.close();
}

int main()
{
    recordTimes("Hamrle1.mtx");
    recordTimes("GD99_b.mtx");
    recordTimes("can_256.mtx");
    recordTimes("dwa512.mtx");
    recordTimes("delaunay_n10.mtx");
}

//! g++ naive.cpp readFile.cpp -o out
