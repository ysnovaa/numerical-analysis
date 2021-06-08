#include <cstdio>
#include <chrono>
#include <vector>
#include <bits/stdc++.h>
#include <iostream>
#include <fstream>
#include "readFile.hpp"

using namespace std;
using namespace std::chrono;

struct hash_pair
{
    template <class T1, class T2>
    size_t operator()(const pair<T1, T2> &p) const
    {
        auto hash1 = hash<T1>{}(p.first);
        auto hash2 = hash<T2>{}(p.second);
        return hash1 ^ hash2;
    }
};

#define sparse_matrix unordered_map<pair<int, int>, double, hash_pair>

sparse_matrix to_sparse(vector<vector<double>> &matrix)
{
    sparse_matrix S;
    for (int i = 0; i < matrix.size(); ++i)
    {
        for (int j = 0; j < matrix[i].size(); ++j)
        {
            if (matrix[i][j])
                S[make_pair(i, j)] = matrix[i][j];
        }
    }
    return S;
}

sparse_matrix multiply(sparse_matrix &A, sparse_matrix &B)
{
    sparse_matrix C;
    for (auto &keyA : A)
    {
        for (auto &keyB : B)
        {
            if (keyA.first.second == keyB.first.first)
            {
                int x = (int)keyA.first.first, y = (int)keyB.first.second;
                C[make_pair(x, y)] += keyA.second * keyB.second;
            }
        }
    }
    return C;
}

void recordTimes(string filename)
{
    string directory = "../matrices/" + filename;
    std::ofstream myfile;
    myfile.open("../benchmarks.csv", std::ios_base::app); // append instead of overwrite
    vector<vector<double>> matrix = fill_matrix(directory);
    vector<vector<double>> res(matrix.size(), vector<double>(matrix[0].size(), 0.));

    sparse_matrix sparse = to_sparse(matrix);
    for (int i = 0; i < 5; ++i)
    {
        auto start = high_resolution_clock::now();
        multiply(sparse, sparse);
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        myfile << "cpp,sparse," << filename << "," << res.size() << "," << duration.count() << "\n";
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

//! g++ sparse.cpp readFile.cpp -o out
// https://www.youtube.com/watch?v=3cdeehULUBo