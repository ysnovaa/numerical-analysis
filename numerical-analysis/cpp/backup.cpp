#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <bits/stdc++.h>
#include <chrono>
#include "readFile.hpp"
// credit: https://github.com/MartinThoma/matrix-multiplication/blob/master/C%2B%2B/strassen-algorithm.cpp

// Set LEAF_SIZE to 1 if you want to the pure strassen algorithm
// otherwise, the ikj-algorithm will be applied when the split
// matrices are as small as LEAF_SIZE x LEAF_SIZE
int leafsize;

using namespace std;
using namespace std::chrono;
/*
 * Implementation of the strassen algorithm, similar to
 * http://en.wikipedia.org/w/index.php?title=Strassen_algorithm&oldid=498910018#Source_code_of_the_Strassen_algorithm_in_C_language
 */

void strassen(vector<vector<double>> &A,
              vector<vector<double>> &B,
              vector<vector<double>> &C, unsigned int tam);
// unsigned int nextPowerOfTwo(double n);
void strassenR(vector<vector<double>> &A,
               vector<vector<double>> &B,
               vector<vector<double>> &C,
               int tam);
void sum(vector<vector<double>> &A,
         vector<vector<double>> &B,
         vector<vector<double>> &C, int tam);
// void subtract(vector<vector<double>> &A,
//               vector<vector<double>> &B,
//               vector<vector<double>> &C, int tam);

void ikjalgorithm(vector<vector<double>> A,
                  vector<vector<double>> B,
                  vector<vector<double>> &C, int n)
{
    for (int i = 0; i < n; i++)
    {
        for (int k = 0; k < n; k++)
        {
            for (int j = 0; j < n; j++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void sum(vector<vector<double>> &A,
         vector<vector<double>> &B,
         vector<vector<double>> &C, int tam)
{
    int i, j;

    for (i = 0; i < tam; i++)
    {
        for (j = 0; j < tam; j++)
        {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}

unsigned int customNextPowerOfTwo(int n)
{
    return pow(2, int(ceil(log2(n))));
}

void subtract(vector<vector<double>> &A,
              vector<vector<double>> &B,
              vector<vector<double>> &C, int tam)
{
    int i, j;

    for (i = 0; i < tam; i++)
    {
        for (j = 0; j < tam; j++)
        {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}
void strassenR(vector<vector<double>> &A,
               vector<vector<double>> &B,
               vector<vector<double>> &C, int tam)
{
    if (tam <= leafsize)
    {
        ikjalgorithm(A, B, C, tam);
        return;
    }

    // other cases are treated here:
    else
    {
        int newTam = tam / 2;
        vector<double> inner(newTam);
        vector<vector<double>>
            a11(newTam, inner), a12(newTam, inner), a21(newTam, inner), a22(newTam, inner),
            b11(newTam, inner), b12(newTam, inner), b21(newTam, inner), b22(newTam, inner),
            c11(newTam, inner), c12(newTam, inner), c21(newTam, inner), c22(newTam, inner),
            p1(newTam, inner), p2(newTam, inner), p3(newTam, inner), p4(newTam, inner),
            p5(newTam, inner), p6(newTam, inner), p7(newTam, inner),
            aResult(newTam, inner), bResult(newTam, inner);

        int i, j;

        //dividing the matrices in 4 sub-matrices:
        for (i = 0; i < newTam; i++)
        {
            for (j = 0; j < newTam; j++)
            {
                a11[i][j] = A[i][j];
                a12[i][j] = A[i][j + newTam];
                a21[i][j] = A[i + newTam][j];
                a22[i][j] = A[i + newTam][j + newTam];

                b11[i][j] = B[i][j];
                b12[i][j] = B[i][j + newTam];
                b21[i][j] = B[i + newTam][j];
                b22[i][j] = B[i + newTam][j + newTam];
            }
        }

        // Calculating p1 to p7:

        sum(a11, a22, aResult, newTam);          // a11 + a22
        sum(b11, b22, bResult, newTam);          // b11 + b22
        strassenR(aResult, bResult, p1, newTam); // p1 = (a11+a22) * (b11+b22)

        sum(a21, a22, aResult, newTam);      // a21 + a22
        strassenR(aResult, b11, p2, newTam); // p2 = (a21+a22) * (b11)

        subtract(b12, b22, bResult, newTam); // b12 - b22
        strassenR(a11, bResult, p3, newTam); // p3 = (a11) * (b12 - b22)

        subtract(b21, b11, bResult, newTam); // b21 - b11
        strassenR(a22, bResult, p4, newTam); // p4 = (a22) * (b21 - b11)

        sum(a11, a12, aResult, newTam);      // a11 + a12
        strassenR(aResult, b22, p5, newTam); // p5 = (a11+a12) * (b22)

        subtract(a21, a11, aResult, newTam);     // a21 - a11
        sum(b11, b12, bResult, newTam);          // b11 + b12
        strassenR(aResult, bResult, p6, newTam); // p6 = (a21-a11) * (b11+b12)

        subtract(a12, a22, aResult, newTam);     // a12 - a22
        sum(b21, b22, bResult, newTam);          // b21 + b22
        strassenR(aResult, bResult, p7, newTam); // p7 = (a12-a22) * (b21+b22)

        // calculating c21, c21, c11 e c22:

        sum(p3, p5, c12, newTam); // c12 = p3 + p5
        sum(p2, p4, c21, newTam); // c21 = p2 + p4

        sum(p1, p4, aResult, newTam);       // p1 + p4
        sum(aResult, p7, bResult, newTam);  // p1 + p4 + p7
        subtract(bResult, p5, c11, newTam); // c11 = p1 + p4 - p5 + p7

        sum(p1, p3, aResult, newTam);       // p1 + p3
        sum(aResult, p6, bResult, newTam);  // p1 + p3 + p6
        subtract(bResult, p2, c22, newTam); // c22 = p1 + p3 - p2 + p6

        // Grouping the results obtained in a single matrix:
        for (i = 0; i < newTam; i++)
        {
            for (j = 0; j < newTam; j++)
            {
                C[i][j] = c11[i][j];
                C[i][j + newTam] = c12[i][j];
                C[i + newTam][j] = c21[i][j];
                C[i + newTam][j + newTam] = c22[i][j];
            }
        }
    }
}

void strassen(vector<vector<double>> &A,
              vector<vector<double>> &B,
              vector<vector<double>> &C, unsigned int n)
{
    //unsigned int n = tam;
    unsigned int m = customNextPowerOfTwo(n);
    vector<double> inner(m);
    vector<vector<double>> APrep(m, inner), BPrep(m, inner), CPrep(m, inner);

    for (unsigned int i = 0; i < n; i++)
    {
        for (unsigned int j = 0; j < n; j++)
        {
            APrep[i][j] = A[i][j];
            BPrep[i][j] = B[i][j];
        }
    }

    strassenR(APrep, BPrep, CPrep, m);
    for (unsigned int i = 0; i < n; i++)
    {
        for (unsigned int j = 0; j < n; j++)
        {
            C[i][j] = CPrep[i][j];
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

    for (int i = 0; i < 1; i++)
    {
        auto start = high_resolution_clock::now();
        strassen(matrix, matrix, res, matrix.size());
        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<milliseconds>(stop - start);
        myfile << "cpp,strassen," << filename << "," << res.size() << "," << duration.count() << "\n";
    }
    myfile.close();
}

int main()
{
    recordTimes("can_256.mtx");
    // recordTimes("bcspwr05.mtx");
    // recordTimes("bcspwr07.mtx");
    // recordTimes("bcspwr08.mtx");
    // recordTimes("bcsstk08.mtx");
}
// ! https://bitsploit.blogspot.com/2017/07/optimizing-matrix-multiplication-using.html
//! g++ strassen.cpp readFile.cpp -o out
