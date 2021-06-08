# credit: https://andrew.gibiansky.com/blog/mathematics/matrix-multiplication/
using SparseArrays
using LinearAlgebra
using SparseMatrixDicts



function naiveMultiply(A, B)
    A = Symmetric(sparse(A), :U)
    B = Symmetric(sparse(B), :U)
    n, m = size(A)
    C = [[0.0 for _ = 1:n] for _ = 1:m]
    for i = 1:n
        for j = 1:n
            for k = 1:n
                C[i][j] += A[i * k] * B[k * j]
            end
        end
    end
    return C
end

# https://discourse.julialang.org/t/how-to-get-the-row-and-column-value-for-sparsematrixcsc-elements/8613/8