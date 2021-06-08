include("./sparse.jl")
include("./strassens.jl")
include("./naive.jl")
using MatrixMarket
using CSV
using DataFrames


function writeTimes(n, file, func, functionName)
    dir = "../matrices/"
    filename = dir * file
    matrix = MatrixMarket.mmread(filename)
    s = size(matrix)[1]
    for _ in 1:n
        t = @timed func(matrix, matrix)
        millis = t[2] * 1000
        dataframe = DataFrame(language=["julia"], algorithm=[functionName], 
        filename=[file], 
        matrix_size=[s], 
        time=[millis])
        CSV.write("../benchmarks.csv", dataframe, append=true)
    end
end

function runBenchmarks()
    # naive
    writeTimes(5, "Hamrle1.mtx", naiveMultiply, "naive")
    writeTimes(5, "GD99_b.mtx", naiveMultiply, "naive")
    writeTimes(5, "can_256.mtx", naiveMultiply, "naive")
    writeTimes(5, "dwa512.mtx", naiveMultiply, "naive")
    writeTimes(5, "delaunay_n10.mtx", naiveMultiply, "naive")
    # # strassens
    writeTimes(5, "Hamrle1.mtx", strassen, "strassen")
    writeTimes(5, "GD99_b.mtx", strassen, "strassen")
    writeTimes(5, "can_256.mtx", strassen, "strassen")
    writeTimes(5, "dwa512.mtx", strassen, "strassen")
    writeTimes(5, "delaunay_n10.mtx", strassen, "strassen")
    # dictionary of keys
    writeTimes(5, "Hamrle1.mtx", sparseMultiplication, "sparse")
    writeTimes(5, "GD99_b.mtx", sparseMultiplication, "sparse")
    writeTimes(5, "can_256.mtx", sparseMultiplication, "sparse")
    writeTimes(5, "dwa512.mtx", sparseMultiplication, "sparse")
    writeTimes(5, "delaunay_n10.mtx", sparseMultiplication, "sparse")
end

runBenchmarks()