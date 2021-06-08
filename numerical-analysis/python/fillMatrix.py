from scipy.io import mmread
import csv
import time


def fillMatrix(filename):
    return mmread(filename).tocsr().toarray()


def recordTimes(filename, function, functionName):
    directory = "../matrices/" + filename
    matrix = fillMatrix(directory)
    with open("../benchmarks.csv", "a") as f:
        writer = csv.writer(f)
        for _ in range(3):
            start = time.time()
            function(matrix, matrix)
            end = time.time()
            writer.writerow([
                "python",
                functionName,
                filename,
                len(matrix),
                (end - start) * 1000
            ])
