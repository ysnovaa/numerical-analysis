from collections import defaultdict
from fillMatrix import fillMatrix
import time
import csv


class SparseMatrix:

    def __init__(self, rows, cols, S):
        self.rows = rows
        self.cols = cols
        self.S = S

    @staticmethod
    def to_sparse(M):
        S = dict()
        for i, row in enumerate(M):
            for j, val in enumerate(row):
                if val:
                    S[(i, j)] = val
        return S

    @classmethod
    def from_dense(cls, M):
        rows, cols = len(M), len(M[0])
        S = cls.to_sparse(M)
        return cls(rows, cols, S)

    def multiply(self, B):
        A = self.S
        B = B.S
        C = defaultdict(int)
        for (a_r, a_c), a_val in A.items():
            for (b_r, b_c), b_val in B.items():
                if a_c == b_r:
                    C[(a_r, b_c)] += a_val * b_val
        return C


def recordTimes(filename, functionName):
    directory = "../matrices/" + filename
    matrix = fillMatrix(directory)
    matrix = SparseMatrix.from_dense(matrix)
    with open("../benchmarks.csv", "a") as f:
        writer = csv.writer(f)
        for _ in range(5):
            start = time.time()
            matrix.multiply(matrix)
            end = time.time()
            writer.writerow([
                "python",
                functionName,
                filename,
                matrix.rows,
                (end - start) * 1000
            ])


if __name__ == '__main__':

    recordTimes("Hamrle1.mtx", "sparse")
    recordTimes("GD99_b.mtx", "sparse")
    recordTimes("can_256.mtx", "sparse")
    recordTimes("dwa512.mtx", "sparse")
    recordTimes("delaunay_n10.mtx", "sparse")
