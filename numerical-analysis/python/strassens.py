# credit: https://github.com/jaeho3690/Matrix_multiplication_python/blob/master/strassen_base.py
"""20161349 JaehoKim
Intro to Algorithms Strassen"""

# Import Libray. I only used numpy for numpy array. None of the functions implemented in numpy has been used during matrix calculations
import numpy as np
from fillMatrix import recordTimes


def strassen(array_first, array_second):
    n = len(array_first)
    if n == 1:
        return array_first * array_second
    else:
        a11 = array_first[:int(len(array_first)/2), :int(len(array_first)/2)]
        a12 = array_first[:int(len(array_first)/2), int(len(array_first)/2):]
        a21 = array_first[int(len(array_first)/2):, :int(len(array_first)/2)]
        a22 = array_first[int(len(array_first)/2):, int(len(array_first)/2):]

        b11 = array_second[:int(len(array_second)/2),
                           :int(len(array_second)/2)]
        b12 = array_second[:int(len(array_second)/2),
                           int(len(array_second)/2):]
        b21 = array_second[int(len(array_second)/2):,
                           :int(len(array_second)/2)]
        b22 = array_second[int(len(array_second)/2):,
                           int(len(array_second)/2):]

        S1 = b12 - b22
        S2 = a11 + a12
        S3 = a21 + a22
        S4 = b21 - b11
        S5 = a11 + a22
        S6 = b11 + b22
        S7 = a12 - a22
        S8 = b21 + b22
        S9 = a11 - a21
        S10 = b11 + b12

        P1 = strassen(a11, S1)
        P2 = strassen(S2, b22)
        P3 = strassen(S3, b11)
        P4 = strassen(a22, S4)
        P5 = strassen(S5, S6)
        P6 = strassen(S7, S8)
        P7 = strassen(S9, S10)

        c11 = P5 + P4 - P2 + P6
        c12 = P1 + P2
        c21 = P3 + P4
        c22 = P5 + P1 - P3 - P7

        result = np.zeros((n, n))
        result[:int(len(result)/2), :int(len(result)/2)] = c11
        result[:int(len(result)/2), int(len(result)/2):] = c12
        result[int(len(result)/2):, :int(len(result)/2)] = c21
        result[int(len(result)/2):, int(len(result)/2):] = c22
        return result


if __name__ == "__main__":

    # recordTimes("Hamrle1.mtx", strassen, "strassen")
    # recordTimes("GD99_b.mtx", strassen, "strassen")
    # recordTimes("can_256.mtx", strassen, "strassen")
    # recordTimes("dwa512.mtx", strassen, "strassen")
    recordTimes("delaunay_n10.mtx", strassen, "strassen")
