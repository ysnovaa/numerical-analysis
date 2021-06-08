from fillMatrix import recordTimes


def naive(A, B):
    result = [[sum(a * b for a, b in zip(A_row, B_col))
               for B_col in zip(*B)]
              for A_row in A]
    return result


if __name__ == '__main__':

    recordTimes("Hamrle1.mtx", naive, "naive")
    recordTimes("GD99_b.mtx", naive, "naive")
    recordTimes("can_256.mtx", naive, "naive")
    recordTimes("dwa512.mtx", naive, "naive")
    recordTimes("delaunay_n10.mtx", naive, "naive")
