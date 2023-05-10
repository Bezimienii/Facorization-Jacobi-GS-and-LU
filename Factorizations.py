import matplotlib.pyplot as pp
import math
import time
import forMatrices as fm


def norma(wektor: list):
    norm = 0
    for el in wektor:
        norm += el ** 2
    norm = math.sqrt(norm)
    return norm


def jacobi(A: list, b: list):
    start = time.time()
    x = []
    xk1 = []
    for _ in range(0, len(b)):
        x.append(0)
        xk1.append(0)
    ilerazy = 0
    while True:
        for i in range(0, len(b)):
            xk1[i] = b[i]
            for j in range(len(A[i])):
                if i != j:
                    xk1[i] -= A[i][j]*x[j]
            xk1[i] /= A[i][i]
        dot = fm.dotpro(A, xk1)
        norm = norma(fm.subvecs(dot, b))
        ilerazy += 1
        if norm < (10**(-9)):
            break
        x = xk1
    end = time.time()
    print("Czas Jacobi: ", end-start)
    print("Ilosc iteracji: ", ilerazy)
    return end-start

def gauss_seidel(A: list, b: list):
    start = time.time()
    x = []
    xk1 = []
    for _ in range(0, len(b)):
        x.append(0)
        xk1.append(0)
    ilerazy = 0
    while True:
        for i in range(0, len(b)):
            xk1[i] = b[i]
            for j in range(len(A[i])):
                if i > j:
                    xk1[i] -= (A[i][j] * xk1[j])
                elif i < j:
                    xk1[i] -= (A[i][j] * x[j])
            xk1[i] /= A[i][i]
        dot = fm.dotpro(A, x)
        norm = norma(fm.subvecs(dot, b))
        ilerazy += 1
        if norm < (10 ** (-9)):
            break
        x = xk1
    end = time.time()
    print("Czas Gaussa-Seidla: ", end - start)
    print("Ilosc iteracji: ", ilerazy)
    return end-start


def tableA(N: int, a1: int, a23: int):
    tablica = []
    for i in range(0, N):
        tablica.append([])
        for j in range(0, N):
            if i == j:
                tablica[i].append(a1)
            elif abs(i-j) < 3:
                tablica[i].append(a23)
            else:
                tablica[i].append(0)
    return tablica


def makeb(N: int):
    b = []
    for i in range(0, N):
        b.append(math.sin(5*i))
    return b

def factorizationLU(A: list, b: list):
    start = time.time()

    U = fm.copymat(A)
    L = [[float(0) if i != j else float(1) for j in range(len(A[0]))] for i in range(len(A))]
    P = [[float(0) if i != j else float(1) for j in range(len(A[0]))] for i in range(len(A))]

    for i in range(len(A) - 1):

        ind = 0
        maxval = 0.0

        for j in range(i, len(A)):
            if abs(U[j][i]) > maxval:
                ind = j
                maxval = abs(U[j][i])

        if i != ind:
            for j in range(i, len(A)):
                U[i][j], U[ind][j] = U[ind][j], U[i][j]

            for j in range(0, i):
                L[i][j], L[ind][j] = L[ind][j], L[i][j]

            for j in range(len(A)):
                P[i][j], P[ind][j] = P[ind][j], P[i][j]

        for j in range(i + 1, len(A)):
            L[j][i] = U[j][i] / U[i][i]
            for k in range(i, len(A)):
                U[j][k] -= (U[i][k] * L[j][i])

    pivotedb = fm.dotpro(P, b)

    for i in L:
        print(i)

    y = [j for j in pivotedb]

    for j in range(len(A)):
        for k in range(j):
            y[j] -= L[j][k] * y[k]
        y[j] /= L[j][j]

    x = [j for j in y]

    for j in range(len(A) - 1, -1, -1):
        for k in range(j + 1, len(A)):
            x[j] -= U[j][k] * x[k]
        x[j] /= U[j][j]

    B = [[float(0) for _ in range(len(A[0]))] for _ in range(len(A))]
    for i in range(len(A)):
        for j in range(len(A[0])):
            B[i][j] += L[i][j]*U[j][i]

    res = fm.subvecs(fm.dotpro(B, x), pivotedb)

    print("Metoda LU")
    timehere = time.time() - start
    print("Norma:", norma(res))
    return timehere

def GSExample1():
    A = tableA(943, 10, -1)
    b = makeb(943)
    jacobi(A, b)
    gauss_seidel(A, b)

def GSExample2():
    A = tableA(943, 3, -1)
    b = makeb(943)
    jacobi(A, b)
    gauss_seidel(A, b)

def LUExample1():
    A = tableA(943, 3, -1)
    b = makeb(943)
    factorizationLU(A, b)

def AllFactorizations():
    newtable = [100, 500, 1000, 2000, 3000]
    y1 = []
    y2 = []
    y3 = []
    for i in newtable:
        A = tableA(i, 10, -1)
        b = makeb(i)
        x1 = jacobi(A, b)
        x2 = gauss_seidel(A, b)
        newval = factorizationLU(A, b)
        y1.append(x1)
        y2.append(x2)
        y3.append(newval)

    pp.plot(newtable, y1, label = "Jacobi")
    pp.plot(newtable, y2, label = "GS")
    pp.plot(newtable, y3, label = "LU")
    pp.ylabel('czas [s]')
    pp.xlabel('dlugosc macierzy')
    pp.legend()
    pp.show()
