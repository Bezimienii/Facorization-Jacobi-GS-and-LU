def copyvec(A: list):
    res = []
    for i in range(len(A)):
        res.append(A[i])
    return res

def copymat(A: list):
    res = []
    for i in range(len(A)):
        res.append([])
        for j in range(len(A[i])):
            res[i].append(A[i][j])
    return res

def dotpro(A: list, x: list):
    final_vec = []
    for i in range(len(A)):
        final_vec.append(0)
        for j in range(len(x)):
            final_vec[i] += A[i][j] * x[j]
    return final_vec


def subvecs(x: list, b: list):
    cpyx = copyvec(x)
    for i in range(len(x)):
        cpyx[i] -= b[i]
    return cpyx
