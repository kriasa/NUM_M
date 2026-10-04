import os
import sys

EPS = 1e-12

def read_data(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"'{filename}' не найден")
    
    with open(filename, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f if line.strip()]
        
    if not lines:
        raise ValueError("пусто")
        
    try:
        n = int(lines[0])
        if len(lines) != n + 2:
            raise ValueError(f"недостаточно данных")
            
        A = []
        for i in range(1, n + 1):
            row = list(map(float, lines[i].split()))
            if len(row) != n:
                raise ValueError(f"ошибка в строке {i}: ожидается {n} элементов")
            A.append(row)
            
        b = list(map(float, lines[n + 1].split()))
        if len(b) != n:
            raise ValueError(f"ошибка в векторе b: ожидается {n} элементов")
            
        return n, A, b
        
    except ValueError as e:
        raise ValueError(f"ERROR: {e}")

def print_matrix(title, matrix):
    print(f"\n{title}")
    for row in matrix:
        formatted_row = [f"{0.0 if abs(val) < EPS else val:10.4f}" for val in row]
        print(" ".join(formatted_row))

def print_vector(title, vector):
    print(f"\n{title}")
    for i, val in enumerate(vector):
        clean_val = 0.0 if abs(val) < EPS else val
        print(f"x{i+1} = {clean_val:10.4f}")


def create_zero_matrix(rows, cols):
    return [[0.0 for _ in range(cols)] for _ in range(rows)]

def create_identity_matrix(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

def copy_matrix(A):
    return [row[:] for row in A]

def matmul(A, B):
    n = len(A)
    C = create_zero_matrix(n, n)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                C[i][j] += A[i][k] * B[k][j]
    return C

def matvec(A, b):
    n = len(A)
    res = [0.0] * n
    for i in range(n):
        for j in range(n):
            res[i] += A[i][j] * b[j]
    return res

def matrix_difference(A, B):
    n = len(A)
    max_diff = 0.0
    for i in range(n):
        for j in range(n):
            diff = abs(A[i][j] - B[i][j])
            if diff > max_diff:
                max_diff = diff
    return max_diff


def lu_decompose(A):

    n = len(A)
    U = copy_matrix(A)
    L = create_identity_matrix(n)
    P = create_identity_matrix(n)
    swaps = 0

    for k in range(n):
        pivot_row = k
        max_val = abs(U[k][k])
        for i in range(k + 1, n):
            if abs(U[i][k]) > max_val:
                max_val = abs(U[i][k])
                pivot_row = i

        if max_val < EPS:
            raise ArithmeticError("det(A) = 0")

        if pivot_row != k:
            U[k], U[pivot_row] = U[pivot_row], U[k]
            P[k], P[pivot_row] = P[pivot_row], P[k]
            for j in range(k):
                L[k][j], L[pivot_row][j] = L[pivot_row][j], L[k][j]
            swaps += 1

        for i in range(k + 1, n):
            mu = U[i][k] / U[k][k]
            L[i][k] = mu
            for j in range(k, n):
                U[i][j] -= mu * U[k][j]

    return L, U, P, swaps

def forward_sub(L, Pb):
    # L*y = Pb
    n = len(L)
    y = [0.0] * n
    for i in range(n):
        s = sum(L[i][j] * y[j] for j in range(i))
        y[i] = Pb[i] - s
    return y

def backward_sub(U, y):
    # U*x = y
    n = len(U)
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = sum(U[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (y[i] - s) / U[i][i]
    return x

def solve_system(L, U, P, b):
    # Ax = b
    Pb = matvec(P, b)
    y = forward_sub(L, Pb)
    x = backward_sub(U, y)
    return x

def calculate_inverse(L, U, P):
    # A*A^-1=E
    n = len(L)
    A_inv = create_zero_matrix(n, n)
    for j in range(n):
        e_j = [1.0 if i == j else 0.0 for i in range(n)]
        x_j = solve_system(L, U, P, e_j)
        for i in range(n):
            A_inv[i][j] = x_j[i]
    return A_inv

def calculate_determinant(U, swaps):
    # det(A)
    det = 1.0
    for i in range(len(U)):
        det *= U[i][i]
    if swaps % 2 != 0:
        det = -det
    return det


def main():
    if len(sys.argv) == 2:
        filename = sys.argv[1]
    elif len(sys.argv) >2:
        print (f"неверное кол-во аргументов")
        return
    else:
        filename = input("ввод имени:").strip()
        if not filename:
            filename = "input.txt"
    try:
        n, A, b = read_data(filename)
    except Exception as e:
        print(f"{e}")
        return
    
    try:
        L, U, P, swaps = lu_decompose(A)

        print_matrix("L", L)
        print_matrix("U", U)

        LU = matmul(L, U)
        print_matrix("L * U", LU)

        x = solve_system(L, U, P, b)
        print_vector("A*x = b", x)

        A_inv = calculate_inverse(L, U, P)
        print_matrix("A^-1", A_inv)

        det = calculate_determinant(U, swaps)
        print(f"det(A) = {det:.4f}")

        E_calculated = matmul(A, A_inv)
        E_true = create_identity_matrix(n)
        diff_E = matrix_difference(E_calculated, E_true)
        print_matrix("проверка A*A^-1", E_calculated)
        print(f"погрешность - {diff_E:.2e}")

        PA = matmul(P, A)
        diff_PA = matrix_difference(LU, PA)
        print_matrix("проверка P*A", PA)
        print(f"погрешность L*U == P*A - {diff_PA:.2e}")

    except ArithmeticError as ae:
        print(f"{ae}")

if __name__ == "__main__":
    main()