import os
import sys
import math

EPS_ZERO = 1e-12
MAX_ITERATIONS = 10000

def read_input(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"'{filename}' не найден")
    
    with open(filename, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f if line.strip()]
        
    try:
        n = int(lines[0])
        if len(lines) != n + 2:
            raise ValueError(f"кол-во строк неверно")
        
        A = []
        for i in range(1, n + 1):
            row = list(map(float, lines[i].split()))
            if len(row) != n:
                raise ValueError(f"некорректна по размеру строка - {i}")
            A.append(row)
            
        eps = float(lines[n + 1])
        if eps <= 0:
            raise ValueError("некорректное {eps}")

        for i in range(n):
            for j in range(i + 1, n):
                if abs(A[i][j] - A[j][i]) > EPS_ZERO:
                    raise ValueError(f"Матрица не симметрична")
            
        return A, eps, n

    except Exception as e:
        raise ValueError(f"ERROR: {e}")

def create_identity(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

def create_zero(n):
    return [[0.0 for _ in range(n)] for _ in range(n)]

def transpose(A):
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]

def matmul(A, B):
    n = len(A)
    C = create_zero(n)
    for i in range(n):
        for j in range(n):
            for k in range(n):
                C[i][j] += A[i][k] * B[k][j]
    return C

def matrix_difference_norm(A, B):
    n = len(A)
    max_row_sum = 0.0
    for i in range(n):
        row_sum = sum(abs(A[i][j] - B[i][j]) for j in range(n))
        if row_sum > max_row_sum:
            max_row_sum = row_sum
    return max_row_sum

def jacobi_eigenvalue_algorithm(A_orig, eps, n):
    A = [row[:] for row in A_orig]
    V = create_identity(n) 
    iterations = 0
    while True:
        sum_sq = 0.0
        for i in range(n):
            for j in range(i + 1, n):
                sum_sq += A[i][j] ** 2
        t = math.sqrt(sum_sq)
        
        if t < eps:
            break
        if iterations >= MAX_ITERATIONS:
            raise ValueError(f"метод не сошёлся за {MAX_ITERATIONS} итераций")

        max_val = -1.0
        p, q = 0, 1
        for i in range(n):
            for j in range(i + 1, n):
                if abs(A[i][j]) > max_val:
                    max_val = abs(A[i][j])
                    p, q = i, j

        if max_val < EPS_ZERO:
            break

        if abs(A[p][p] - A[q][q]) < EPS_ZERO:
            phi = math.pi / 4.0
        else:
            phi = 0.5 * math.atan((2.0 * A[p][q]) / (A[p][p] - A[q][q]))
            
        c = math.cos(phi)
        s = math.sin(phi)

        U = create_identity(n)
        U[p][p] = c
        U[q][q] = c
        U[p][q] = -s
        U[q][p] = s

        U_T = transpose(U)
        A = matmul(matmul(U_T, A), U)

        V = matmul(V, U)
        iterations += 1

    eigenvalues = [A[i][i] for i in range(n)]
    
    return eigenvalues, V, iterations

def print_matrix(title, matrix):
    print(f"\n{title}:")
    for row in matrix:
        formatted = [f"{0.0 if abs(val) < EPS_ZERO else val:10.4f}" for val in row]
        print(" ".join(formatted))

def print_vector(title, vector):
    print(f"\n{title}:")
    formatted = [f"{0.0 if abs(val) < EPS_ZERO else val:10.4f}" for val in vector]
    print(" ".join(formatted))

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
        A_orig, eps, n = read_input(filename)

    except Exception as e:
        print(f"ERROR:\n{e}")
        return

    print(f"\n n = {n}")
    print(f" eps = {eps}")

    try:
        eigenvalues, V, iters = jacobi_eigenvalue_algorithm(A_orig, eps, n)
    except ValueError as e:
        print(f"\nERROR:\n{e}")
        return
        
    print(f"k итераций метода: {iters}")
    print_vector("1) СЗ", eigenvalues)
    print("\n2) СВ")
    for j in range(n):
        h_j = [V[i][j] for i in range(n)]
        clean_hj = [f"{0.0 if abs(x) < EPS_ZERO else x:.4f}" for x in h_j]
        print(f"   h_{j+1} = ( {', '.join(clean_hj)} )^T")
        
    print_matrix("Матрица СВ", V)

    Lambda = create_identity(n)
    for i in range(n):
        Lambda[i][i] = eigenvalues[i]

    AV = matmul(A_orig, V)
    VL = matmul(V, Lambda)
    
    print("\n3) проверка (A · V = V · Λ)")
    print_matrix("Матрица A · V", AV)
    print_matrix("Матрица V · Λ", VL)
    
    allowable_error = eps * n 
    diff = matrix_difference_norm(AV, VL)
    print(f"\nмаксимальное расхождение между матрицами- {diff:.2e}")
    if diff <= allowable_error:
        print("=> Проверка +")
    else:
        print("=> Проверка -")

if __name__ == "__main__":
    main()