import os
import sys
import cmath

EPS_ZERO = 1e-12
MAX_ITERATIONS = 5000

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
            raise ValueError(f"некорректное {eps}")
            
        return A, eps, n
    except Exception as e:
        raise ValueError(f"ERROR: {e}")

def identity(n):
    return [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]

def copy_matrix(M):
    return [row[:] for row in M]

def create_zero_matrix(rows, cols):
    return [[0.0 for _ in range(cols)] for _ in range(rows)]

def matmul(A, B):
    rows_A = len(A)
    cols_A = len(A[0])
    cols_B = len(B[0])
    
    C = create_zero_matrix(rows_A, cols_B)
    
    for i in range(rows_A):
        for j in range(cols_B):
            for k in range(cols_A):
                C[i][j] += A[i][k] * B[k][j]
    return C

def max_abs_diff_matrix(A, B):
    n = len(A)
    max_row_sum = 0.0
    for i in range(n):
        row_sum = sum(abs(A[i][j] - B[i][j]) for j in range(n))
        if row_sum > max_row_sum:
            max_row_sum = row_sum
    return max_row_sum

def householder_qr(A, n):
    R = copy_matrix(A)
    Q = identity(n)

    for k in range(n - 1):
        norm = sum(R[i][k] ** 2 for i in range(k, n)) ** 0.5
        if norm < EPS_ZERO:
            continue

        sign = 1.0 if R[k][k] >= 0 else -1.0
        v = [0.0] * n
        v[k] = R[k][k] + sign * norm
        for i in range(k + 1, n):
            v[i] = R[i][k]

        vtv = sum(x * x for x in v)
        if vtv < EPS_ZERO:
            continue

        H = identity(n)
        for i in range(n):
            if abs(v[i]) < EPS_ZERO:
                continue
            for j in range(n):
                H[i][j] -= 2.0 * v[i] * v[j] / vtv

        R = matmul(H, R)
        Q = matmul(Q, H)

    return Q, R

def extract_eigen_blocks(A, n, eps):
    eigenvalues = []
    is_block2 = []
    i = 0
    while i < n:
        if i == n - 1:
            eigenvalues.append(complex(A[i][i]))
            is_block2.append(False)
            i += 1
            continue

        subdiag = sum(abs(A[j][i]) for j in range(i + 1, n))
        if subdiag < eps:
            eigenvalues.append(complex(A[i][i]))
            is_block2.append(False)
            i += 1
        else:
            a, b = A[i][i], A[i][i + 1]
            c, d = A[i + 1][i], A[i + 1][i + 1]
            tr = a + d
            det = a * d - b * c
            disc = tr * tr - 4 * det
            sq = cmath.sqrt(disc)
            
            l1 = (tr + sq) / 2
            l2 = (tr - sq) / 2
            eigenvalues.append(l1)
            eigenvalues.append(l2)
            is_block2.append(True)
            is_block2.append(True)
            i += 2

    return eigenvalues, is_block2
    
def qr_algorithm(A, n, eps):
    A_k = copy_matrix(A)
    prev_eigenvalues = None
    prev_is_block2 = None

    if n == 1:
        return [complex(A_k[0][0])], A_k, identity(1), A_k

    for iteration in range(1, MAX_ITERATIONS):
        Q, R = householder_qr(A_k, n)
        A_k = matmul(R, Q)

        eigenvalues, is_block2 = extract_eigen_blocks(A_k, n, eps)

        converged = True

        for m in range(n - 1):
            if not is_block2[m]:
                sub_sum_sq = sum(A_k[l][m] ** 2 for l in range(m + 1, n))
                if math.sqrt(sub_sum_sq) > eps:
                    converged = False
                    break

        if converged:
            if prev_eigenvalues is not None:
                if is_block2 != prev_is_block2:
                    converged = False
                else:
                    for k in range(n):
                        if is_block2[k]:
                            if abs(eigenvalues[k] - prev_eigenvalues[k]) > eps:
                                converged = False
                                break
            else:
                if any(is_block2):
                    converged = False

        prev_eigenvalues = eigenvalues
        prev_is_block2 = is_block2

        if converged:
            Q_final, R_final = householder_qr(A_k, n)
            return eigenvalues, A_k, Q_final, R_final

    raise ValueError(f"QR-алгоритм не сошёлся за {MAX_ITERATIONS} итераций.")

def format_complex(z):
    if abs(z.imag) < 1e-9:
        return f"{z.real:10.4f}"
    sign = "+" if z.imag >= 0 else "-"
    return f"{z.real:10.4f} {sign} {abs(z.imag):.4f}i"

def print_matrix(title, matrix):
    print(f"\n{title}:")
    for row in matrix:
        formatted = [f"{0.0 if abs(val) < EPS_ZERO else val:10.4f}" for val in row]
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
        A, eps, n = read_input(filename)
    except Exception as e:
        print(f"\nERROR:\n{e}")
        return

    print(f"n = {n}")
    print(f"eps = {eps}")
    print_matrix("A начальная", A)

    try:
        Q, R = householder_qr(A, n)
        print_matrix("Q", Q)
        print_matrix("R", R)

        QR = matmul(Q, R)
        print_matrix("Q х R", QR)
        diff_qr = max_abs_diff_matrix(QR, A)
        print(f"проверка - отклонение Q х R от A: {diff_qr:.2e}")

        eigenvalues, A_final, Q_final, R_final = qr_algorithm(A, n, eps)
        print_matrix("Матрица A (последняя)", A_final)
        print_matrix("Матрица Q (последняя)", Q_final)
        print_matrix("Матрица R (последняя)", R_final)
        QR_f = matmul(Q_final, R_final)
        print_matrix("Q х R", QR_f)
        diff_qr_f = max_abs_diff_matrix(QR_f, A_final)
        print(f"проверка - отклонение Q х R от A: {diff_qr_f:.2e}")
        
        
        print("\nСЗ")
        for i, val in enumerate(eigenvalues):
            print(f"   h_{i+1} = {format_complex(val)}")

    except Exception as e:
        print(f"\nERROR:\n{e}")
        return

    has_complex = any(abs(z.imag) > 1e-6 for z in eigenvalues)
    if has_complex:
        print("Проверка комплексных корней")
        try:
            import numpy as np
            np_eigs = np.linalg.eigvals(np.array(A))
            print("СЗ через NumPy")
            for idx, z in enumerate(np_eigs):
                sign = "+" if z.imag >= 0 else "-"
                print(f"  сз_{idx+1} = {z.real:.4f} {sign} {abs(z.imag):.4f}i")
        except ImportError:
            print("ERROR")

if __name__ == "__main__":
    main()