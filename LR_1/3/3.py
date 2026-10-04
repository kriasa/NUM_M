import os
import sys

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
            if len(row) != n + 1:
                raise ValueError(f"некорректна по размеру строка - {i} ")
            A.append(row[:n])
            if i == 1:
                b = [row[n]]
            else:
                b.append(row[n])
                
        eps = float(lines[n + 1])
        if eps <= 0:
            raise ValueError(f"некорректное {eps}")
            
        return A, b, eps, n

    except Exception as e:
        raise ValueError(f"ERROR: {e}")

def build_iteration_form(A, b, n):
    A_work = [row[:] for row in A]
    b_work = list(b)

    for i in range(n):
        if abs(A_work[i][i]) < EPS_ZERO:
            swap_row = -1
            max_val = 0.0
            for r in range(i + 1, n):
                if abs(A_work[r][i]) > max_val:
                    max_val = abs(A_work[r][i])
                    swap_row = r
            
            if swap_row != -1 and max_val >= EPS_ZERO:
                A_work[i], A_work[swap_row] = A_work[swap_row], A_work[i]
                b_work[i], b_work[swap_row] = b_work[swap_row], b_work[i]
            else:
                raise ValueError(
                    f"система вырождена, метод неприменим."
                )

    alpha = [[0.0] * n for _ in range(n)]
    beta = [0.0] * n

    for i in range(n):
        diag = A_work[i][i]
        beta[i] = b_work[i] / diag
        for j in range(n):
            if j != i:
                alpha[i][j] = -A_work[i][j] / diag

    return alpha, beta

def compute_norms(alpha, n):
    norm_alpha = 0.0
    norm_c = 0.0
    
    for i in range(n):
        sum_alpha = sum(abs(alpha[i][j]) for j in range(n))
        sum_c = sum(abs(alpha[i][j]) for j in range(i, n))
        
        norm_alpha = max(norm_alpha, sum_alpha)
        norm_c = max(norm_c, sum_c)
            
    return norm_alpha, norm_c

def vector_inf_diff(x_new, x_old, n):
    return max(abs(x_new[i] - x_old[i]) for i in range(n))

def simple_iteration_method(alpha, beta, n, eps, norm_alpha):
    x_old = list(beta)
    iterations = 0
    error_factor = norm_alpha / (1.0 - norm_alpha) if norm_alpha < 1.0 else 1.0
    while True:
        x_new = [0.0] * n
        for i in range(n):
            s = beta[i]
            for j in range(n):
                s += alpha[i][j] * x_old[j]
            x_new[i] = s

        iterations += 1
        diff = vector_inf_diff(x_new, x_old, n)
        estimated_error = diff * error_factor if norm_alpha < 1.0 else diff
        x_old = list(x_new)

        if estimated_error < eps:
            break
        if iterations >= MAX_ITERATIONS:
            raise ValueError(f"МПИ не сходится")

    return x_old, iterations

def seidel_method(alpha, beta, n, eps, norm_alpha, norm_c):
    x = list(beta)
    iterations = 0
    error_factor = norm_c / (1.0 - norm_alpha) if norm_alpha < 1.0 else 1.0
    while True:
        x_old = list(x)

        for i in range(n):
            s = beta[i]
            for j in range(i):
                s += alpha[i][j] * x[j]
            for j in range(i + 1, n):
                s += alpha[i][j] * x_old[j]
                
            x[i] = s

        iterations += 1
        diff = vector_inf_diff(x, x_old, n)
        estimated_error = diff * error_factor if norm_alpha < 1.0 else diff

        if estimated_error < eps:
            break
        if iterations >= MAX_ITERATIONS:
            raise ValueError(f"Метод Зейделя не сошёлся")

    return x, iterations

def gauss_reference_solution(A, b, n):
    M = [A[i][:] + [b[i]] for i in range(n)]

    for k in range(n):
        pivot_row = k
        max_val = abs(M[k][k])
        for i in range(k + 1, n):
            if abs(M[i][k]) > max_val:
                max_val = abs(M[i][k])
                pivot_row = i

        if max_val < EPS_ZERO:
            raise ValueError("det = 0 вырождена")

        M[k], M[pivot_row] = M[pivot_row], M[k]
        for i in range(k + 1, n):
            factor = M[i][k] / M[k][k]
            for j in range(k, n + 1):
                M[i][j] -= factor * M[k][j]

    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = M[i][n]
        for j in range(i + 1, n):
            s -= M[i][j] * x[j]
        x[i] = s / M[i][i]

    return x

def print_vector(v, title):
    print(title)
    for i, val in enumerate(v):
        print(f"  x_{i+1} = {val:10.6f}")

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
        A, b, eps, n = read_input(filename)
    except Exception as e:
        print(f"ошибка при чтении входных данных {e}")
        return

    print(f"\nn = {n}")
    print(f"eps = {eps}\n")

    try:
        alpha, beta = build_iteration_form(A, b, n)
        norm_alpha, norm_c = compute_norms(alpha, n)

        print(f"||alpha|| = {norm_alpha:.4f}")

        if norm_alpha < 1.0:
            print("=> сходится\n")
        else:
            print("=> не факт что сходится\n")

        x_simple, iter_simple = simple_iteration_method(alpha, beta, n, eps, norm_alpha)
        x_seidel, iter_seidel = seidel_method(alpha, beta, n, eps, norm_alpha, norm_c)
        x_exact = gauss_reference_solution(A, b, n)
        
    except Exception as e:
        print(f"ERROR: {e}")
        return

    print(f"1) k итераций:")
    print(f"   - МПИ (Якоби) {iter_simple}")
    print(f"   - Метод Зейделя {iter_seidel}\n")
    
    print_vector(x_exact, "2) эталонное решение системы (метод Гаусса):")
    print_vector(x_simple, "МПИ:")
    print_vector(x_seidel, "Методом Зейделя:")

    print("3) Сравнение по скорости сходимости:")
    if iter_seidel < iter_simple:
        print(f"   Метод Зейделя быстрее сходится на {iter_simple - iter_seidel} итераций")
    elif iter_simple < iter_seidel:
        print(f"   МПИ быстрее сходится на {iter_seidel - iter_simple} итераций")
    else:
        print(f"   Сошлись за одинаковое число итераций ({iter_simple})")

if __name__ == "__main__":
    main()