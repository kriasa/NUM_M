import os
import sys

EPS = 1e-12

def read_diagonals(filename):
    if not os.path.exists(filename):
        raise FileNotFoundError(f"'{filename}' не найден")
    
    with open(filename, 'r', encoding='utf-8') as f:
        lines = [line.strip() for line in f if line.strip()]
        
    try:
        n = int(lines[0])

        if len(lines) != n:
            raise ValueError(f"кол-во строк неверно")

        a = list(map(float, lines[1].split()))
        b = list(map(float, lines[2].split()))
        c = list(map(float, lines[3].split()))
        d = list(map(float, lines[4].split()))
        
        if not (len(a) == len(b) == len(c) == len(d) == n):
            raise ValueError("неправильная длина векторов")

        if abs(a[0]) > EPS:
            raise ValueError("A1 = 0!")
        if abs(c[n-1]) > EPS:
            raise ValueError("Cn = 0!")
        return n, a, b, c, d
    except Exception as e:
        raise ValueError(f"ERROR: {e}")

def check_stability(a, b, c):
    strict_flag = False

    for i in range(len(b)):
        sum_ac = abs(a[i]) + abs(c[i])
        if abs(b[i]) < sum_ac - EPS:
            return False, f"устойчивость нарушена - строка {i+1}: |{b[i]}| < |{a[i]}| + |{c[i]}|"
        if abs(b[i]) > sum_ac + EPS:
            strict_flag = True
            
    if strict_flag:
        return True, "метод устойчив"
    return False, "устойчивость не гарантируется"

def thomas_algorithm(n, a, b, c, d):

    P = [0.0] * n
    Q = [0.0] * n
    x = [0.0] * n

    denom = b[0]
    if abs(denom) < EPS:
        raise ZeroDivisionError("B1 = 0")
        
    P[0] = -c[0] / denom
    Q[0] = d[0] / denom
    
    for i in range(1, n):
        denom = b[i] + a[i] * P[i - 1]

        if abs(denom) < EPS:
            raise ZeroDivisionError(f"деление на ноль на шаге {i+1}.")
            
        P[i] = -c[i] / denom
        Q[i] = (d[i] - a[i] * Q[i - 1]) / denom

    P[n - 1] = 0.0
    x[n - 1] = Q[n - 1]
    
    for i in range(n - 2, -1, -1):
        x[i] = P[i] * x[i + 1] + Q[i]
        
    return P, Q, x

def print_coefficients_table(P, Q, n):
    print("\nпрогоночные коэффициенты")
    print(f"{'i':>3} | {'P_i':>12} | {'Q_i':>12}\n")

    for i in range(n):
        clean_P = 0.0 if abs(P[i]) < EPS else P[i]
        clean_Q = 0.0 if abs(Q[i]) < EPS else Q[i]
        print(f"{i+1:>3} | {clean_P:>12.6f} | {clean_Q:>12.6f}")

def print_solution(x, n):
    for i in range(n):
        clean_x = 0.0 if abs(x[i]) < EPS else x[i]
        print(f"x[{i+1}] = {clean_x:.6f}")

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
        n, a, b, c, d = read_diagonals(filename)
    except Exception as e:
        print(f"ERROR:{e}")
        return
    
    is_stable, msg = check_stability(a, b, c)
    print(f"|Bi|>=|Ai|+|Ci|: {msg}")

    try:
        P, Q, x = thomas_algorithm(n, a, b, c, d)
    except ZeroDivisionError as e:
        print(f"\nERROR:{e}")
        return

    print_coefficients_table(P, Q, n)
    print_solution(x, n)

if __name__ == "__main__":
    main()