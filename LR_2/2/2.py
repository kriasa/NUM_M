import sys
import math
import os

MAX_ITERATIONS = 1000

def f1(x1, x2): 
    return x1 - math.cos(x2) - 3.0
def f2(x1, x2): 
    return x2 - math.sin(x1) - 3.0

def phi1(x1, x2): 
    return math.cos(x2) + 3.0
def phi2(x1, x2): 
    return math.sin(x1) + 3.0

def df1_dx1(x1, x2): 
    return 1.0
def df1_dx2(x1, x2): 
    return math.sin(x2)
def df2_dx1(x1, x2): 
    return -math.cos(x1)
def df2_dx2(x1, x2): 
    return 1.0

def dphi1_dx1(x1, x2): 
    return 0.0
def dphi1_dx2(x1, x2): 
    return -math.sin(x2)
def dphi2_dx1(x1, x2): 
    return math.cos(x1)
def dphi2_dx2(x1, x2): 
    return 0.0

def vector_norm_inf(v1, v2):
    return max(abs(v1), abs(v2))

def estimate_q(a1, b1, a2, b2, samples=50):
    q_max = 0.0
    step1 = (b1 - a1) / samples
    step2 = (b2 - a2) / samples
    
    for i in range(samples + 1):
        x1 = a1 + i * step1
        for j in range(samples + 1):
            x2 = a2 + j * step2
            
            row1_sum = abs(dphi1_dx1(x1, x2)) + abs(dphi1_dx2(x1, x2))
            row2_sum = abs(dphi2_dx1(x1, x2)) + abs(dphi2_dx2(x1, x2))
            current_q = max(row1_sum, row2_sum)
            
            if current_q > q_max:
                q_max = current_q
                
    return q_max

def read_input(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"'{path}' не найден.")

    with open(path, "r", encoding="utf-8") as fh:
        lines = [ln.strip() for ln in fh if ln.strip() != ""]

    if len(lines) < 2:
        raise ValueError("некорректный ввод")

    eps = float(lines[0])
    if eps <= 0:
        raise ValueError("eps > 0!")

    all_intervals = []
    for ln in lines[1:]:
        parts = ln.replace(",", ".").split()
        if len(parts) != 4:
            continue
        x1_lo, x1_hi, x2_lo, x2_hi = map(float, parts)
        all_intervals.append((x1_lo, x1_hi, x2_lo, x2_hi))

    positive_intervals = [
        iv for iv in all_intervals if all(v > 0 for v in iv)
    ]

    return eps, all_intervals, positive_intervals


def newton_system(x1_0, x2_0, eps):
    history = [(0, x1_0, x2_0, float('nan'))]
    x1_k, x2_k = x1_0, x2_0

    for k in range(1, MAX_ITERATIONS):
        F1 = f1(x1_k, x2_k)
        F2 = f2(x1_k, x2_k)
        
        J11, J12 = df1_dx1(x1_k, x2_k), df1_dx2(x1_k, x2_k)
        J21, J22 = df2_dx1(x1_k, x2_k), df2_dx2(x1_k, x2_k)
        
        det_J = J11 * J22 - J12 * J21
        if abs(det_J) < 1e-12:
            raise ValueError(f"detJ = {det_J:.2e} близок к 0")
            
        det_A1 = F1 * J22 - J12 * F2
        det_A2 = J11 * F2 - F1 * J21
        
        x1_next = x1_k - det_A1 / det_J
        x2_next = x2_k - det_A2 / det_J
        
        diff = vector_norm_inf(x1_next - x1_k, x2_next - x2_k)
        history.append((k, x1_next, x2_next, diff))
        
        if diff <= eps:
            return x1_next, x2_next, history
            
        x1_k, x2_k = x1_next, x2_next

    raise ValueError("Метод Ньютона не сошёлся")

def simple_iteration_system(a1, b1, a2, b2, x1_0, x2_0, eps):
    q = estimate_q(a1, b1, a2, b2)
    if q >= 1:
        raise ValueError(f"q = {q:.4f} >= 1.")
        
    history = [(0, x1_0, x2_0, float('nan'))]
    x1_k, x2_k = x1_0, x2_0

    error_multiplier = q / (1.0 - q)

    for k in range(1, MAX_ITERATIONS):
        x1_next = phi1(x1_k, x2_k)
        x2_next = phi2(x1_k, x2_k)
        
        diff = vector_norm_inf(x1_next - x1_k, x2_next - x2_k)
        history.append((k, x1_next, x2_next, diff))
        
        estimated_error = error_multiplier * diff
        if estimated_error <= eps:
            return x1_next, x2_next, history, q
            
        x1_k, x2_k = x1_next, x2_next

    raise ValueError("Метод простых итераций не сошёлся.")

def print_table(history, title):
    print(f"\n{title}")
    print(f"{'k':>3}  {'x1^(k)':>12}  {'x2^(k)':>12}  {'||x(k) - x(k-1)||':>20}")
    for k, x1, x2, diff in history:
        diff_str = f"{diff:.4e}" if k > 0 else "-"
        print(f"{k:>3}  {x1:>12.6f}  {x2:>12.6f}  {diff_str:>20}")

def main():
    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        path = input("файл: ").strip()
        if not path:
            path = "input.txt"

    try:
        eps, all_intervals, positive_intervals = read_input(path)
    except Exception as e:
        print(f"ERROR: {e}")
        return

    print(" f1: x1 - cos(x2) - 3 = 0")
    print(" f2: x2 - sin(x1) - 3 = 0")
    print(f"eps = {eps}")

    if not positive_intervals:
        print("Подходящих промежутков нет")
        return

    for a1, b1, a2, b2 in positive_intervals:
        print(f"x1 - [{a1}, {b1}], x2 - [{a2}, {b2}]")

        x1_0 = (a1 + b1) / 2.0
        x2_0 = (a2 + b2) / 2.0
        print(f"x0 = ({x1_0:.4f}, {x2_0:.4f})")

        J11, J12 = df1_dx1(x1_0, x2_0), df1_dx2(x1_0, x2_0)
        J21, J22 = df2_dx1(x1_0, x2_0), df2_dx2(x1_0, x2_0)
        detJ0 = J11 * J22 - J12 * J21
        print(f"det J(x0) = {detJ0:.4f} -> невырождена.\n")

        try:
            res_n1, res_n2, hist_n = newton_system(x1_0, x2_0, eps)
            print_table(hist_n, "Метод Ньютона:")
            print(f"x1* = {res_n1:.6f}, x2* = {res_n2:.6f}")
        except Exception as e:
            print(f"Ошибка метода Ньютона: {e}")

        try:
            res_s1, res_s2, hist_s, q = simple_iteration_system(a1, b1, a2, b2, x1_0, x2_0, eps)
            print_table(hist_s, f"МПИ (q = {q:.4f}):")
            print(f"x1* = {res_s1:.6f}, x2* = {res_s2:.6f}")
        except Exception as e:
            print(f"Ошибка МПИ: {e}")

if __name__ == "__main__":
    main()