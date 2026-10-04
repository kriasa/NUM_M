import sys
import math
import os

EPS_ZERO = 1e-12
MAX_ITERATIONS = 1000

def f(x):
    return 3 ** x - 5 * x ** 2 + 1

def f_prime(x):
    return (3 ** x) * math.log(3) - 10 * x

def f_second(x):
    return (3 ** x) * (math.log(3) ** 2) - 10


def find_optimal_lambda(a, b):
    mid = (a + b) / 2.0
    df_mid = f_prime(mid)

    lambd = 1.0 / df_mid
    return lambd

def phi(x, lambd):
    return x - lambd * f(x)

def phi_prime(x, lambd):
    return 1.0 - lambd * f_prime(x)

def estimate_q(a, b, lambd, samples=500):
    q = 0.0
    for i in range(samples + 1):
        x = a + (b - a) * i / samples
        q = max(q, abs(phi_prime(x, lambd)))
    return q

def read_input(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"'{path}' не найден.")

    with open(path, "r", encoding="utf-8") as fh:
        lines = [line.strip() for line in fh if line.strip()]

    if len(lines) < 2:
        raise ValueError("некорректный ввод")

    try:
        eps = float(lines[0])
        if eps <= 0:
            raise ValueError("eps > 0!")
            
        intervals = []
        for line in lines[1:]:
            parts = line.split()
            if len(parts) != 2:
                raise ValueError(f"'{line}' - [a,b]")
            a, b = float(parts[0]), float(parts[1])
            if a >= b:
                raise ValueError(f"неверный отрезок [{a}, {b}]")
            if a>0 and b>0:
                intervals.append((a, b))
            
        return eps, intervals
    except ValueError as e:
        raise ValueError(f"ERROR: {e}")

def newton_method(x0, eps):
    history = [(0, x0, f(x0), 0.0)]
    x_prev = x0

    for k in range(1, MAX_ITERATIONS):
        fx = f(x_prev)
        dfx = f_prime(x_prev)

        if abs(dfx) < EPS_ZERO:
            raise ValueError(f"f'(x)=0 на итерации {k}")

        x_next = x_prev - fx / dfx
        diff = abs(x_next - x_prev)
        history.append((k, x_next, f(x_next), diff))

        if diff < eps:
            return x_next, k, history

        x_prev = x_next

    raise ValueError(f"Метод Ньютона не сошёлся за {MAX_ITERATIONS} итераций")

def simple_iteration_method(a, b, eps):
    lambd = find_optimal_lambda(a, b)
    q = estimate_q(a, b, lambd)
    
    if q >= 1:
        raise ValueError(f"q > 1")

    x0 = (a + b) / 2.0
    history = [(0, x0, f(x0), 0.0)]
    x_prev = x0
    error_multiplier = q / (1.0 - q)

    for k in range(1, MAX_ITERATIONS):
        x_next = phi(x_prev, lambd)
        diff = abs(x_next - x_prev)
        history.append((k, x_next, f(x_next), diff))

        bound = error_multiplier * diff

        if bound <= eps:
            return x_next, k, history, q

        x_prev = x_next

    raise ValueError(f"Метод простой итерации не сошёлся за {MAX_ITERATIONS} итераций")

def print_table(history):
    print(f"{'k':>3} {'x(k)':>14} {'f(x(k))':>14} {'|x(k) - x(k-1)|':>18}")
    for k, x, fx, diff in history:
        diff_str = f"{diff:.2e}" if k > 0 else "---"
        print(f"{k:>3} {x:>14.8f} {fx:>14.8f} {diff_str:>18}")
    print()

def main():
    if len(sys.argv) > 1:
        path = sys.argv[1]
    else:
        path = input("файл: ").strip()
        if not path:
            path = "input.txt"

    try:
        eps, intervals = read_input(path)
    except Exception as e:
        print(f"ERROR: {e}")
        return

    print(f"3^x - 5x^2 + 1 = 0")
    print(f"eps = {eps}\n")

    for a, b in intervals:
        print(f"[{a}, {b}]")

        fa, fb = f(a), f(b)
        print(f"f(a) = {fa:.6f}, f(b) = {fb:.6f}")
        if fa * fb > 0:
            print(f"на [{a}, {b}] - f(a)*f(b) > 0\n")
            continue
        print("проверка - f(a)*f(b) < 0 ==> +")

        x0_newton = None
        for candidate in [a, b]:
            if f(candidate) * f_second(candidate) > 0:
                x0_newton = candidate
                break
        
        if x0_newton is None:
            print("ни один конец отрезка не удовлетворяет f(x0)*f''(x0) > 0\n")
            continue
        print(f"x0 = {x0_newton}\n")

        try:
            x_newton, iter_newton, hist_newton = newton_method(x0_newton, eps)
            print(f"Метод Ньютона (x0 = {x0_newton})")
            print_table(hist_newton)
            print(f"x* = {x_newton:.8f} (итераций: {iter_newton})")
        except Exception as e:
            print(f"Ошибка метода Ньютона: {e}")

        try:
            x_simple, iter_simple, hist_simple, q = simple_iteration_method(a, b, eps)
            print(f"\n МПИ (q = {q:.4f})")
            print_table(hist_simple)
            print(f"x* = {x_simple:.8f} (итераций: {iter_simple})")
        except Exception as e:
            print(f"Ошибка МПИ: {e}")

        print("\n")

if __name__ == "__main__":
    main()