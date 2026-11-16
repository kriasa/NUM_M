import math
import sys
import os

def f(x):
    return x / (x ** 3 + 8)

def F_exact(x):
    sqrt3 = math.sqrt(3)
    return (-math.log(x + 2) / 6.0 
            + math.log(x * x - 2 * x + 4) / 12.0 
            + (sqrt3 / 6.0) * math.atan((x - 1) / sqrt3))

HAS_EXACT = True

def read_input(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"'{path}' не найден")

    with open(path, "r", encoding="utf-8") as fh:
        tokens = fh.read().split()

    if len(tokens) < 4:
        raise ValueError("некорректное кол-во параметров")

    a = float(tokens[0])
    b = float(tokens[1])
    h1 = float(tokens[2])
    h2 = float(tokens[3])

    if a >= b:
        raise ValueError("a < b")
    if h1 <= 0 or h2 <= 0:
        raise ValueError("h > 0")
    if h1 <= h2:
        raise ValueError("h1 > h2")

    return a, b, h1, h2

def n_intervals(a, b, h):
    n = round((b - a) / h)
    if abs(n * h - (b - a)) > 1e-9:
        raise ValueError(f"{h} не делит[{a}, {b}] нацело")
    return n

def left_rect(a, b, h):
    n = n_intervals(a, b, h)
    return h * sum(f(a + i * h) for i in range(n))

def right_rect(a, b, h):
    n = n_intervals(a, b, h)
    return h * sum(f(a + i * h) for i in range(1, n + 1))

def mid_rect(a, b, h):
    n = n_intervals(a, b, h)
    return h * sum(f(a + (i + 0.5) * h) for i in range(n))

def trapezoid(a, b, h):
    n = n_intervals(a, b, h)
    return h * (f(a) / 2.0 + sum(f(a + i * h) for i in range(1, n)) + f(b) / 2.0)

def simpson(a, b, h):
    n = n_intervals(a, b, h)
    if n % 2 != 0:
        raise ValueError("для метода Симпсона число интервалов должно быть чётным")
    s = f(a) + f(b)
    for i in range(1, n):
        coef = 4.0 if i % 2 != 0 else 2.0
        s += coef * f(a + i * h)
    return (h / 3.0) * s

METHODS = [
    ("Левых прямоугольников", left_rect, 1),
    ("Правых прямоугольников", right_rect, 1),
    ("Средних прямоугольников", mid_rect, 2),
    ("Трапеций", trapezoid, 2),
    ("Симпсона", simpson, 4),
]

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "input.txt"
    try:
        a, b, h1, h2 = read_input(path)
    except Exception as e:
        print(f"ERROR: {e}")
        return

    k = h1 / h2
    print(f"y = x / (x^3 + 8) на [{a:g}, {b:g}]")
    print(f"h1 = {h1:g} (N1 = {n_intervals(a, b, h1)}), h2 = {h2:g} (N2 = {n_intervals(a, b, h2)}),  k = h1/h2 = {k:g}")

    rows = []
    for name, method, p in METHODS:
        I1 = method(a, b, h1)
        I2 = method(a, b, h2)
        corr = (I2 - I1) / (k ** p - 1.0)
        I_refined = I2 + corr
        est_err = abs(corr)
        rows.append((name, I1, I2, p, I_refined, est_err))

    print("Таблица численного интегрирования и уточнения по Рунге-Ромбергу:")
    head = f"{'Метод':<25} | {'I(h1)':>10} | {'I(h2)':>10} | {'p':>3} | {'Уточнённое':>12} | {'Оценка погр.(h2)':>18}"
    print(head)
    for name, I1, I2, p, Ir, est in rows:
        print(f"{name:<25} | {I1:10.6f} | {I2:10.6f} | {p:>3} | {Ir:12.6f} | {est:18.2e}")

    if HAS_EXACT:
        exact = F_exact(b) - F_exact(a)
        print(f"Точное значение интеграла (по формуле Ньютона-Лейбница): {exact:.10f}")
        
        head2 = f"{'Метод':<25} | {'|погр.| h1':>12} | {'|погр.| h2':>12} | {'|погр.| уточн.':>15}"
        print(head2)
        
        err2 = {}
        for name, I1, I2, p, Ir, est in rows:
            err2[name] = abs(exact - I2)
            print(f"{name:<25} | {abs(exact - I1):12.2e} | {abs(exact - I2):12.2e} | {abs(exact - Ir):15.2e}")

        best_name = None
        min_err = float('inf')
        for name, err in err2.items():
            if err < min_err:
                min_err = err
                best_name = name
                
        print(f"\nНаиболее точный базовый метод (с шагом h2): {best_name} ({min_err:.2e})")

        best_r_name = None
        min_r_err = float('inf')
        for r in rows:
            err_r = abs(exact - r[4])
            if err_r < min_r_err:
                min_r_err = err_r
                best_r_name = r[0]
                
        print(f"После уточнения по Рунге-Ромбергу лучший результат: {best_r_name} ({min_r_err:.2e})")

if __name__ == "__main__":
    main()