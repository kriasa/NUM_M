import sys
import os

def read_input(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"'{path}' не найден.")

    with open(path, "r", encoding="utf-8") as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]

    if len(lines) < 3:
        raise ValueError("неверное кол-во аргументов")

    try:
        xstar = float(lines[0].replace(",", "."))
        xs = [float(t.replace(",", ".")) for t in lines[1].split()]
        ys = [float(t.replace(",", ".")) for t in lines[2].split()]
    except ValueError:
        raise ValueError("error format file")

    if len(xs) != len(ys) or len(xs) < 3:
        raise ValueError("count x = count y")
    if any(xs[i] >= xs[i + 1] for i in range(len(xs) - 1)):
        raise ValueError("x_i должны возрастать.")
    if not (xs[0] <= xstar <= xs[-1]):
        raise ValueError("x* лежит вне диапазона")
    return xs, ys, xstar


def choose_start(xs, xstar):
    M = len(xs) - 1
    eps = 1e-12
    for k in range(M + 1):
        if abs(xs[k] - xstar) < eps:
            return min(max(k - 1, 0), M - 2)
    for i in range(M):
        if xs[i] < xstar < xs[i + 1]:
            return min(i, M - 2)
    return 0

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "input.txt"
    try:
        xs, ys, xstar = read_input(path)
    except Exception as e:
        print(f"error: {e}")
        return

    print("1 Исходная таблица:")
    print(f"{'i':>2} {'x_i':>10} {'y_i':>12}")
    print("-" * 26)
    for i, (x, y) in enumerate(zip(xs, ys)):
        print(f"{i:>2} {x:>10.4f} {y:>12.5f}")

    print(f"\n2 x* = {xstar}")

    i = choose_start(xs, xstar)
    x0, x1, x2 = xs[i], xs[i + 1], xs[i + 2]
    y0, y1, y2 = ys[i], ys[i + 1], ys[i + 2]

    print("\n3 Узлы таблицы, использованные при вычислении (многочлен 2-й степени):")
    for k in (i, i + 1, i + 2):
        print(f"   i = {k}: x = {xs[k]:.4f}, y = {ys[k]:.5f}")

    f01 = (y1 - y0) / (x1 - x0)
    f12 = (y2 - y1) / (x2 - x1)
    f012 = (f12 - f01) / (x2 - x0)

    d1 = f01 + f012 * (2 * xstar - x0 - x1)
    d2 = 2 * f012

    print(f"f(x_i, x_i+1)= {f01:.6f}")
    print(f"f(x_i+1, x_i+2)= {f12:.6f}")
    print(f"f(x_i, x_i+1, x_i+2)= {f012:.6f}")
    print(f"4 y'(x*)  = {d1:.6f}")
    print(f"5 y''(x*) = {d2:.6f}")

    if abs(xstar - x1) < 1e-12:
        print(f"Левосторонняя производная: {f01:.6f}")
        print(f"Правосторонняя производная: {f12:.6f}")
        print(f"Их полусумма:{(f01 + f12) / 2:.6f}")


if __name__ == "__main__":
    main()