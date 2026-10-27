import sys
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_input(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"'{path}' не найден.")

    with open(path, "r", encoding="utf-8") as fh:
        lines = [line.strip() for line in fh if line.strip()]

    if len(lines) < 2:
        raise ValueError("неверный ввод в файле")

    xs = [float(t.replace(",", ".")) for t in lines[0].split()]
    ys = [float(t.replace(",", ".")) for t in lines[1].split()]

    if len(xs) != len(ys):
        raise ValueError("count x = count y")
    return xs, ys


def gauss_solve(A, b):
    n = len(b)
    M = [A[i][:] + [b[i]] for i in range(n)]

    for k in range(n):
        pivot_row = k
        max_val = abs(M[k][k])
        for i in range(k + 1, n):
            if abs(M[i][k]) > max_val:
                max_val = abs(M[i][k])
                pivot_row = i

        if max_val < 1e-14: raise ValueError("det = 0")

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

def fit_poly(xs, ys, deg):
    A = [[sum(x ** (k + i) for x in xs) for i in range(5)] for k in range(5)]
    b = [sum(y * (x ** k) for x, y in zip(xs, ys)) for k in range(5)]
    return gauss_solve(A, b), A, b

def poly_val(a, x):
    r = 0.0
    for i in range(len(a)):
        r += a[i] * (x ** i)
    return r

def poly_str(a, name):
    s = f"{name}(x) = {a[0]:.6f}"
    for k in range(1, len(a)):
        sign = "+" if a[k] >= 0 else "-"
        s += f" {sign} {abs(a[k]):.6f}" + ("*x" if k == 1 else f"*x^{k}")
    return s

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "input.txt"
    try:
        xs, ys = read_input(path)
    except Exception as e:
        print(f"error: {e}")
        return

    N1 = len(xs)

    print(f"{'j':>2} {'x_j':>10} {'y_j':>10}")
    for j in range(N1):
        print(f"{j:>2} {xs[j]:>10.5f} {ys[j]:>10.5f}")

    coefs, mse, sse = {}, {}, {}
    for deg in (1, 2, 3):
        a, A, b = fit_poly(xs, ys, deg)
        coefs[deg] = a
        sse[deg] = sum((poly_val(a, x) - y) ** 2 for x, y in zip(xs, ys))
        mse[deg] = sse[deg] / N1
        
        print(f"\nприближающий многочлен степени {deg}")
        print("МНК (расширенная матрица):")

        print("" + ", ".join(f"a{i} = {c:.6f}" for i, c in enumerate(a)))
        print("" + poly_str(a, f"P{deg}"))
        print(f"Φ = {sse[deg]:.6f}")
        print(f"MSE{deg} = {mse[deg]:.6f}")
        
    best = min(mse, key=mse.get)
    print(f"\nлучше данные приближает многочлен P{best}(x) (min MSE = {mse[best]:.6f})")

    lo, hi = min(xs), max(xs)
    pad = 0.05 * (hi - lo)
    g = [lo - pad + (hi - lo + 2 * pad) * i / 400 for i in range(401)]
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(xs, ys, "ko", ms=9, zorder=5, label="исходные данные")
    
    styles = {1: ("tab:red", "-"), 2: ("tab:blue", "--"), 3: ("tab:green", "-.")}
    for deg in (1, 2, 3):
        col, ls = styles[deg]
        ax.plot(g, [poly_val(coefs[deg], t) for t in g], color=col, ls=ls, lw=2.2,
                label=f"$P_{deg}(x)$, MSE = {mse[deg]:.4f}")
                
    ax.set_title("МНК", fontsize=13)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.grid(True, linestyle="--", alpha=0.7)
    ax.legend(fontsize=11)
    fig.tight_layout()
    fig.savefig("plot_lsq.png", dpi=130)
    plt.close(fig)

if __name__ == "__main__":
    main()