import sys
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

def read_input(path):
    if not os.path.exists(path):
        raise FileNotFoundError(f"'{path}' не найден")

    with open(path, "r", encoding="utf-8") as fh:
        lines = [ln.strip() for ln in fh if ln.strip()]

    if len(lines) < 3:
        raise ValueError("неверные аргументы")

    try:
        xstar = float(lines[0].replace(",", "."))
        xs = [float(t.replace(",", ".")) for t in lines[1].split()]
        fs = [float(t.replace(",", ".")) for t in lines[2].split()]
    except ValueError:
        raise ValueError("ошибка чтения чисел")

    if len(xs) != len(fs):
        raise ValueError("count x != count y")
    if any(xs[i] >= xs[i + 1] for i in range(len(xs) - 1)):
        raise ValueError("x_i должны возрастать")
    if len(xs) < 3:
        raise ValueError("неверное количество узлов")

    return xs, fs, xstar

def thomas(lower, diag, upper, rhs):
    m = len(diag)
    p, q = [0.0] * m, [0.0] * m
    p[0] = -upper[0] / diag[0]
    q[0] = rhs[0] / diag[0]
    for k in range(1, m):
        den = diag[k] + lower[k] * p[k - 1]
        p[k] = -upper[k] / den if k < m - 1 else 0.0
        q[k] = (rhs[k] - lower[k] * q[k - 1]) / den
    sol = [0.0] * m
    sol[-1] = q[-1]
    for k in range(m - 2, -1, -1):
        sol[k] = p[k] * sol[k + 1] + q[k]
    return sol

def build_spline(xs, fs):
    n = len(xs) - 1
    h = [xs[i] - xs[i - 1] for i in range(1, n + 1)]

    m = n - 1
    lower, diag, upper, rhs = [0.0] * m, [0.0] * m, [0.0] * m, [0.0] * m
    for k in range(m):
        i = k + 2
        hi_1, hi = h[i - 2], h[i - 1]
        lower[k] = hi_1
        diag[k] = 2.0 * (hi_1 + hi)
        upper[k] = hi
        rhs[k] = 3.0 * ((fs[i] - fs[i - 1]) / hi - (fs[i - 1] - fs[i - 2]) / hi_1)

    c_inner = thomas(lower, diag, upper, rhs)
    c = [0.0] + c_inner + [0.0]

    a, b, d = [], [], []
    for i in range(1, n + 1):
        hi = h[i - 1]
        ci, ci1 = c[i - 1], c[i]
        a.append(fs[i - 1])
        b.append((fs[i] - fs[i - 1]) / hi - hi * (ci1 + 2.0 * ci) / 3.0)
        d.append((ci1 - ci) / (3.0 * hi))

    return a, b, c[:n], d, h

def seg_value(xs, coef, i, x):
    a, b, c, d = coef
    t = x - xs[i - 1]
    return a[i - 1] + t * (b[i - 1] + t * (c[i - 1] + t * d[i - 1]))

def seg_d1(xs, coef, i, x):
    a, b, c, d = coef
    t = x - xs[i - 1]
    return b[i - 1] + t * (2.0 * c[i - 1] + 3.0 * d[i - 1] * t)

def seg_d2(xs, coef, i, x):
    a, b, c, d = coef
    t = x - xs[i - 1]
    return 2.0 * c[i - 1] + 6.0 * d[i - 1] * t

def find_interval(xs, x):
    for i in range(1, len(xs)):
        if xs[i - 1] <= x <= xs[i]:
            return i
    return None

COLORS = ["tab:red", "tab:blue", "tab:green", "tab:orange", "tab:purple"]

def plot(xs, fs, coef, xstar, i_star, s_star, fname):
    n = len(xs) - 1
    fig, ax = plt.subplots(figsize=(10, 6))
    
    for i in range(1, n + 1):
        g = [xs[i - 1] + (xs[i] - xs[i - 1]) * k / 100 for k in range(101)]
        ax.plot(g, [seg_value(xs, coef, i, t) for t in g], color=COLORS[(i - 1) % len(COLORS)],
                lw=2.5, label=f"$S_{i}(x)$ на [{xs[i-1]:g}; {xs[i]:g}]")
        
    for x in xs:
        ax.axvline(x, color="0.85", lw=0.8, ls="--")
        
    ax.plot(xs, fs, "ko", ms=8, zorder=5, label="узлы $(x_i, f_i)$")
    ax.plot([xstar], [s_star], "*", color="m", ms=16, zorder=6,
            label=f"$(x^*, S(x^*))$ = ({xstar:g}; {s_star:.5f})")
    ax.annotate(f"Интервал {i_star}", (xstar, s_star), textcoords="offset points",
                xytext=(20, -30), arrowprops=dict(arrowstyle="->"))
    
    ax.set_title("Кубический сплайн (Вариант 12)", fontsize=13)
    ax.set_xlabel("x")
    ax.set_ylabel("S(x)")
    ax.grid(True, linestyle="--", alpha=0.7)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(fname, dpi=130)
    plt.close(fig)

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "input.txt"
    try:
        xs, fs, xstar = read_input(path)
    except Exception as e:
        print(f"error: {e}")
        return

    n = len(xs) - 1
    a, b, c, d, h = build_spline(xs, fs)
    coef = (a, b, c, d)

    print("\n1 исходные узлы и значения функции:")
    print(f"{'i':>2} {'x_i':>10} {'f_i':>12}")
    for i in range(n + 1):
        print(f"{i:>2} {xs[i]:>10.4f} {fs[i]:>12.5f}")

    print("\n2 коэффициенты сплайна:")
    print(f"{'i':>2} | {'[x_{i-1}, x_i]':>14} | {'a_i':>10} | {'b_i':>10} | {'c_i':>10} | {'d_i':>10}")
    for i in range(1, n + 1):
        inv_str = f"[{xs[i-1]:g}, {xs[i]:g}]"
        print(f"{i:>2} | {inv_str:>14} | {a[i-1]:10.5f} | {b[i-1]:10.5f} | {c[i-1]:10.5f} | {d[i-1]:10.5f}")

    print("\n3 аналитический вид:")
    for i in range(1, n + 1):
        print(f"  S_{i}(x) = {a[i-1]:.5f} {b[i-1]:+.5f}*(x-{xs[i-1]:g}) "
              f"{c[i-1]:+.5f}*(x-{xs[i-1]:g})^2 {d[i-1]:+.5f}*(x-{xs[i-1]:g})^3,  "
              f"x in [{xs[i-1]:g}, {xs[i]:g}]")

    i_star = find_interval(xs, xstar)
    if i_star is None:
        print(f"\nx* = {xstar} вне отрезка интерполяции")
        return

    s_star = seg_value(xs, coef, i_star, xstar)
    print(f"\n4 в x* = {xstar}:")
    print(f"   интервал{i_star}: [{xs[i_star-1]:g}, {xs[i_star]:g}]")
    print(f"   S(x*) = {s_star:.6f}")

    worst = max(abs(fs[i] - (seg_value(xs, coef, 1, xs[0]) if i == 0 else seg_value(xs, coef, i, xs[i]))) for i in range(n + 1))
    print(f"точность прохождения S(x_i) = f_i: {worst:.1e}")

    max_j0 = max(abs(seg_value(xs, coef, i, xs[i]) - seg_value(xs, coef, i + 1, xs[i])) for i in range(1, n))
    max_j1 = max(abs(seg_d1(xs, coef, i, xs[i]) - seg_d1(xs, coef, i + 1, xs[i])) for i in range(1, n))
    max_j2 = max(abs(seg_d2(xs, coef, i, xs[i]) - seg_d2(xs, coef, i + 1, xs[i])) for i in range(1, n))
    print(f"гладкость на стыках участков S: {max_j0:.1e}, S': {max_j1:.1e}, S'': {max_j2:.1e}")

    s2_0 = abs(seg_d2(xs, coef, 1, xs[0]))
    s2_n = abs(seg_d2(xs, coef, n, xs[n]))
    print(f"нулевая кривизна на концах: S''(x_0) = {s2_0:.1e}, S''(x_n) = {s2_n:.1e}")

    plot(xs, fs, coef, xstar, i_star, s_star, "plot_spline.png")

if __name__ == "__main__":
    main()