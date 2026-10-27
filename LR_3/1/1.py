import math
import sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from fractions import Fraction

def f(x):
    return math.sin(x) + x

def f_deriv(k, x):
    if k == 1:
        return math.cos(x) + 1
    if k == 2:
        return -math.sin(x)
    if k == 3:
        return -math.cos(x)
    if k == 4:
        return math.sin(x)
    if k > 4:
        return f_deriv(k - 4, x)
    else:
        return 0.0

def parse_expr(s):
    return float(eval(s.strip(), {"__builtins__": {}}, {"pi": math.pi}))

def read_input(path):
    with open(path, "r", encoding="utf-8") as fh: lines = [line.strip() for line in fh if line.strip()]

    if len(lines) < 2: raise ValueError("некорректный формат файла")

    xstar = parse_expr(lines[0])
    sets = {}
    for idx, line in enumerate(lines[1:], start=1):
        name = "а" if idx == 1 else ("б" if idx == 2 else str(idx))
        nodes = [parse_expr(t) for t in line.split()]
        sets[name] = nodes

    return sets, xstar

def lagrange_weights(xs, ys):
    w = []
    for i in range(len(xs)):
        d = 1.0
        for j in range(len(xs)):
            if j != i:
                d *= xs[i] - xs[j]
        w.append(ys[i] / d)
    return w

def lagrange_value(xs, ys, x):
    total = 0.0
    for i in range(len(xs)):
        li = 1.0
        for j in range(len(xs)):
            if j != i:
                li *= (x - xs[j]) / (xs[i] - xs[j])
        total += ys[i] * li
    return total

def lagrange_str(xs, ys):
    w = lagrange_weights(xs, ys)
    out = []
    for i in range(len(xs)):
        fac = "".join(f"(x - {xs[j]:.6f})" for j in range(len(xs)) if j != i)
        out.append(f"{w[i]:+.6f}*{fac}")
    return "L(x) = " + " ".join(out)

def divided_differences(xs, ys):
    n = len(xs)
    table = [list(ys)]
    for k in range(1, n):
        prev = table[-1]
        table.append([(prev[i + 1] - prev[i]) / (xs[i + k] - xs[i]) for i in range(n - k)])
    return table

def newton_value(xs, table, x):
    n = len(xs)
    total = 0.0
    for i in range(n):
        prod = 1.0
        for j in range(i):
            prod *= (x - xs[j])
        total += table[i][0] * prod
    return total

def newton_str(xs, table):
    s = f"P(x) = {table[0][0]:.6f}"
    for k in range(1, len(xs)):
        fac = "".join(f"(x - {xs[j]:.6f})" for j in range(k))
        s += f" {table[k][0]:+.6f}*{fac}"
    return s

def print_table(xs, table):
    n = len(xs)
    head = f"{'i':>2} {'x_i':>10} {'f_i':>11}" + "".join(f" {'f['+str(k)+']':>11}" for k in range(1, n))
    print(head)
    for r in range(2 * n - 1):
        line = ""
        if r % 2 == 0:
            i = r // 2
            line += f"{i:>2} {xs[i]:>10.6f} {table[0][i]:>11.6f}"
        else:
            line += f"{'':>2} {'':>10} {'':>11}"
        for k in range(1, n):
            i = (r - k) // 2
            if r >= k and (r - k) % 2 == 0 and i < n - k:
                line += f" {table[k][i]:>11.6f}"
            else:
                line += f" {'':>11}"
        print(line)

def error_estimate(xs, x, n):
    a, b = min(xs), max(xs)
    m = 0.0
    steps = 10000
    for i in range(steps + 1):
        m = max(m, abs(f_deriv(n + 1, a + (b - a) * i / steps)))
    omega = 1.0
    for xi in xs:
        omega *= x - xi
    return m / math.factorial(n + 1) * abs(omega)

def pi_label(x):
    r = Fraction(x / math.pi).limit_denominator(12)
    if r == 0:
        return "0"
    sign = "-" if r < 0 else ""
    r = abs(r)
    num = "pi" if r.numerator == 1 else f"{r.numerator}pi"
    return sign + (num if r.denominator == 1 else f"{num}/{r.denominator}")


def plot(xs, ys, xstar, ptitle, fname, table):
    lo, hi = min(xs), max(xs)
    fig, ax1 = plt.subplots(figsize=(10, 6))

    def set_pi_ticks(ax, a, b):
        step = math.pi / 12
        k0, k1 = math.ceil(a / step), math.floor(b / step)
        ticks = [k * step for k in range(k0, k1 + 1)]
        ax.set_xticks(ticks)
        ax.set_xticklabels([pi_label(t) for t in ticks], fontsize=9)
        ax.set_xlim(a, b)

    a1, b1 = lo - 0.6, hi + 0.6
    g1 = [a1 + (b1 - a1) * i / 600 for i in range(601)]
    ax1.plot(g1, [f(t) for t in g1], color="0.55", lw=6, alpha=0.6, label="f(x) = sin x + x")
    ax1.plot(g1, [lagrange_value(xs, ys, t) for t in g1], "r--", lw=2, label="Лагранж $L(x)$")
    ax1.plot(g1, [newton_value(xs, table, t) for t in g1], "b:", lw=2.5, label="Ньютон $P(x)$")
    ax1.axvspan(lo, hi, color="green", alpha=0.07, label="отрезок узлов [x₀, xₙ]")
    ax1.plot(xs, ys, "o", color="green", ms=9, zorder=5, label="узлы интерполяции")
    ax1.axvline(xstar, color="m", ls="-.", lw=1)
    ax1.plot([xstar], [f(xstar)], "*", color="m", ms=15, zorder=6,
             label=f"$x^*$ = {xstar}, f($x^*$) = {f(xstar):.5f}")
    ax1.annotate(f"L($x^*$) = {lagrange_value(xs, ys, xstar):.5f}\n"
                 f"f($x^*$) = {f(xstar):.5f}",
                 (xstar, f(xstar)), textcoords="offset points", xytext=(25, -45),
                 arrowprops=dict(arrowstyle="->"), fontsize=9)
    set_pi_ticks(ax1, a1, b1)
    ax1.set_title(ptitle)
    ax1.set_ylabel("y")
    ax1.grid(True)
    ax1.legend(loc="upper left", fontsize=9)

    fig.tight_layout()
    fig.savefig(fname, dpi=130)
    plt.close(fig)

def run_set(name, xs, xstar):
    n = len(xs) - 1
    ys = [f(x) for x in xs]
    print(f"'{name})'")

    print("\nузлы интерполяции и значения функции")
    print(f"{'i':>2} {'x_i':>12} {'y_i':>12}")
    for i, (x, y) in enumerate(zip(xs, ys)):
        print(f"{i:>2} {x:>12.6f} {y:>12.6f}")

    print("\nЛагранж")
    print("   " + lagrange_str(xs, ys))

    print("\nразделённые разности")
    table = divided_differences(xs, ys)
    print_table(xs, table)

    print("\nНьютон")
    print("   " + newton_str(xs, table))

    L = lagrange_value(xs, ys, xstar)
    P = newton_value(xs, table, xstar)
    exact = f(xstar)
    print(f"\nx* = {xstar}")
    print(f"L(x*)          = {L:.10f}")
    print(f"P(x*)          = {P:.10f}")
    print(f"f(x*)          = {exact:.10f}")
    print(f"|f(x*)-L(x*)|  = {abs(exact - L):.3e}   (Лагранж)")
    print(f"|f(x*)-P(x*)|  = {abs(exact - P):.3e}   (Ньютон)")
    print(f"|L(x*)-P(x*)|  = {abs(L - P):.3e}")
    est = error_estimate(xs, xstar, n)
    print(f" погрешность: {est:.3e}  "
          f"({'выполняется' if abs(exact - L) <= est + 1e-15 else 'НЕ выполняется'})")

    fname = f"plot_{name}.png"
    plot(xs, ys, xstar, f"Вариант 12, y = sin(x) + x, набор узлов ({name})", fname,table)
    return xs, ys, table, L, P

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "input.txt"
    sets, xstar = read_input(path)
    for name, xs in sets.items():
        xs_, ys, table, L, P = run_set(name, xs, xstar)
        print()


if __name__ == "__main__":
    main()