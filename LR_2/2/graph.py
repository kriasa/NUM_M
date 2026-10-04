import numpy as np
import matplotlib.pyplot as plt

a = 3.0
X1_MIN, X1_MAX = -3.0, 7.0
X2_MIN, X2_MAX = -3.0, 7.0
N = 2000


def main():
    x2_grid = np.linspace(X2_MIN, X2_MAX, N)
    x1_of_f1 = a + np.cos(x2_grid)
    x1_grid = np.linspace(X1_MIN, X1_MAX, N)
    x2_of_f2 = a + np.sin(x1_grid)

    fig, ax = plt.subplots(figsize=(7, 7))

    ax.plot(x1_of_f1, x2_grid, 'b-', lw=1.6,
            label=r'$f_1(x_1,x_2)=0:\ \ x_1-\cos x_2-a=0$')
    ax.plot(x1_grid, x2_of_f2, 'r-', lw=1.6,
            label=r'$f_2(x_1,x_2)=0:\ \ x_2-\sin x_1-a=0$')

    ax.axhline(0, color='gray', lw=0.8)
    ax.axvline(0, color='gray', lw=0.8)
    ax.axhspan(0, X2_MAX, xmin=(0 - X1_MIN) / (X1_MAX - X1_MIN),
               color='green', alpha=0.04)

    ax.set_xlim(X1_MIN, X1_MAX)
    ax.set_ylim(X2_MIN, X2_MAX)
    ax.set_xlabel('$x_1$')
    ax.set_ylabel('$x_2$')
    ax.set_title(f'Вариант 12 (a = {a:g}): локализация корней системы')
    ax.grid(True, lw=0.4)
    ax.legend(loc='upper left')
    ax.set_aspect('equal', adjustable='box')

    out_png = 'graph.png'
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    try:
        plt.show()
    except Exception:
        pass


if __name__ == "__main__":
    main()