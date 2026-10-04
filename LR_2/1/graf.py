import matplotlib.pyplot as plt
import numpy as np

def main():
    x = np.linspace(-2, 2, 500)
    y1 = 3**x
    y2 = 5 * x**2 - 1
    plt.figure(figsize=(10, 6))

    plt.plot(x, y1, label="y = 3^x", color="#d62728", linewidth=2.5)
    plt.plot(x, y2, label="y = 5x^2 - 1", color="#1f77b4", linewidth=2.5)

    plt.scatter([-0.547, 0.826], [3**-0.547, 3**0.826], 
                color='gray', s=80, zorder=5, label='корни уравнения')

    plt.axhline(0, color='black', linewidth=1.2)
    plt.axvline(0, color='black', linewidth=1.2)
    plt.grid(True, which='both', linestyle='-', linewidth=0.5, color='gray', alpha=0.5)

    plt.xlim(-2, 2)
    plt.ylim(-1.5, 4.5) 

    plt.legend(fontsize=12)
    plt.title("3^x = 5x^2 - 1", fontsize=14)

    out_png = 'graf.png'
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    plt.show()

if __name__ == "__main__":
    main()