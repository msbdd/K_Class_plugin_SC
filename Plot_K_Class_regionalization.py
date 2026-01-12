import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter


def format_ticks(x, p):
    if x >= 1:
        return f'{int(x)}'
    return f'{x:.1f}'


def B(R, l1=75.0, l2=264.0, l3=800.0,
      a1=2.11, a2=1.1, a3=2.98, a4=0.0,
      b1=1.32, b2=3.21, b3=-1.34, b4=8.0):

    R = np.asarray(R)
    result = np.zeros_like(R, dtype=float)

    mask1 = R <= l1
    mask2 = (R > l1) & (R <= l2)
    mask3 = (R > l2) & (R <= l3)
    mask4 = R > l3

    result[mask1] = a1 * np.log10(R[mask1]) + b1
    result[mask2] = a2 * np.log10(R[mask2]) + b2
    result[mask3] = a3 * np.log10(R[mask3]) + b3
    result[mask4] = a4 * np.log10(R[mask4]) + b4

    return result


def plot_B(R_min=1, R_max=1000, num_points=1000,
           filename="./K_Class/descriptions/B_R_plot.png",
           **kwargs):

    R = np.linspace(R_min, R_max, num_points)

    B_values = B(R, **kwargs)

    l1 = kwargs.get('l1', 75.0)
    l2 = kwargs.get('l2', 264.0)
    l3 = kwargs.get('l3', 800.0)

    plt.figure(figsize=(10, 6))
    plt.plot(R, B_values, 'b-', linewidth=2, label='B(R)')

    plt.axvline(x=l1, color='r', linestyle='--', alpha=0.7,
                label=f'$l_1$ = {l1} km')
    plt.axvline(x=l2, color='g', linestyle='--', alpha=0.7,
                label=f'$l_2$ = {l2} km')
    plt.axvline(x=l3, color='orange', linestyle='--', alpha=0.7,
                label=f'$l_3$ = {l3} km')

    plt.xlabel('Epicentral Distance R (km)', fontsize=12)
    plt.ylabel('B(R)', fontsize=12)
    plt.title('Piecewise distance-correction term B(R)', fontsize=14,
              fontweight='bold')
    plt.xscale('log')
    plt.grid(True, alpha=0.3)
    plt.legend(loc='best')

    ax = plt.gca()
    if ax.get_xscale() == 'log':
        ax.xaxis.set_major_formatter(FuncFormatter(format_ticks))

    plt.tight_layout()

    plt.savefig(filename, dpi=300, bbox_inches='tight')
    print(f"Plot saved as '{filename}'")

    plt.show()


if __name__ == "__main__":
    # ?TODO: Use parameters from the SeisComP config?
    plot_B()
