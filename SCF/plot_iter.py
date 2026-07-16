"""
Function makes two plots next to each other: Change of dE and dRMS
over the SCF cycles (left) and change of the total E for n last 
iterations (right).
"""

import matplotlib.pyplot as plt
import numpy as np

def plot_iter(iterations, e_conv, d_conv, e_list, last_iter):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    dist = np.linspace(1, iterations, iterations)
    ax1.plot(dist, e_conv, color='red', linestyle='-', label='$\\log{\\Delta E}$')
    ax1.plot(dist, d_conv, color='blue', linestyle='-', label='$\\log{\\Delta\\text{RMS}}$')
    ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.7)
    ax1.legend()
    ax1.set_title("Convergence Log")
    ax1.set_xlabel("Iterations")

    ax2.plot(dist[-last_iter:], e_list[-last_iter:], color='red', linestyle='-', label='$E$')
    ax2.legend()
    ax2.set_title("Total Energy (Last Iterations)")
    ax2.set_xlabel("Iterations")

    plt.tight_layout()
    return plt.show()

