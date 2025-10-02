import matplotlib.pyplot as plt
import numpy as np

FIGURES_DIR = "../figures"
RESULTS_DIR = "../results"
SIMULATIONS_DIR = "../../simulations/results"

plt.rcParams.update({'legend.fontsize': 'x-large',
'legend.title_fontsize': 'x-large',
'axes.labelsize': 'xx-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

# ----------------------------------------------------------------------------------------------------
def plot_spectrum():

    E = np.loadtxt(f"{SIMULATIONS_DIR}/shock_downstream.txt")[:,2] * 1e18

    bin_edges = 10**np.linspace(15, 21)
    bin_width = bin_edges[1:] - bin_edges[:-1]
    bin_center = bin_edges[:-1] + 0.5 * bin_width

    H = np.histogram(E, bins = bin_edges)
    J = H[0] / bin_width
    dJ = J / np.sqrt(H[0])

    plt.errorbar(bin_center, J * bin_center**2, xerr = bin_width / 2, yerr = dJ * bin_center**2, c = 'k', ls = 'None')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1e15, 1e21])
    plt.ylim([1e19, 1e21])
    plt.xlabel(r'$\rm Energy \: [eV]$')
    plt.ylabel(r'$E^2 \times {\rm Flux} \: \rm [a.u.]$')
    plt.savefig(f"{FIGURES_DIR}/spectrum.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/spectrum.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
def plot_spectrum_splitting():

    data = np.loadtxt(f"{SIMULATIONS_DIR}/shock_downstream_splitting.txt")
    E = data[:,2] * 1e18
    w = data[:,17]

    bin_edges = 10**np.linspace(15, 21)
    bin_width = bin_edges[1:] - bin_edges[:-1]
    bin_center = bin_edges[:-1] + 0.5 * bin_width

    H = np.histogram(E, bins = bin_edges, weights = w)
    H_count = np.histogram(E, bins = bin_edges)
    J = H[0] / bin_width
    dJ = J / np.sqrt(H_count[0])

    plt.errorbar(bin_center, J * bin_center**2, xerr = bin_width / 2, yerr = dJ * bin_center**2, c = 'k', ls = 'None')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim([1e15, 1e21])
    plt.ylim([1e19, 1e21])
    plt.xlabel(r'$\rm Energy \: [eV]$')
    plt.ylabel(r'$E^2 \times {\rm Flux} \: \rm [a.u.]$')
    plt.savefig(f"{FIGURES_DIR}/spectrum_splitting.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/spectrum_splitting.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_spectrum()
    plot_spectrum_splitting()

# ----------------------------------------------------------------------------------------------------
