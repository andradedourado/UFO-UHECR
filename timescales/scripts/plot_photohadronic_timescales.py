from matplotlib import lines
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'legend.fontsize': 'x-large',
'legend.title_fontsize': 'x-large',
'axes.labelsize': 'xx-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

FIGURES_DIR = "../figures"
REFERENCES_DIR = "../references"
RESULTS_DIR = "../results"

PARTICLES = ['1H', '4He', '14N', '28Si', '56Fe']
ZS = [1, 2, 7, 14, 26]

Mpc_to_cm = 3.0857e24
s_to_yr = 3.1688e-8 # yr

c = 2.99792458e10 # cm / s

# ----------------------------------------------------------------------------------------------------
def iZ(Z):

    try:
        return ZS.index(Z)
    except ValueError:
        raise ValueError(f"Z ({Z}) not found in ZS.")

# ----------------------------------------------------------------------------------------------------
def get_region_color(region):

    if region == 'torus':
        return '#f17b74'

    elif region == 'disk':
        return '#625b95'

    elif region == 'corona':
        return '#f5b237'
    
# ----------------------------------------------------------------------------------------------------
def get_region_label(region):

    if region == 'torus':
        return 'Torus'

    elif region == 'disk':
        return 'Disk'

    elif region == 'corona':
        return 'Corona'

# ----------------------------------------------------------------------------------------------------
def plot_pairproduction_timescale(Z):

    for region in ['corona', 'disk', 'torus']:
        LAD_data = np.loadtxt(f"{RESULTS_DIR}/pairproduction_timescale_{PARTICLES[iZ(Z)]}_{region}.dat")
        plt.plot(np.log10(LAD_data[:,0]), LAD_data[:,1], c = get_region_color(region), label = get_region_label(region))

        if Z == 7: 
            if region != 'corona': 
                DE_data = np.loadtxt(f"{REFERENCES_DIR}/Ehlert2025_pairproduction_timescale_{region}.dat")
                plt.plot(DE_data[:,0], DE_data[:,1] * Mpc_to_cm / c * s_to_yr, c = get_region_color(region), ls = '--')

    LAD_label = lines.Line2D([], [], c = 'k', label = 'LAD')
    DE_label = lines.Line2D([], [], c = 'k', ls = '--', label = 'DE')
    lgnd = plt.legend(handles = [LAD_label, DE_label], frameon = True, loc = 'lower left')
    plt.gca().add_artist(lgnd)

    plt.yscale('log')
    # plt.ylim(top = 1e5)
    plt.xlabel(r'$\rm \log_{10}{(Energy / eV)}$')
    plt.ylabel(r'$\rm Timescales \: [yr]$')
    plt.legend()
    plt.savefig(f"{FIGURES_DIR}/pairproduction_timescale.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/pairproduction_timescale.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_pairproduction_timescale(7)

# ----------------------------------------------------------------------------------------------------