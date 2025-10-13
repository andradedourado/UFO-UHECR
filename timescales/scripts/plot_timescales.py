from matplotlib import lines
import matplotlib.pyplot as plt
import numpy as np 

FIGURES_DIR = "../figures"
REFERENCES_DIR = "../references"
RESULTS_DIR = "../results"

Mpc_to_cm = 3.0857e24
s_to_yr = 3.1688e-8 # yr

c = 2.99792458e10 # cm / s

plt.rcParams.update({'legend.fontsize': 'x-large',
'legend.title_fontsize': 'x-large',
'axes.labelsize': 'xx-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

# ----------------------------------------------------------------------------------------------------
def get_timescale_color(timescale):
    
    if timescale == 'acceleration':
        return '#1984C5'

    elif timescale == 'advection':
        return '#B3BFD1'

    elif timescale == 'diffusion':
        return '#FFD700'

    elif timescale == 'ufo':
        return '#2E8B57'

# ----------------------------------------------------------------------------------------------------
def get_timescale_label(timescale):
    
    if timescale == 'acceleration':
        return 'Acceleration'

    elif timescale == 'advection':
        return 'Advection'

    elif timescale == 'diffusion':
        return 'Diffusion'

    elif timescale == 'ufo':
        return 'UFO size'
    
# ----------------------------------------------------------------------------------------------------
def plot_timescales():

    for timescale in ['acceleration', 'advection', 'diffusion', 'ufo']:
        LAD_data = np.loadtxt(f"{RESULTS_DIR}/{timescale}_timescale.dat")
        DE_data = np.loadtxt(f"{REFERENCES_DIR}/Ehlert2025_{timescale}_timescale.dat")

        plt.plot(np.log10(LAD_data[:,0]), LAD_data[:,1], c = get_timescale_color(timescale), label = get_timescale_label(timescale))
        plt.plot(DE_data[:,0], DE_data[:,1] * Mpc_to_cm / c * s_to_yr, c = get_timescale_color(timescale), ls = '--')

    LAD_label = lines.Line2D([], [], c = 'k', label = 'LAD')
    DE_label = lines.Line2D([], [], c = 'k', ls = '--', label = 'DE')
    lgnd = plt.legend(handles = [LAD_label, DE_label], frameon = True, loc = 'lower left', bbox_to_anchor = (0, 0.36))
    plt.gca().add_artist(lgnd)
    
    plt.yscale('log')
    plt.ylim([1e-1, 1e5])
    plt.xlabel(r'$\rm \log_{10}{(Energy / eV)}$')
    plt.ylabel(r'$\rm Timescales \: [yr]$')
    plt.legend(loc = 'lower left')
    plt.savefig(f"{FIGURES_DIR}/timescales.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/timescales.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_timescales()

# ----------------------------------------------------------------------------------------------------
