import matplotlib.pyplot as plt
import numpy as np

FIGURES_DIR = "../figures"
REFERENCES_DIR = "../references"
RESULTS_DIR = "../results"

plt.rcParams.update({'legend.fontsize': 'x-large',
'legend.title_fontsize': 'x-large',
'axes.labelsize': 'xx-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

# ----------------------------------------------------------------------------------------------------
def get_region_color(region):

    if region == 'torus':
        return '#f17b74'

    elif region == 'disk':
        return '#625b95'

    elif region == 'corona':
        return '#f5b237'

# ----------------------------------------------------------------------------------------------------
def plot_photon_fields(region):

    DE_data = np.loadtxt(f"{REFERENCES_DIR}/Ehlert2025_photon_field_{region}.dat")
    LAD_data = np.loadtxt(f"{RESULTS_DIR}/photon_field_{region}.dat")

    plt.plot(np.log10(LAD_data[:,0]), LAD_data[:,0] * LAD_data[:,1], c = get_region_color(region), label = 'LAD')
    plt.plot(DE_data[:,0], DE_data[:,1], c = get_region_color(region), ls = '--', label = 'DE')
    plt.yscale('log')
    plt.ylim([1e41, 1e45])    
    plt.xlabel(r'$\rm \log_{10}{(Frequency / Hz)}$')
    plt.ylabel(r'$\nu L_{\nu} \: \rm [erg / s]$')
    plt.legend(title = 'Results')
    plt.savefig(f"{FIGURES_DIR}/photon_field_{region}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/photon_field_{region}.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    for region in ['disk', 'corona', 'torus']:
        plot_photon_fields(region)

# ----------------------------------------------------------------------------------------------------