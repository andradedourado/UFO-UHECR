import matplotlib.pyplot as plt
import numpy as np

FIGURES_DIR = "../figures"
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
def get_region_label(region):

    if region == 'torus':
        return 'Dust torus'

    elif region == 'disk':
        return 'Accretion disk'

    elif region == 'corona':
        return 'Corona'

# ----------------------------------------------------------------------------------------------------
def plot_photon_fields():

    for region in ['torus', 'disk', 'corona']:
        data = np.loadtxt(f"{RESULTS_DIR}/photon_field_density_{region}.dat")
        plt.plot(np.log10(data[:,0]), data[:,1], c = get_region_color(region), label = get_region_label(region))
    
    plt.yscale('log')
    plt.ylim([1e2, 1e12])
    plt.xlabel(r'$\rm \log_{10}$(Photon energy/eV)')
    plt.ylabel(r'Number density$\rm \: [cm^{-3} eV^{-1}]$')
    plt.legend()
    plt.savefig(f"{FIGURES_DIR}/photon_fields_density.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/photon_fields_density.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_photon_fields()

# ----------------------------------------------------------------------------------------------------