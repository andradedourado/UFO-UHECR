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
def plot_radii():

    plt.axhline(y = 299792.458, color = 'gray', linestyle = ':')

    data_v_fs = np.loadtxt(f"{RESULTS_DIR}/velocity_forward_shock.dat")
    data_v_sh = np.loadtxt(f"{RESULTS_DIR}/velocity_termination_shock.dat")

    plt.plot(data_v_fs[:,0], data_v_fs[:,1], label = 'FS')
    plt.plot(data_v_sh[:,0], data_v_sh[:,1], label = 'SH')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$\rm Time \: [yr]$')
    plt.ylabel(r'$\rm Velocity \: [km/s]$')
    plt.legend(title = 'Shock')
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_radii()

# ----------------------------------------------------------------------------------------------------
