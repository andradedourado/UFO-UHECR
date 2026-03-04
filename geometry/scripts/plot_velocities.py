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

cm_to_km = 1e-5
s_to_yr = 3.17098e-8
solar_mass_to_g = 1.989e33
yr_to_s = 3.15576e7

c = 2.99792458e10 # cm / s

mu_H = 2.34e-24 # Mean mass per hydrogen nucleus; g

# Benchmark model
M_dot_w = 0.1 * solar_mass_to_g / yr_to_s # g / s 
n_ISM = 1e4 # cm^-3
v_w = 0.2 * c 

# ----------------------------------------------------------------------------------------------------
def kinetic_luminosity_wind(): # L_kin

    return M_dot_w * v_w**2 / 2 # erg / s

# ----------------------------------------------------------------------------------------------------
def free_expansion_time():

    return np.sqrt(kinetic_luminosity_wind() / (2 * np.pi * mu_H * n_ISM * v_w**5)) * s_to_yr 

# ----------------------------------------------------------------------------------------------------
def plot_velocities():

    data_v_fs = np.loadtxt(f"{RESULTS_DIR}/velocity_forward_shock.dat")
    data_v_ts = np.loadtxt(f"{RESULTS_DIR}/velocity_termination_shock.dat")

    plt.axhline(y = c * cm_to_km, color = 'gray', linewidth = 0.75)
    plt.text(x = 1.375e4, y = c * cm_to_km  / 1.25 / 1.25, s = 'Speed of light', color = 'gray', fontsize = 'large', ha = 'right')

    ax = plt.gca()
    plt.axvspan(ax.get_xlim()[0], 100 * free_expansion_time(), color = 'gray', alpha = 0.25)
    plt.text(100 * free_expansion_time() * 0.95, ax.get_ylim()[1] * 6.5e-4, r'$100\,t_{\rm free}$', color = 'gray', ha = 'right', va = 'top')

    plt.plot(data_v_fs[:,0], data_v_fs[:,1], c = '#0FE200', label = 'Forward shock')
    plt.plot(data_v_ts[:,0], data_v_ts[:,1], c = '#E20088', label = 'Termination shock')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$\rm Time \: [yr]$')
    plt.ylabel(r'$\rm Velocity \: [km/s]$')
    plt.legend(loc = 'upper left', frameon = False)
    plt.savefig(f"{FIGURES_DIR}/shock_velocities.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/shock_velocities.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
def plot_upstream_velocity(): # Termination shock

    data_v_ts = np.loadtxt(f"{RESULTS_DIR}/velocity_termination_shock.dat")

    v_up = (v_w * cm_to_km - data_v_ts[:,1]) / (1 - (v_w * cm_to_km * data_v_ts[:,1]) / (c * cm_to_km)**2) 

    plt.plot(data_v_ts[:,0], v_up, c = '#E20088')

    ax = plt.gca()
    plt.axvspan(ax.get_xlim()[0], 100 * free_expansion_time(), color = 'gray', alpha = 0.25)
    plt.text(100 * free_expansion_time() * 0.95, ax.get_ylim()[0] * 1.04, r'$100\,t_{\rm free}$', color = 'gray', ha = 'right', va = 'top')

    plt.xscale('log')
    plt.xlabel(r'$\rm Time \: [yr]$')
    plt.ylabel(r'Upstream velocity$\rm \: [km/s]$')
    plt.savefig(f"{FIGURES_DIR}/upstream_velocity.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/upstream_velocity.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_velocities()
    plot_upstream_velocity()

# ----------------------------------------------------------------------------------------------------
