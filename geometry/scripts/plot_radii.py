import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'legend.fontsize': 'x-large',
'legend.title_fontsize': 'x-large',
'axes.labelsize': 'xx-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

FIGURES_DIR = "../figures"
RESULTS_DIR = "../results"

regions = ['corona', 'disk', 'torus']

cm_to_pc = 3.24078e-19
s_to_yr = 3.17098e-8
solar_mass_to_g = 1.989e33
yr_to_s = 3.15576e7

c = 2.99792458e10 # cm / s
G = 6.67430e-8 # cm^3 / g / s^2

mu_H = 2.34e-24 # Mean mass per hydrogen nucleus; g

# Benchmark model 
L_bol = 1e45 # erg / s
L_disk = 0.5 * L_bol # erg / s
M = 1e8 * solar_mass_to_g
M_dot_w = 0.1 * solar_mass_to_g / yr_to_s # g / s 
n_ISM = 1e4 # cm^-3
v_w = 0.2 * c 

# ----------------------------------------------------------------------------------------------------
def region_radius_pc(region):

    if region == 'corona': # Where does this value come from?
        return 30 * schwarzschild_radius() * cm_to_pc
    
    elif region == 'disk': # Outer radius
        return 500 * schwarzschild_radius() * cm_to_pc
    
    elif region == 'torus': # Inner boundary
        return 2.5e18 * np.sqrt(L_disk / 1e45) * cm_to_pc # R_IR, Equation (C7)

# ----------------------------------------------------------------------------------------------------
def schwarzschild_radius():

    return 2 * G * M / c**2 # cm 

# ----------------------------------------------------------------------------------------------------
def kinetic_luminosity_wind(): # L_kin

    return M_dot_w * v_w**2 / 2 # erg / s

# ----------------------------------------------------------------------------------------------------
def free_expansion_time():

    return np.sqrt(kinetic_luminosity_wind() / (2 * np.pi * mu_H * n_ISM * v_w**5)) * s_to_yr 

# ----------------------------------------------------------------------------------------------------
def plot_radii():

    data_R_fs = np.loadtxt(f"{RESULTS_DIR}/radius_forward_shock.dat")
    data_R_ts = np.loadtxt(f"{RESULTS_DIR}/radius_termination_shock.dat")

    plt.axhline(y = region_radius_pc('corona'), color = 'gray', linewidth = 0.75)
    plt.axhline(y = region_radius_pc('disk'), color = 'gray', linewidth = 0.75)
    plt.axhline(y = region_radius_pc('torus'), color = 'gray', linewidth = 0.75)

    plt.text(x = 1.375e4, y = region_radius_pc('corona') * 1.25, s = 'Corona', color = 'gray', fontsize = 'large', ha = 'right')
    plt.text(x = 1.375e4, y = region_radius_pc('disk') * 1.25, s = 'Accretion disk', color = 'gray', fontsize = 'large', ha = 'right')
    plt.text(x = 1.375e4, y = region_radius_pc('torus') * 1.25, s = 'Dust torus', color = 'gray', fontsize = 'large', ha = 'right')

    plt.plot(data_R_fs[:,0], data_R_fs[:,1], c = '#0FE200', label = 'Forward shock')
    plt.plot(data_R_ts[:,0], data_R_ts[:,1], c = '#E20088', label = 'Termination shock')

    ax = plt.gca()
    plt.axvspan(ax.get_xlim()[0], 100 * free_expansion_time(), color = 'gray', alpha = 0.25)
    plt.text(100 * free_expansion_time() * 0.95, ax.get_ylim()[1] * 1.25, r'$100\,t_{\rm free}$', color = 'gray', ha = 'right', va = 'top')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$\rm Time \: [yr]$')
    plt.ylabel(r'$\rm Radius \: [pc]$')
    plt.legend(loc = 'lower left',frameon = False)
    plt.savefig(f"{FIGURES_DIR}/shock_radii.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/shock_radii.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_radii()

# ----------------------------------------------------------------------------------------------------
