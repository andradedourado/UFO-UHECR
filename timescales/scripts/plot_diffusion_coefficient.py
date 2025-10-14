import matplotlib.pyplot as plt
import numpy as np

REFERENCES_DIR = "../references"

plt.rcParams.update({'legend.fontsize': 'x-large',
'legend.title_fontsize': 'x-large',
'axes.labelsize': 'xx-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

mG_to_nG = 1e6
pc_to_cm = 3.085677581e18
yr_to_s = 3.15576e7 

c = 2.99792458e10 # cm / s

# Magnetic field
B_2 = 85 * mG_to_nG
l_c = 0.01 # pc
dlt = 5/3

Z = 7

# ----------------------------------------------------------------------------------------------------
def characteristic_diffusion_energy(): # B in nG, l_c in pc 

    return (l_c / (2 * np.pi)) * (Z * B_2 / 1.081e6) # EeV

# ----------------------------------------------------------------------------------------------------
def diffusion_coefficient(E):

    E_diff = characteristic_diffusion_energy()
    return c * l_c * pc_to_cm / (6 * np.pi) * ((E / E_diff)**(2 - dlt) + (E / E_diff) / 2 + (2/3) * (E / E_diff)**2) # cm^2 / s

# ----------------------------------------------------------------------------------------------------
def plot_diffusion_coefficent():

    E = np.logspace(-4, 4, num = 100)
    plt.plot(np.log10(E * 1e18), diffusion_coefficient(E), label = 'LAD')

    DE_data = np.loadtxt(f"{REFERENCES_DIR}/DE_down_diff_coeff_protons.dat") # pc^2 / yr
    plt.plot(DE_data[:,0], DE_data[:,1] * pc_to_cm**2 / yr_to_s, label = 'DE') 

    plt.yscale('log')
    plt.xlabel(r'$\rm \log_{10}{(Energy / eV)}$')
    plt.ylabel(r'Diffusion coefficient$\rm \: [cm^2/s]$')
    plt.legend(title = 'Results')
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_diffusion_coefficent()

# ----------------------------------------------------------------------------------------------------