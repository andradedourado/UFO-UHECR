from scipy.integrate import quad
from scipy.interpolate import interp1d
import numpy as np 

RESULTS_DIR = "../results"

keV_to_eV = 1e3

h = 4.136e-15 # eV * s

b = 0.9
E_break = 1 * keV_to_eV # eV
E_crit = 500 * keV_to_eV # eV

L_X = pow(10, 43.8) # erg / s

# ----------------------------------------------------------------------------------------------------
def photon_field_disk(E):

    data = np.loadtxt(f"{RESULTS_DIR}/photon_field_luminosity_disk.dat")
    return np.interp(E, data[:,0] * h, data[:,1])

# ----------------------------------------------------------------------------------------------------
def integrand_photon_field(E):

    if E >= E_break:
        return (E / E_break)**-b * np.exp(-E / E_crit)
    else:
        return 0

# ----------------------------------------------------------------------------------------------------
def compute_photon_field(E):

    L = np.zeros_like(E)
    mask = E >= E_break

    # Hot corona
    L[mask] = (E[mask] / E_break)**-b * np.exp(-E[mask] / E_crit)
    L = L_X * L / quad(integrand_photon_field, 2 * keV_to_eV, 10 * keV_to_eV)[0] * h # erg / s / Hz

    # Soft X-ray regime: L(E) = const * (E / E_break)**-index
    E_soft_start = 25 # eV
    L_soft_25eV = photon_field_disk(E_soft_start)
    L_Ebreak = L[mask][0]

    index = -np.log(L_soft_25eV / L_Ebreak) / np.log(E_soft_start / E_break)
    L[~mask] = L_Ebreak * (E[~mask] / E_break)**-index - photon_field_disk(E[~mask])
    return L

# ----------------------------------------------------------------------------------------------------
def write_photon_field():

    E_soft_01 = np.logspace(np.log10(25), 2, num = 50, endpoint = False)
    E_soft_02 = np.logspace(2, np.log10(E_break), num = 50, endpoint = False)
    E_hot = np.logspace(np.log10(E_break), 7, num = 50)
    E = np.concatenate((E_soft_01, E_soft_02, E_hot)) # eV

    np.savetxt(f"{RESULTS_DIR}/photon_field_luminosity_corona.dat", np.column_stack((E / h, compute_photon_field(E))), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    write_photon_field()

# ----------------------------------------------------------------------------------------------------