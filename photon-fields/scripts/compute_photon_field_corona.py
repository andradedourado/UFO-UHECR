from scipy.integrate import quad
import numpy as np 

RESULTS_DIR = "../results"

keV_to_eV = 1e3

h = 4.136e-15 # eV * s

b = 0.9
E_break = 1 * keV_to_eV # eV
E_crit = 500 * keV_to_eV # eV

L_X = pow(10, 43.8) # erg / s

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
    L[mask] = (E[mask] / E_break)**-b * np.exp(-E[mask] / E_crit)

    return L_X * L / quad(integrand_photon_field, 2 * keV_to_eV, 10 * keV_to_eV)[0] * h # erg / s / Hz

# ----------------------------------------------------------------------------------------------------
def write_photon_field():

    E = np.logspace(3, 7, num = 100) # eV
    np.savetxt(f"{RESULTS_DIR}/photon_field_corona.dat", np.column_stack((E / h, compute_photon_field(E))), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    write_photon_field()

# ----------------------------------------------------------------------------------------------------