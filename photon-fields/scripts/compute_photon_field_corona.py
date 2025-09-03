import numpy as np 

RESULTS_DIR = "../results"

keV_to_eV = 1e3

b = 0.9
E_break = 1 * keV_to_eV # eV
E_crit = 500 * keV_to_eV # eV

# ----------------------------------------------------------------------------------------------------
def compute_photon_field(E):

    L = np.zeros_like(E)

    mask = E >= E_break
    L[mask] = (E[mask] / E_break)**-b * np.exp(-E[mask] / E_crit)

    return L

# ----------------------------------------------------------------------------------------------------
def write_photon_field():

    E = np.logspace(3, 7, num = 100) # eV
    np.savetxt(f"{RESULTS_DIR}/photon_field_corona.dat", np.column_stack((E, compute_photon_field(E))), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    write_photon_field()

# ----------------------------------------------------------------------------------------------------