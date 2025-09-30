from scipy.integrate import quad
import numpy as np 

RESULTS_DIR = "../results"

eV_to_J = 1.602176634e-19

c = 2.99792458e10 # cm / s
h = 6.62607015e-27 # erg * s
k = 1.380649e-16 # erg / K

# L_bol = 1e45 # erg / s
# L_disk = 0.5 * L_bol # erg / s
L_X = pow(10, 43.8) # erg / s
T_IR = 200 # K

# ----------------------------------------------------------------------------------------------------
def normalization_correction_factor(): # Mullaney et al. (2011)

    return (np.exp(0.53) + (L_X / 1e43)**1.1) * 1e43 # erg / s 

# ----------------------------------------------------------------------------------------------------
def blackbody_radiance_nu(nu): # Spectral radiance: energy emitted per unit area, per unit solid angle, per unit frequency, and per unit time

    return 2 * h * nu**3 / c**2 / (np.exp((h * nu) / (k * T_IR)) - 1) # erg / cm^2 / sr^-1 / Hz / s

# ----------------------------------------------------------------------------------------------------
def compute_photon_field(nu):

    return normalization_correction_factor() * blackbody_radiance_nu(nu) / quad(blackbody_radiance_nu, 1e10, 1e15)[0]

# ----------------------------------------------------------------------------------------------------
def write_photon_field():

    nu = np.logspace(11, 14, num = 100) # Hz
    return np.savetxt(f"{RESULTS_DIR}/photon_field_torus.dat", np.column_stack((nu, compute_photon_field(nu))), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    write_photon_field()

# ----------------------------------------------------------------------------------------------------