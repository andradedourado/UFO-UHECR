from scipy.integrate import quad
from scipy.interpolate import interp1d
import numpy as np

RESULTS_DIR = "../results"

eV_to_erg = 1.602176634e-12
solar_mass_to_g = 1.989e33

c = 2.99792458e10 # cm / s
G = 6.67430e-8 # m^3 / g / s^2
h = 6.62607015e-27 # erg * s
k = 1.380649e-16 # erg / K
sigma = 5.670374419e-5 # erg / s / cm^2 / K^4
eta = 0.1

# Benchmark model 
L_bol = 1e45 # erg / s
L_disk = 0.5 * L_bol # erg / s
M = 1e8 * solar_mass_to_g

# ----------------------------------------------------------------------------------------------------
def schwarzschild_radius():
    
    return 2 * G * M / c**2 # cm 

# ----------------------------------------------------------------------------------------------------
def temperature(R):

    R_S = schwarzschild_radius()
    return (3 * R_S * L_disk / (16 * np.pi * eta * sigma * R**3) * (1 - (3 * R_S / R)**(1/2)))**(1/4) # K

# ----------------------------------------------------------------------------------------------------
def blackbody_radiance_nu(nu, R): # Spectral radiance: energy emitted per unit area, per unit solid angle, per unit frequency, and per unit time

    T = temperature(R)
    return 2 * h * nu**3 / c**2 / (np.exp((h * nu) / (k * T)) - 1) # erg / s / cm^2 / sr / Hz

# ----------------------------------------------------------------------------------------------------
def integrand_mono_disk_lum(R, nu): # Equation (C3)

    return 4 * np.pi**2 * R * blackbody_radiance_nu(nu, R)
 
# ----------------------------------------------------------------------------------------------------
def compute_photon_field(nu_array):

    L_disk = np.zeros_like(nu_array)
    R_S = schwarzschild_radius()

    for inu, nu in enumerate(nu_array):
        L_disk[inu] = quad(integrand_mono_disk_lum, 3 * R_S, 500 * R_S, args = (nu))[0]
    
    f = interp1d(nu_array, L_disk, kind = 'linear')
    nu_B = 6.813e14 # Hz, λ = 4400 Å, B band (blue)
    L_nu_B = f(nu_B)

    return (L_bol / 5.13) * L_disk / (nu_B * L_nu_B)

# ----------------------------------------------------------------------------------------------------
def write_photon_field():

    nu_array = np.logspace(13, 17, num = 100) # Hz    
    np.savetxt(f"{RESULTS_DIR}/photon_field_disk.dat", np.column_stack((nu_array, compute_photon_field(nu_array))), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    write_photon_field()

# ----------------------------------------------------------------------------------------------------