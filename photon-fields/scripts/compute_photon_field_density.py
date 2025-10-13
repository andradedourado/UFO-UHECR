import numpy as np

RESULTS_DIR = "../results"

erg_to_eV = 6.241509074e11
solar_mass_to_g = 1.989e33

c = 2.99792458e10 # cm / s
G = 6.67430e-8 # m^3 / g / s^2
h = 6.62607015e-27  # erg * s

# Benchmark model 
L_bol = 1e45 # erg / s
L_disk = 0.5 * L_bol # erg / s
M = 1e8 * solar_mass_to_g

# ----------------------------------------------------------------------------------------------------
def schwarzschild_radius():
    
    return 2 * G * M / c**2 # cm 

# ----------------------------------------------------------------------------------------------------
def region_radius(region):

    if region == 'corona':
        return 30 * schwarzschild_radius()
    
    elif region == 'disk':
        return 500 * schwarzschild_radius()
    
    elif region == 'torus':
        return 2.5e18 * (L_disk / 1e45)**(1/2) # cm 

# ----------------------------------------------------------------------------------------------------
def write_photon_field_density(region): # [Number of photons per unit volume and energy] = cm^-3 eV^1
    
    data = np.loadtxt(f"{RESULTS_DIR}/photon_field_luminosity_{region}.dat")
        
    E = h * data[:,0] # erg
    nu_L_nu = data[:,0] * data[:,1]
    photon_number_density = nu_L_nu / (4 * np.pi * region_radius(region)**2 * c * E**2 * erg_to_eV)
        
    np.savetxt(f"{RESULTS_DIR}/photon_field_density_{region}.dat", np.column_stack((E * erg_to_eV, photon_number_density)), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    for region in ['corona', 'disk', 'torus']:
        write_photon_field_density(region)

# ----------------------------------------------------------------------------------------------------


