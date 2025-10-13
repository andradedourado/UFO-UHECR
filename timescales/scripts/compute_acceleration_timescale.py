from scipy.interpolate import interp1d
import numpy as np

GEOMETRY_RESULTS_DIR = "../../geometry/results"
RESULTS_DIR = "../results"

G_to_nG = 1e9
pc_to_cm = 3.085677581e18
s_to_yr = 3.168808781e-8
solar_mass_to_g = 1.98847e33
yr_to_s = 3.15576e7

c = 2.99792458e10 # cm / s

# UFO 
M_dot_w = 0.1 * solar_mass_to_g / yr_to_s # g / s 
v_w = 0.2 * c
t_age = 1e3 # yr

# Magnetic field
eps_B = 0.05 
l_c = 0.01 # pc
dlt = 5/3

s = 4 # Spectral index for strong shocks

Z = 7 # Nitrogen

# ----------------------------------------------------------------------------------------------------
def termination_shock_radius(t): # t in yr

    data = np.loadtxt(f"{GEOMETRY_RESULTS_DIR}/radius_termination_shock.dat")
    return interp1d(data[:,0], data[:,1], kind = 'linear')(t) * pc_to_cm

# ----------------------------------------------------------------------------------------------------
def mass_density(R):

    return M_dot_w / (4 * np.pi * R**2 * v_w) # g / cm^3

# ----------------------------------------------------------------------------------------------------
def magnetic_field_density(R):

    return eps_B * mass_density(R) * v_w**2 # erg / cm^3

# ----------------------------------------------------------------------------------------------------
def magnetic_field(R):

    return np.sqrt(8 * np.pi * magnetic_field_density(R)) * G_to_nG 

# ----------------------------------------------------------------------------------------------------
def characteristic_diffusion_energy(R): # B in nG, l_c in pc

    return (l_c / (2 * np.pi)) * (Z * magnetic_field(R) / 1.081e6) # EeV

# ----------------------------------------------------------------------------------------------------
def diffusion_coefficient(E, R):

    E_diff = characteristic_diffusion_energy(R)
    return c * l_c * pc_to_cm / (6 * np.pi) * ((E / E_diff)**(2 - dlt) + (E / E_diff) / 2 + (2/3) * (E / E_diff)**2) # cm^2 / s

# ----------------------------------------------------------------------------------------------------
def compute_acceleration_timescale(E, R):

    return s * diffusion_coefficient(E, R) / v_w**2 * s_to_yr

# ----------------------------------------------------------------------------------------------------
def write_acceleration_timescale(R):

    E = np.logspace(-3, 3, num = 100)
    np.savetxt(f"{RESULTS_DIR}/acceleration_timescale.dat", np.column_stack((E * 1e18, compute_acceleration_timescale(E, R))), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    R_sh = termination_shock_radius(t_age)
    write_acceleration_timescale(R_sh)

# ----------------------------------------------------------------------------------------------------