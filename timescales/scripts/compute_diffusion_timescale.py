from scipy.interpolate import interp1d
import numpy as np

GEOMETRY_RESULTS_DIR = "../../geometry/results"
RESULTS_DIR = "../results"

mG_to_nG = 1e6
pc_to_cm = 3.085677581e18
s_to_yr = 3.168808781e-8

c = 2.99792458e10 # cm / s

# UFO
t_age = 1e3 # yr

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
def forward_shock_radius(t): # t in yr

    data = np.loadtxt(f"{GEOMETRY_RESULTS_DIR}/radius_forward_shock.dat")
    return interp1d(data[:,0], data[:,1], kind = 'linear')(t) * pc_to_cm

# ----------------------------------------------------------------------------------------------------
def termination_shock_radius(t): # t in yr

    data = np.loadtxt(f"{GEOMETRY_RESULTS_DIR}/radius_termination_shock.dat")
    return interp1d(data[:,0], data[:,1], kind = 'linear')(t) * pc_to_cm

# ----------------------------------------------------------------------------------------------------
def compute_diffusion_timescales(E):

    R_fs = forward_shock_radius(t_age)
    R_sh = termination_shock_radius(t_age)

    return (R_fs - R_sh)**2 / diffusion_coefficient(E) * s_to_yr

# ----------------------------------------------------------------------------------------------------
def write_diffusion_timescales():

    E = np.logspace(-3, 3, num = 100)
    np.savetxt(f"{RESULTS_DIR}/diffusion_timescale.dat", np.column_stack((E * 1e18, compute_diffusion_timescales(E))), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    write_diffusion_timescales()

# ----------------------------------------------------------------------------------------------------