from scipy.interpolate import interp1d
import numpy as np

GEOMETRY_RESULTS_DIR = "../../geometry/results"
RESULTS_DIR = "../results"

pc_to_cm = 3.085677581e18
s_to_yr = 3.168808781e-8

c = 2.99792458e10 # cm / s

t_age = 1e3 # yr

# ----------------------------------------------------------------------------------------------------
def forward_shock_radius(t): # t in yr

    data = np.loadtxt(f"{GEOMETRY_RESULTS_DIR}/radius_forward_shock.dat")
    return interp1d(data[:,0], data[:,1], kind = 'linear')(t) * pc_to_cm

# ----------------------------------------------------------------------------------------------------
def compute_ufo_size(E):

    return np.full_like(E, forward_shock_radius(t_age) / c * s_to_yr, dtype = float)  

# ----------------------------------------------------------------------------------------------------
def write_ufo_size():

    E = np.logspace(15, 21, num = 100)
    np.savetxt(f"{RESULTS_DIR}/ufo_timescale.dat", np.column_stack((E, compute_ufo_size(E))), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    write_ufo_size()

# ----------------------------------------------------------------------------------------------------