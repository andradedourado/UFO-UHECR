import numpy as np

RESULTS_DIR = "../results"

km_to_cm = 1e5 
solar_mass_to_g = 1.98847e33
Myr_to_yr = 1e6
pc_to_km = 3.0857e13
yr_to_s = 3.15576e7

c = 2.99792458e10 # cm / s

# Benchmark model
M_dot_w = 0.1 * solar_mass_to_g / yr_to_s # g / s 
n_ISM = 1e4 # cm^-3
v_w = 0.2 * c 

# ----------------------------------------------------------------------------------------------------
def kinetic_energy_wind(): # L_kin 

    return M_dot_w * v_w**2 / 2 # erg / s

# ----------------------------------------------------------------------------------------------------
def compute_v_fs(t_age): # Forward shock, [t_age] = yr 

    L_kin = kinetic_energy_wind()
    return 76 * pc_to_km / (Myr_to_yr * yr_to_s) * pow(t_age / Myr_to_yr, -2/5) * pow(L_kin / 1e38, 1/5) * pow(n_ISM, -1/5) # km / s

# ----------------------------------------------------------------------------------------------------
def compute_v_sh(t_age): # Termination shock, [t_age] = yr 

    L_kin = kinetic_energy_wind()
    return 23 * pc_to_km / (Myr_to_yr * yr_to_s) * pow(t_age / Myr_to_yr, -3/5) * pow(L_kin / 1e38, 3/10)* pow(n_ISM, -3/10) * pow(v_w / (1e3 * km_to_cm), -1/2) # km / s
    
# ----------------------------------------------------------------------------------------------------
def write_velocities():

    t_age = np.logspace(0, 4, num = 100)
    np.savetxt(f"{RESULTS_DIR}/velocity_forward_shock.dat", np.column_stack((t_age, compute_v_fs(t_age))), fmt = "%.15e")
    np.savetxt(f"{RESULTS_DIR}/velocity_termination_shock.dat", np.column_stack((t_age, compute_v_sh(t_age))), fmt = "%.15e")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    write_velocities()

# ----------------------------------------------------------------------------------------------------