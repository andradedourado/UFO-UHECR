import numpy as np

s_to_yr = 3.17098e-8
solar_mass_to_g = 1.98847e33
yr_to_s = 3.15576e7

c = 2.99792458e10 # cm / s
m_p = 1.6726e-24 # g

# Benchmark model
M_dot_w = 0.1 * solar_mass_to_g / yr_to_s # g / s
mu = 1 # Mean molecular weight
n_ISM = 1e4 # cm^-3
v_w = 0.2 * c 

# ----------------------------------------------------------------------------------------------------
def deceleration_time(): # s

    rho_ISM = mu * m_p * n_ISM
    return np.sqrt(3 * M_dot_w / (4 * np.pi * v_w**3 * rho_ISM))

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    print(deceleration_time() * s_to_yr)

# ----------------------------------------------------------------------------------------------------