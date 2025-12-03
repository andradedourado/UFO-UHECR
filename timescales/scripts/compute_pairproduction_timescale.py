import scipy.constants as const
from scipy.integrate import quad
from scipy.interpolate import interp1d
from scipy.special import zeta
import numpy as np

PHOTON_FIELDS_RESULTS_DIR = "../../photon-fields/results"
RESULTS_DIR = "../results"

PARTICLES = ['1H', '4He', '14N', '28Si', '56Fe']
ZS = [1, 2, 7, 14, 26]

cm_to_m = 1e-2
eV_to_J = 1.60218e-19
s_to_yr = 1 / (60 * 60 * 24 * 365.25)

m_p = 0.9383e9 # eV
r0 = (const.e**2) / (4 * const.pi * const.epsilon_0 * const.m_e * const.c**2)

# ----------------------------------------------------------------------------------------------------
def iZ(Z):

    try:
        return ZS.index(Z)
    except ValueError:
        raise ValueError(f"Z ({Z}) not found in ZS.")

# ----------------------------------------------------------------------------------------------------
def photon_density(eps, region): # Number of photons per unit volume and energy

    data = np.loadtxt(f"{PHOTON_FIELDS_RESULTS_DIR}/photon_field_density_{region}.dat") # eV, cm^-3 eV^-1 

    E = data[:,0] * eV_to_J # SI 
    dn_dE = data[:,1] / cm_to_m**3 / eV_to_J # SI

    return interp1d(E, dn_dE, kind = 'linear', bounds_error = False, fill_value = 0.0)(eps)

# ----------------------------------------------------------------------------------------------------
def phi(kappa):

    c1 = 0.8048
    c2 = 0.1459
    c3 = 1.137 * 1.e-3
    c4 = -3.879 * 1.e-6
    c = [c1, c2, c3, c4]

    d0 = -170 + 84*np.log(2) - 16*np.log(2)**2 + np.pi**2/3 * (10 - 4*np.log(2)) + 8*zeta(3)
    d1 = 88 - 40*np.log(2) + 8*np.log(2)**2 - 4/3*np.pi**2
    d2 = -20 + 8*np.log(2)
    d3 = 8/3
    d = [d0, d1, d2, d3]

    f1 = 2.910
    f2 = 78.35
    f3 = 1837
    f = [f1, f2, f3]

    if kappa < 25:

        sum_c_term = 0

        for i in range(len(c)):
            sum_c_term += c[i] * (kappa - 2)**(i + 1)

        return np.pi / 12 * (kappa - 2)**4 / (1 + sum_c_term)
    
    elif kappa >= 25:

        sum_d_term = 0
        sum_f_term = 0 

        for i in range(len(d)):
            sum_d_term += d[i] * np.log(kappa)**i 
        
        for i in range(len(f)):
            sum_f_term += f[i] * kappa**-(i + 1)

        return kappa * sum_d_term / (1 - sum_f_term)

# ----------------------------------------------------------------------------------------------------
def integrand_interaction_rate(kappa, Gmm, region):

    eps = kappa * const.m_e * const.c**2 / (2 * Gmm)
    return photon_density(eps, region) * phi(kappa) / kappa**2

# ----------------------------------------------------------------------------------------------------
def compute_pairproduction_timescales(A, Z, Gmms, region):

    timescales = []

    for Gmm in Gmms:
        interaction_rate = const.alpha * r0**2 * const.c * Z**2 * const.m_e / (A * const.m_p) / Gmm
        interaction_rate = interaction_rate * quad(integrand_interaction_rate, 2, 1e4, args = (Gmm, region))[0] * const.m_e * const.c**2
        timescales.append(interaction_rate**-1 * s_to_yr)
    
    return np.array(timescales)

# ----------------------------------------------------------------------------------------------------
def write_pairproduction_timescales(A, Z, region):

    E = np.logspace(15, 21, num = 100)
    Gmms = E / (A * m_p) 
    timescales = compute_pairproduction_timescales(A, Z, Gmms, region)
    np.savetxt(f"{RESULTS_DIR}/pairproduction_timescale_{PARTICLES[iZ(Z)]}_{region}.dat", np.column_stack((E, timescales)), fmt = "%.15e", delimiter = "\t")

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    for region in ['corona', 'disk', 'torus']:
        write_pairproduction_timescales(14, 7, region)

# ----------------------------------------------------------------------------------------------------