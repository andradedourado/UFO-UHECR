import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({'legend.fontsize': 'x-large',
'legend.title_fontsize': 'x-large',
'axes.labelsize': 'xx-large',
'axes.titlesize': 'xx-large',
'xtick.labelsize': 'x-large',
'ytick.labelsize': 'x-large'})

FIGURES_DIR = "../figures"
SIMULARIONS_DIR = "../../simulations/DiffusiveSDE/results"

pc_to_m = 3.086e16
s_to_yr = 3.1688e-8
yr_to_s = 31536000

c = 299792458 # m / s

# ----------------------------------------------------------------------------------------------------
def get_position_mask(sim):

    if sim == '3D_sph_Dconst':
        return np.inf # m
    
# ----------------------------------------------------------------------------------------------------
def get_simulation_parameters(sim):

    if sim == '3D_sph_Dconst':
        T_min = 1e-4 * yr_to_s
        T_max = 1e1 * yr_to_s
        n_time = 200
        is_log_scale = True

    return T_min, T_max, n_time, is_log_scale

# ----------------------------------------------------------------------------------------------------
def get_time_intervals(sim):

    T_min, T_max, n_time, is_log_scale = get_simulation_parameters(sim)
    t_i = np.zeros(n_time)

    if is_log_scale:
        for i in range(n_time):
            t_i[i] = T_min * (T_max / T_min)**(i / (n_time - 1))
    else:
        for i in range(n_time):
            t_i[i] = T_min + i * (T_max - T_min) / (n_time - 1)

    return t_i

# ----------------------------------------------------------------------------------------------------
def compute_energy_spectrum(candidates, E_bins, is_weighted, is_stationary):

    E = candidates[:,1]

    if is_stationary == False:
        w_T = np.ones(len(candidates))
    elif is_stationary == True:
        w_T = candidates[:,0] / c

    if is_weighted:
        weights = candidates[:,5]
        hist_counts_w, E_edges = np.histogram(E, bins = E_bins, weights = weights * w_T)
        hist_counts, _ = np.histogram(E, bins = E_bins, weights = w_T)
        spectrum = hist_counts_w / (E_edges[1:] - E_edges[:-1])
        sigma = spectrum / np.sqrt(hist_counts)
    else:
        hist_counts, E_edges = np.histogram(E, bins = E_bins, weights = w_T)
        E_width = E_edges[1:] - E_edges[:-1]
        spectrum = hist_counts / E_width
        sigma = np.sqrt(hist_counts) / E_width

    E_center = E_edges[:-1] + 0.5 * (E_edges[1:] - E_edges[:-1])

    return E_center, spectrum, sigma

# ----------------------------------------------------------------------------------------------------
def plot_energy_spectrum_over_time(sim):

    data = np.loadtxt(f"{SIMULARIONS_DIR}/candidates_{sim}.txt")

    D, E, X, Y, Z = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]
    R = np.sqrt(X**2 + Y**2 + Z**2)
    T = D / c

    E_bins = np.logspace(np.log10(min(E)), np.log10(max(E)), num = 50) 
    t_i = get_time_intervals(sim)

    for i in range(118, len(t_i) - 1, 10):
        mask = (T > t_i[i]) & (T < t_i[i+1]) # & (R > 0) & (R < get_position_mask(sim))
        candidates = data[mask,:]
        E_center, spectrum, sigma = compute_energy_spectrum(candidates, E_bins, True, False)
        plt.errorbar(E_center, spectrum * E_center**2, sigma * E_center**2, ls = 'None', marker = '.', label = f'{float(t_i[i + 1] * s_to_yr):.1e}')

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$\rm Energy \: [TeV]$')
    plt.ylabel(r'$E^2 \times J(E) \rm \: [arb. units]$')
    plt.legend(title = r'Time values$\rm \: [yr]$', loc = 'best', ncol = 3)
    plt.savefig(f"{FIGURES_DIR}/energy_spectrum_over_time_{sim}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/energy_spectrum_over_time_{sim}.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
def plot_energy_spectrum_stationary(sim):

    data = np.loadtxt(f"{SIMULARIONS_DIR}/candidates_{sim}.txt")

    D, E, X, Y, Z = data[:,0], data[:,1], data[:,2], data[:,3], data[:,4]
    R = np.sqrt(X**2 + Y**2 + Z**2)
    T = D / c

    E_bins = np.logspace(np.log10(min(E)), np.log10(max(E)), num = 50) 
    t_i = get_time_intervals(sim)

    for i in range(118, len(t_i) - 1, 10):
        mask = (T <= t_i[i+1]) # & (R > 0) & (R < get_position_mask(sim))
        candidates = data[mask,:]
        E_center, spectrum, sigma = compute_energy_spectrum(candidates, E_bins, True, True)
        plt.errorbar(E_center, spectrum * E_center**2, sigma * E_center**2, ls = 'None', marker = '.', label = f'{float(t_i[i + 1] * s_to_yr):.1e}')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$\rm Energy \: [TeV]$')
    plt.ylabel(r'$E^2 \times J(E) \rm \: [arb. units]$')
    plt.legend(title = r'Time values$\rm \: [yr]$', loc = 'best', ncol = 3)
    plt.savefig(f"{FIGURES_DIR}/energy_spectrum_stationary_{sim}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/energy_spectrum_stationary_{sim}.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_energy_spectrum_over_time('3D_sph_Dconst')
    plot_energy_spectrum_stationary('3D_sph_Dconst')

# ----------------------------------------------------------------------------------------------------