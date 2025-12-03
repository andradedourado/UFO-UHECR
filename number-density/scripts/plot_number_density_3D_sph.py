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

m_to_pc = 3.2407793e-17 
s_to_yr = 3.1688e-8
yr_to_s = 31536000

c = 299792458 # m / s

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
def compute_number_density(candidates, bins, is_weighted, is_stationary):

    X, Y, Z = candidates[:,2], candidates[:,3], candidates[:,4] 
    R = np.sqrt(X**2 + Y**2 + Z**2)

    if is_stationary == False:
        w_T = np.ones(len(candidates))
    elif is_stationary == True:
        w_T = candidates[:,0] / c

    if is_weighted:
        weights = candidates[:,5]
        hist_counts_w, R_edges = np.histogram(R, bins = bins, weights = weights * w_T)
        hist_counts, _ = np.histogram(R, bins = bins, weights = w_T)
        num_density = hist_counts_w / (R_edges[1:] - R_edges[:-1])
        sigma = num_density / np.sqrt(hist_counts)
    
    else:
        hist_counts, R_edges = np.histogram(R, bins = bins, weights = w_T)
        R_width = R_edges[1:] - R_edges[:-1]
        num_density = hist_counts / R_width
        sigma = np.sqrt(hist_counts) / R_width

    R_center = R_edges[:-1] + 0.5 * (R_edges[1:] - R_edges[:-1])

    return R_center, num_density, sigma

# ----------------------------------------------------------------------------------------------------
def plot_number_density_over_time(sim):

    data = np.loadtxt(f"{SIMULARIONS_DIR}/candidates_{sim}.txt")

    D, X, Y, Z = data[:,0], data[:,2], data[:,3], data[:,4]
    R = np.sqrt(X**2 + Y**2 + Z**2)
    T = D / c

    num_density_bins = np.linspace(min(R), max(R), num = 100) 
    t_i = get_time_intervals(sim)

    for i in range(118, len(t_i) - 1, 10):
        mask = (T > t_i[i]) & (T < t_i[i+1])
        candidates = data[mask,:]
        R_center, num_density, sigma = compute_number_density(candidates, num_density_bins, True, False)
        plt.errorbar(R_center * m_to_pc, num_density, sigma, ls = 'None', marker = '.', label = f'{float(t_i[i + 1] * s_to_yr):.1e}')

    plt.axvline(x = 2.39e16 * m_to_pc, color = 'gray', linestyle = ':')
    plt.yscale('log')
    plt.xlabel(r'$X \rm \: [pc]$')
    plt.ylabel(r'Number density$\rm \: [arb. units]$')
    plt.legend(title = r'Time values$\rm \: [yr]$', loc = 'best')
    plt.savefig(f"{FIGURES_DIR}/num_density_over_time_{sim}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/num_density_over_time_{sim}.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
def plot_number_density_stationary(sim):

    data = np.loadtxt(f"{SIMULARIONS_DIR}/candidates_{sim}.txt")
    
    D, X, Y, Z = data[:,0], data[:,2], data[:,3], data[:,4]
    R = np.sqrt(X**2 + Y**2 + Z**2)
    T = D / c

    num_density_bins = np.linspace(min(R), max(R), num = 100) 
    t_i = get_time_intervals(sim)

    for i in range(118, len(t_i) - 1, 10):
        mask = (T <= t_i[i+1])
        candidates = data[mask,:]
        R_center, num_density, sigma = compute_number_density(candidates, num_density_bins, True, True)
        plt.errorbar(R_center * m_to_pc, num_density, sigma, ls = 'None', marker = '.', label = f'{float(t_i[i + 1] * s_to_yr):.1e}')

    plt.axvline(x = 2.39e16 * m_to_pc, color = 'gray', linestyle = ':')
    plt.yscale('log')
    plt.xlabel(r'$X \rm \: [pc]$')
    plt.ylabel(r'Number density$\rm \: [arb. units]$')
    plt.legend(title = r'Time values$\rm \: [yr]$', loc = 'best')
    plt.savefig(f"{FIGURES_DIR}/num_density_stationary_{sim}.pdf", bbox_inches = 'tight')
    plt.savefig(f"{FIGURES_DIR}/num_density_stationary_{sim}.png", bbox_inches = 'tight', dpi = 300)
    plt.show()

# ----------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    plot_number_density_over_time('3D_sph_Dconst')
    plot_number_density_stationary('3D_sph_Dconst')

# ----------------------------------------------------------------------------------------------------