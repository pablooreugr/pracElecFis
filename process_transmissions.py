import os
import re
import numpy as np
from scipy.integrate import simpson
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import csv

# ==========================================
# 1. PARAMETERS & FUNCTIONS
# ==========================================
T_temp = 300.0         # Temperature in Kelvin
mu_S   = 0.35          # Source Fermi level (eV)
mu_D   = 0.15          # Drain Fermi level (eV)
kB     = 8.61733326e-5 # Boltzmann constant (eV/K)

def fermi_dirac(E, mu, T):
    if T == 0:
        return np.heaviside(mu - E, 0.5)
    exponent = np.clip((E - mu) / (kB * T), -100, 100)
    return 1.0 / (1.0 + np.exp(exponent))

def parse_mass_str(mass_str):
    if mass_str == "11":
        return 1.1
    try:
        if len(mass_str) == 3:
            return float(mass_str) / 100.0
        return float(mass_str)
    except ValueError:
        return 0.0

data_dir = 'data/MOSFET'
file_pattern = re.compile(r'(\d+(?:_\d+)?)\s*nm\.txt$', re.IGNORECASE)
folder_pattern = re.compile(r'^(\d+)-06eV$')

results = [] # list of dicts: {'mass': ..., 'mass_str': ..., 'folder': ..., 'L': [], 'T_avg': [], 'Lc': ..., 'Lc_err': ...}

# ==========================================
# 2. PROCESS ALL DATA
# ==========================================
if not os.path.exists(data_dir):
    print(f"Error: Directory {data_dir} not found.")
    exit(1)

for folder in os.listdir(data_dir):
    folder_path = os.path.join(data_dir, folder)
    if not os.path.isdir(folder_path):
        continue
    
    match = folder_pattern.match(folder)
    if not match:
        continue
        
    mass_str = match.group(1)
    mass = parse_mass_str(mass_str)
    
    files = []
    for filename in os.listdir(folder_path):
        fmatch = file_pattern.search(filename)
        if fmatch:
            length_str = fmatch.group(1).replace('_', '.')
            length = float(length_str)
            files.append((length, os.path.join(folder_path, filename)))
            
    files.sort()
    
    lengths = []
    T_avgs = []
    
    for L, filepath in files:
        try:
            data = np.genfromtxt(filepath, delimiter=',', skip_header=4)
            transmission = data[:, 0]
            energy_eV    = data[:, 1]
            
            f_S = fermi_dirac(energy_eV, mu_S, T_temp)
            f_D = fermi_dirac(energy_eV, mu_D, T_temp)
            weight_function = f_S - f_D
            
            numerator = simpson(transmission * weight_function, x=energy_eV)
            denominator = simpson(weight_function, x=energy_eV)
            
            T_avg = numerator / denominator
            lengths.append(L)
            T_avgs.append(T_avg)
        except Exception as e:
            print(f"Error processing {filepath}: {e}")
            
    if len(lengths) > 2:
        results.append({
            'mass': mass,
            'mass_str': mass_str,
            'folder': folder,
            'L': np.array(lengths),
            'T_avg': np.array(T_avgs)
        })

# Sort results by mass for nice plotting
results.sort(key=lambda x: x['mass'])

# ==========================================
# 3. CURVE FITTING & CRITICAL LENGTH
# ==========================================
target_T = 1e-4
target_y = np.log(target_T)

for res in results:
    L = res['L']
    T = res['T_avg']
    
    # We want to fit ln(T) = a * L + b
    # Filter out the numerical noise floor at ~3e-6 which flattens the curves
    valid = T > 5e-6
    if np.sum(valid) < 2:
        valid = T > min(T) # Fallback to all non-minimum points
        
    L_fit = L[valid]
    y_fit = np.log(T[valid])
    
    if len(L_fit) == 2:
        # Cannot compute meaningful covariance for 2 points
        coeffs = np.polyfit(L_fit, y_fit, 1)
        a, b = coeffs[0], coeffs[1]
        var_a = var_b = cov_ab = 0.0
    else:
        coeffs, cov = np.polyfit(L_fit, y_fit, 1, cov=True)
        a, b = coeffs[0], coeffs[1]
        var_a = cov[0, 0]
        var_b = cov[1, 1]
        cov_ab = cov[0, 1]
    
    # Lc = (target_y - b) / a
    Lc = (target_y - b) / a
    
    # Error propagation
    if var_a > 0 and var_b > 0:
        dLc_db = -1.0 / a
        dLc_da = -(target_y - b) / (a**2)
        var_Lc = (dLc_db**2) * var_b + (dLc_da**2) * var_a + 2 * dLc_da * dLc_db * cov_ab
        Lc_err = np.sqrt(max(0, var_Lc))
    else:
        Lc_err = 0.0
    
    res['a'] = a
    res['b'] = b
    res['Lc'] = Lc
    res['Lc_err'] = Lc_err

# ==========================================
# 4. PLOTS & CSV EXPORT
# ==========================================
print("Writing critical_distances.csv ...")
with open('critical_distances.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(['Effective Mass (m0)', 'Folder', 'Critical Distance Lc (nm)', 'Error dLc (nm)', 'Fit Slope (a)', 'Fit Intercept (b)'])
    for res in results:
        writer.writerow([res['mass'], res['folder'], res['Lc'], res['Lc_err'], res['a'], res['b']])

# Plot 1: Combined Transmission vs Length
plt.figure(figsize=(10, 6))
cmap = plt.get_cmap('viridis')
colors = [cmap(i) for i in np.linspace(0, 1, len(results))]

for i, res in enumerate(results):
    L = res['L']
    T = res['T_avg']
    plt.plot(L, T, marker='o', linestyle='', color=colors[i], markersize=5, label=f"m* = {res['mass']}")
    
    # Plot the fit line as well over a slightly wider range
    L_fit = np.linspace(min(L), max(max(L), res['Lc']*1.1), 100)
    T_fit = np.exp(res['a'] * L_fit + res['b'])
    plt.plot(L_fit, T_fit, linestyle='-', color=colors[i], alpha=0.5)

plt.axhline(y=target_T, color='red', linestyle='--', label=f'Threshold (10^-4)')
plt.title('Average Transmission vs Transistor Length for various Effective Masses')
plt.xlabel('Length (nm)')
plt.ylabel('Average Transmission (T_avg)')
plt.yscale('log')
plt.grid(True, linestyle='--', alpha=0.6)
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()
plt.savefig('combined_transmission.png', dpi=300)
plt.close()

# Plot 2: Critical Distance vs Effective Mass
masses = [res['mass'] for res in results]
Lcs = [res['Lc'] for res in results]
Lc_errs = [res['Lc_err'] for res in results]

plt.figure(figsize=(8, 5))
plt.errorbar(masses, Lcs, yerr=Lc_errs, fmt='o-', color='navy', capsize=5, capthick=1.5, markerfacecolor='red')
plt.title('Critical Transistor Length vs Effective Mass')
plt.xlabel('Effective Mass (m_0)')
plt.ylabel('Critical Distance L_c (nm)')
plt.grid(True, linestyle='--', alpha=0.6)
plt.tight_layout()
plt.savefig('critical_distance_vs_mass.png', dpi=300)
plt.close()

print("Done! Check combined_transmission.png, critical_distance_vs_mass.png, and critical_distances.csv")
