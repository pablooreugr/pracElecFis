import numpy as np
from scipy.integrate import simpson
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import os
import re

# ==========================================
# 1. ADJUSTABLE PARAMETERS
# ==========================================
T_temp = 300.0         # Temperature in Kelvin
mu_S   = 0.35          # Source Fermi level (eV)
mu_D   = 0.15          # Drain Fermi level (eV)
kB     = 8.61733326e-5 # Boltzmann constant (eV/K)

def fermi_dirac(E, mu, T):
    """Calculates the Fermi-Dirac probability."""
    if T == 0:
        return np.heaviside(mu - E, 0.5)
    exponent = np.clip((E - mu) / (kB * T), -100, 100)
    return 1.0 / (1.0 + np.exp(exponent))

# ==========================================
# 2. DATA FILE DISCOVERY
# ==========================================
data_dir = 'data/MOSFET/110-06eV'
file_pattern = re.compile(r'(\d+(?:_\d+)?)\s*nm\.txt$', re.IGNORECASE)

lengths = []
weighted_averages = []

# ==========================================
# 3. LOOP THROUGH FILES & CALCULATE
# ==========================================
print("Starting calculations...")
print("-" * 40)

# Find files like 10nm.txt and extract the length from the file name
files = []
for filename in os.listdir(data_dir):
    match = file_pattern.search(filename)
    if not match:
        continue
    length_str = match.group(1).replace('_', '.')
    length = float(length_str)
    filepath = os.path.join(data_dir, filename)
    files.append((length, filepath))

for length, filename in sorted(files):
    if not os.path.exists(filename):
        print(f"Warning: '{filename}' not found in the directory. Skipping...")
        continue
        
    try:
        # Load the data, skipping the 4 header lines
        data = np.genfromtxt(filename, delimiter=',', skip_header=4)
        transmission = data[:, 0]
        energy_eV    = data[:, 1]
        
        # Calculate the weight function (Fermi window)
        f_S = fermi_dirac(energy_eV, mu_S, T_temp)
        f_D = fermi_dirac(energy_eV, mu_D, T_temp)
        weight_function = f_S - f_D
        
        # Calculate weighted average integrals
        numerator = simpson(transmission * weight_function, x=energy_eV)
        denominator = simpson(weight_function, x=energy_eV)
        
        # Calculate the final weighted average for this length
        T_weighted_avg = numerator / denominator
        
        # Store results for plotting
        lengths.append(length)
        weighted_averages.append(T_weighted_avg)
        
        print(f"L = {length} nm | File: {filename} | T_avg = {T_weighted_avg:.4e}")
        
    except Exception as e:
        print(f"Error processing '{filename}': {e}")

print("-" * 40)

# ==========================================
# 4. PLOTTING THE RESULTS
# ==========================================
if lengths and weighted_averages:
    # Define output directory and folder name from data path
    output_dir = os.path.dirname(os.path.abspath(__file__))
    folder_name = os.path.basename(os.path.normpath(data_dir))

    # --- Save data to a .txt file ---
    output_data_filename = f"L_vs_Tavg_{folder_name}.txt"
    output_data_path = os.path.join(output_dir, output_data_filename)

    # Combine data and save
    output_data = np.array([lengths, weighted_averages]).T
    np.savetxt(output_data_path, output_data, header="Length (nm)\tT_avg", delimiter='\t', fmt='%.18e')

    print(f"Saved data to: {output_data_path}")
    print("-" * 40)

    plt.figure(figsize=(8, 5))
    
    # Plot Lengths vs Weighted Averages
    plt.plot(lengths, weighted_averages, linestyle='dashed', color='indigo', marker='o', linewidth=2, markersize=4,)
    
    # Añadir línea discontinua horizontal en y=1e-4
    plt.axhline(y=1e-4, color='red', linestyle='--')
    
    # Formatting the plot
    plt.title('Weighted Average Transmission vs. Transistor Length')
    plt.xlabel('Transistor Length (nm)')
    plt.ylabel('Weighted Average Transmission ($T_{avg}$)')
    
    # Ensure only the actual provided lengths appear on the X-axis
    plt.xticks(lengths) 
    
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.yscale('log')
    ax = plt.gca()
    ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0))
    ax.yaxis.set_major_formatter(ticker.LogFormatterMathtext())
    plt.tight_layout()

    output_filename = f"transmission_{folder_name}.png"
    output_path = os.path.join(output_dir, output_filename)
    plt.savefig(output_path, dpi=300)
    print(f"Saved plot as: {output_path}")
    plt.show()
else:
    print("No valid data was processed, so no plot will be generated. Please check your filenames!")