import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. SETUP
# ==========================================
# Specify the data file to use
# Make sure this file is in the same directory as the script,
# or provide the full path.
data_filename = 'L_vs_Tavg_098-06eV.txt'

# Target average transmission value
target_T_avg = 1e-4

# ==========================================
# 2. LOAD AND PREPARE DATA
# ==========================================
try:
    # Load data from the text file, skipping the header row
    data = np.loadtxt(data_filename, skiprows=1)
    
    # Extract length (L) and average transmission (T_avg)
    L = data[:, 0]
    T_avg = data[:, 1]
    
    # For interpolation, it's better to work in log space for T_avg
    # since it spans many orders of magnitude and the relationship
    # with L is approximately log-linear.
    log_T_avg = np.log10(T_avg)
    
    print(f"Successfully loaded data from '{data_filename}'")
    
except FileNotFoundError:
    print(f"Error: The file '{data_filename}' was not found.")
    exit()
except Exception as e:
    print(f"An error occurred while reading the file: {e}")
    exit()

# ==========================================
# 3. INTERPOLATION
# ==========================================
# The goal is to find the length L where T_avg = 1e-4.
# This is equivalent to finding L where log10(T_avg) = -4.

# We will use numpy's interpolation function `np.interp`.
# The function signature is: np.interp(x_new, x_values, y_values)
# It finds the y_new value for a given x_new.

# In our case, we want to find L (our 'y') for a given log_T_avg (our 'x').
# IMPORTANT: `np.interp` requires the x_values to be monotonically increasing.
# Our log_T_avg is decreasing as L increases. So, we must reverse
# both arrays before passing them to the function.
log_T_avg_increasing = log_T_avg[::-1]
L_for_increasing_log_T = L[::-1]

# The target value for interpolation
target_log_T = np.log10(target_T_avg) # This is -4.0

# Perform the linear interpolation
interpolated_L = np.interp(target_log_T, log_T_avg_increasing, L_for_increasing_log_T)

# ==========================================
# 4. RESULTS & VISUALIZATION
# ==========================================
print("-" * 40)
print(f"Target T_avg: {target_T_avg:.1e}")
print(f"Interpolated Length (L): {interpolated_L:.4f} nm")
print("-" * 40)

plt.figure(figsize=(10, 6))
plt.plot(L, T_avg, 'o-', label='Original Data', color='indigo')
plt.axhline(y=target_T_avg, color='red', linestyle='--', label=f'Target T_avg = {target_T_avg:.1e}')
plt.plot(interpolated_L, target_T_avg, 'X', color='lime', markersize=10, mew=2, markeredgecolor='black', label=f'Interpolated Point\nL = {interpolated_L:.4f} nm')
plt.yscale('log')
plt.grid(True, which="both", ls="--", alpha=0.6)
plt.xlabel('Transistor Length (nm)')
plt.ylabel('Weighted Average Transmission ($T_{avg}$)')
plt.title('Interpolation to find Length for a Target Transmission')
plt.legend()
plt.tight_layout()
plt.show()
