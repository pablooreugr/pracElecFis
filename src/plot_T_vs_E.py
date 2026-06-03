import numpy as np
import matplotlib.pyplot as plt
import os

def plot_transmission(folder, mass_label, length_file="5nm.txt", out_file="out.png"):
    filepath = os.path.join("../data/MOSFET", folder, length_file)
    if not os.path.exists(filepath):
        print(f"File not found: {filepath}")
        return
        
    # Read data (skip first 4 header lines)
    data = np.genfromtxt(filepath, delimiter=',', skip_header=4)
    transmission = data[:, 0]
    energy_eV = data[:, 1]
    
    plt.figure(figsize=(8, 5))
    plt.plot(energy_eV, transmission, color='blue', linewidth=2)
    plt.title(f'Transmisión vs Energía\n(Masa Efectiva = {mass_label}, L = 5 nm)')
    plt.xlabel('Energía (eV)')
    plt.ylabel('Transmisión ($T$)')
    plt.yscale('log')
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.tight_layout()
    plt.savefig(f"../graficas/{out_file}", dpi=300)
    plt.close()
    print(f"Saved {out_file}")

# Extremes
plot_transmission("006-06eV", "0.06 $m_0$", "5nm.txt", "T_vs_E_m006.png")
plot_transmission("11-06eV", "1.10 $m_0$", "5nm.txt", "T_vs_E_m11.png")
