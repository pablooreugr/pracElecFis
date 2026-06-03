import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import os

# Load data
csv_path = 'critical_distances.csv'
df = pd.read_csv(csv_path)
df.columns = [col.strip() for col in df.columns]

x = df['Effective Mass (m0)'].values
y = df['Critical Distance Lc (nm)'].values
y_err = df['Error dLc (nm)'].values

# Define models
def power_law(x, a, b):
    return a * np.power(x, b)

def exponential(x, a, b):
    return a * np.exp(b * x)

def inverse_sqrt(x, a, b):
    return a / np.sqrt(x) + b

def inverse(x, a, b):
    return a / x + b

models = {
    'Power Law ($y = a \cdot x^b$)': power_law,
    'Exponential ($y = a \cdot e^{bx}$)': exponential,
    'Inverse Sqrt ($y = a/\sqrt{x} + b$)': inverse_sqrt,
    'Inverse ($y = a/x + b$)': inverse
}

results = {}
x_fit = np.linspace(min(x) * 0.9, max(x) * 1.1, 100)

plt.figure(figsize=(10, 7))
plt.errorbar(x, y, yerr=y_err, fmt='o', label='Data', color='black', capsize=4)

colors = ['red', 'blue', 'green', 'purple']
for i, (name, func) in enumerate(models.items()):
    try:
        # Initial guesses
        if func == power_law:
            p0 = [1, -0.5]
        elif func == exponential:
            p0 = [1, -1]
        else:
            p0 = [1, 0]
            
        popt, pcov = curve_fit(func, x, y, sigma=y_err, absolute_sigma=True, p0=p0, maxfev=10000)
        
        # Calculate R-squared
        residuals = y - func(x, *popt)
        ss_res = np.sum(residuals**2)
        ss_tot = np.sum((y - np.mean(y))**2)
        r_squared = 1 - (ss_res / ss_tot)
        
        results[name] = {
            'params': popt,
            'r_squared': r_squared,
            'errors': np.sqrt(np.diag(pcov))
        }
        
        y_fit = func(x_fit, *popt)
        label = f"{name}\n$R^2 = {r_squared:.4f}$"
        plt.plot(x_fit, y_fit, color=colors[i], label=label, linewidth=2)
        
        print(f"--- {name} ---")
        if func == power_law:
            print(f"a = {popt[0]:.4f} ± {results[name]['errors'][0]:.4f}")
            print(f"b = {popt[1]:.4f} ± {results[name]['errors'][1]:.4f}")
        elif func == exponential:
            print(f"a = {popt[0]:.4f} ± {results[name]['errors'][0]:.4f}")
            print(f"b = {popt[1]:.4f} ± {results[name]['errors'][1]:.4f}")
        else:
            print(f"a = {popt[0]:.4f} ± {results[name]['errors'][0]:.4f}")
            print(f"b = {popt[1]:.4f} ± {results[name]['errors'][1]:.4f}")
        print(f"R^2 = {r_squared:.4f}\n")
        
    except Exception as e:
        print(f"Could not fit {name}: {e}")

# Save results to CSV
csv_out = 'regression_results.csv'
with open(csv_out, 'w') as f:
    f.write('Model,Parameter a,Error a,Parameter b,Error b,R_squared\n')
    for name, res in results.items():
        a = res['params'][0]
        err_a = res['errors'][0]
        b = res['params'][1]
        err_b = res['errors'][1]
        r2 = res['r_squared']
        # Clean name for CSV
        clean_name = name.split(' (')[0]
        f.write(f'{clean_name},{a},{err_a},{b},{err_b},{r2}\n')
print(f"Regression results saved to {csv_out}")

plt.title('Critical Distance ($L_c$) vs Effective Mass ($m^*$)\nRegression Analysis', fontsize=14)
plt.xlabel('Effective Mass ($m_0$)', fontsize=12)
plt.ylabel('Critical Distance $L_c$ (nm)', fontsize=12)
plt.legend(fontsize=10)
plt.grid(True, linestyle='--', alpha=0.7)
plt.tight_layout()
plt.savefig('regression_analysis.png', dpi=300)
print("Plot saved to regression_analysis.png")
