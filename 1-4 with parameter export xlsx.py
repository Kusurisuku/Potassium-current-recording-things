# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 14:01:19 2025

@author: Junting
"""

import numpy as np
import pandas as pd
import logging
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, differential_evolution
from scipy.stats import f

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# Load Data
file_path = 
Time_name    = ''
Current_name = ''
pF_name      = ''

def load_data(file_path):
    df = pd.read_csv(file_path).dropna(subset=[Time_name, Current_name])
    time = df[Time_name].values
    current_raw = df[Current_name].values
    pF = df[pF_name].values[0]
    return time, current_raw / pF, pF  # Normalize current by pF

time, current, pF = load_data(file_path)
upper_bound_current = max(current) + 1


# Define exponential models
def single_exponential(t, A1r, tau1, As):
    return A1r * np.exp(-t / tau1) + As

def double_exponential(t, A1r, tau1, A2r, tau2, As):
    return A1r * np.exp(-t / tau1) + A2r * np.exp(-t / tau2) + As

def triple_exponential(t, A1r, tau1, A2r, tau2, A3r, tau3, As):
    return A1r * np.exp(-t / tau1) + A2r * np.exp(-t / tau2) + A3r * np.exp(-t / tau3) + As

def quadruple_exponential(t, A1r, tau1, A2r, tau2, A3r, tau3, A4r, tau4, As):
    return (A1r * np.exp(-t / tau1) + A2r * np.exp(-t / tau2) +
            A3r * np.exp(-t / tau3) + A4r * np.exp(-t / tau4) + As)

models = [single_exponential, double_exponential, triple_exponential, quadruple_exponential]
#models = [single_exponential, double_exponential, triple_exponential]
bounds_list = [
    [(0.001, upper_bound_current), (0.001, 5000), (0.001, upper_bound_current)],  # Single
    [(0.001, upper_bound_current), (0.001, 200), (0.001, upper_bound_current), (0.001, 4000), (0.001, upper_bound_current)],  # Double
    [(0.001, upper_bound_current), (0.001, 200), (0.001, upper_bound_current), (0.001, 1500), (0.001, upper_bound_current), (0.001, 7000), (0.001, upper_bound_current)],  # Triple
    [(0.001, upper_bound_current), (0.001, 200), (0.001, upper_bound_current), (0.001, 1500), (0.001, upper_bound_current), (0.001, 5000), (0.001, upper_bound_current), (0.001, 7000), (0.001, upper_bound_current)]  # Quadruple
]

# Fitting process
results = []
plt.figure(figsize=(10, 6))
for i, model in enumerate(models, start=1):
    logging.info(f"Fitting {i}-exponential model...")
    
    def objective_function(params):
        return np.sum((current - model(time, *params)) ** 2)
    
    result = differential_evolution(objective_function, bounds_list[i-1])
    optimized_params = result.x
    params_final, _ = curve_fit(model, time, current, p0=optimized_params, maxfev=20000)


#%%
    # Calculate metrics
    SSE = np.sum((current - model(time, *params_final)) ** 2)
    dof = len(time) - len(params_final)
    reduced_chi_square = SSE / dof
    R_value = np.corrcoef(current, model(time, *params_final))[0, 1]
    Ikv_total = max(current)
    
    #A_values = [params_final[2 * j] / sum(params_final[2 * k] for k in range(i)) * 100 for j in range(i)]
    A_values = [params_final[2 * j] for j in range(i)]
    taus = [params_final[2 * j + 1] for j in range(i)]
    As = params_final[-1]

    results.append({
        "Model": f"{i}-Exponential",
        "A1 (pA/pF)": A_values[0] if i > 0 else np.nan,
        "A2 (pA/pF)": A_values[1] if i > 1 else np.nan,
        "A3 (pA/pF)": A_values[2] if i > 2 else np.nan,
        "A4 (pA/pF)": A_values[3] if i > 3 else np.nan,
        "Tau1 (ms)": taus[0] if i > 0 else np.nan,
        "Tau2 (ms)": taus[1] if i > 1 else np.nan,
        "Tau3 (ms)": taus[2] if i > 2 else np.nan,
        "Tau4 (ms)": taus[3] if i > 3 else np.nan,
        "As (pA)": As,
        "SSE": SSE,
        "Reduced Chi-Square": reduced_chi_square,
        "R Value": R_value,
        "Degrees of Freedom": dof,
    })

    plt.plot(time, model(time, *params_final), label=f"{i}-Exponential Fit")

# Pairwise F-Test for model comparison
for i in range(len(results) - 1):
    SSE_i, SSE_ip1 = results[i]["SSE"], results[i + 1]["SSE"]
    k_i, k_ip1 = len(models[i].__code__.co_varnames) - 1, len(models[i + 1].__code__.co_varnames) - 1
    n = len(time)
    
    F_stat = ((SSE_i - SSE_ip1) / (k_ip1 - k_i)) / (SSE_ip1 / (n - k_ip1))
    p_value = 1 - f.cdf(F_stat, k_ip1 - k_i, n - k_ip1)
        
    results[i + 1]["F-stat"] = F_stat
    results[i + 1]["p-value"] = p_value
    results[i + 1]["Significant?"] = "Yes" if p_value < 0.05 else "No"
    results[i + 1]["Total IK current"] = Ikv_total

# Select best model
best_model = "1-Exponential"
for i in range(len(results) - 1):
    if results[i + 1]["Significant?"] == "Yes":
        best_model = results[i + 1]["Model"]

#%%
# Save results
pd.DataFrame(results).to_excel("Fitting_Results1.xlsx", index=False)

# Plot settings
plt.scatter(time, current, color='black', label='Data', s=10)
plt.xlabel("Time (ms)")
plt.ylabel("Normalized Current (pA/pF)")
plt.legend()
plt.title("Fitted Exponential Models")
plt.show()

logging.info(f"Best Model Selected: {best_model}")
print(f"Best Model: {best_model}")

