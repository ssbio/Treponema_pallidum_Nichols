import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import font_manager

# Read data from kcat.txt
with open('Kcat_TP.txt', 'r') as file:
    kcat_data = np.array([float(line.strip()) for line in file.readlines() if float(line.strip()) <= 100])

# Create a DataFrame for the data
df_kcat = pd.DataFrame({'Kcat': kcat_data})

# Plotting function with formatting
def plot_density(data, variable, color, title, xlabel, ylabel):
    plt.figure(figsize=(16, 10))
    sns.set_style("white")
    sns.kdeplot(data=data, x=variable, color=color, fill=True, linewidth=2)
    plt.xlabel(xlabel, fontsize=32, fontweight='bold')
    plt.ylabel(ylabel, fontsize=32, fontweight='bold')
    plt.xticks(fontsize=22)
    plt.yticks(fontsize=22)
    plt.title(title)
    plt.tight_layout()
    plt.show()

# Plot original density plot for Kcat
plot_density(df_kcat, 'Kcat', 'skyblue', '', 'Kcat (1/s)', 'Density')

# Perform Monte Carlo simulation for Kcat
num_samples_kcat = 80
num_simulations_kcat = 100
monte_carlo_results_kcat = []

for _ in range(num_simulations_kcat):
    kcat_samples = np.random.choice(kcat_data, num_samples_kcat, replace=True)
    monte_carlo_results_kcat.append(kcat_samples)

# Write Kcat Monte Carlo simulation results to a text file
with open('Kcat_MC_TP.txt', 'w') as file:
    for i, result in enumerate(monte_carlo_results_kcat):
        file.write(f"Simulation {i+1}:\n")
        for kcat in result:
            file.write(f"{kcat}\n")
        file.write("\n")

# Define the font properties
legend_font_props = font_manager.FontProperties(style='normal', size=22, weight='bold')

# Plot the original and Monte Carlo simulation results for Kcat
plt.figure(figsize=(16, 10))
sns.set_style("white")
sns.kdeplot(data=kcat_data, color='blue', fill=True, linewidth=2, label='Original Kcat Distribution')
for i in range(num_simulations_kcat):
    sns.kdeplot(data=monte_carlo_results_kcat[i], color='orange', alpha=0.3, linewidth=1)
plt.xlabel('Kcat (1/s)', fontsize=32, fontweight='bold')
plt.ylabel('Density', fontsize=32, fontweight='bold')
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.legend(['Original Kcat Distribution', 'Monte Carlo Simulation'], loc='upper right', prop=legend_font_props)
plt.tight_layout()
plt.savefig('Kcat_distribution_TP.png')
plt.show()
