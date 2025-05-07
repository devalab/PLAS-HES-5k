import re
import subprocess
import pandas as pd
import mdtraj as md
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

BASH = "/bin/bash" 

def extract(data, start_str, end_str):
    pattern = re.compile(f'{re.escape(start_str)}(.*?){re.escape(end_str)}', re.DOTALL)
    match = pattern.search(data)
    if match:
        return match.group(1).strip()
    else:
        return None

with open('plumed.dat', 'r') as file:
    dat_content = file.read()

p_r = extract(dat_content, 'c1: CENTER ATOMS=', 'MASS')
l_r = extract(dat_content, 'c2: CENTER ATOMS=', 'MASS')

with open('cpptraj_dist.in', 'w') as cpptraj_file:
    cpptraj_file.write('parm complex_solvated.prmtop\n')
    cpptraj_file.write('trajin plum.dcd\n')
    cpptraj_file.write(f'distance @{p_r}  @{l_r} out P-L.dat\n')
    cpptraj_file.write('run\n')

command1 = f"""
cpptraj  -i cpptraj_dist.in
"""
subprocess.run(command1, shell=True, check=True, executable=BASH)


df = pd.read_csv('P-L.dat', delim_whitespace=True)

plt.plot(df['#Frame'], df['Dis_00001'])
plt.xlabel('Frame')
plt.ylabel('Distance binding_pocket-LIG(in Å)')

initial_dist = df['Dis_00001'].iloc[0]
final_dist = df['Dis_00001'].iloc[-1]

plt.text(0.5, 0.9, f'Initial Dist: {initial_dist:.2f} Å', transform=plt.gca().transAxes, color='blue', fontsize=10, verticalalignment='top')
plt.text(0.5, 0.85, f'Final Dist: {final_dist:.2f} Å', transform=plt.gca().transAxes, color='red', fontsize=10, verticalalignment='top')

plt.show()

# Sigmoid function
def sigmoid(x, a, b, c, d):
    return a / (1 + np.exp(-c * (x - d))) + b

def fit_sigmoid(x, y):
    initial_guess = [np.max(y) - np.min(y), np.min(y), 0.01, (np.max(x) + np.min(x)) / 2]
    bounds = ([0, -np.inf, 0, -np.inf], [np.inf, np.inf, np.inf, np.inf])
    popt, _ = curve_fit(sigmoid, x, y, p0=initial_guess, bounds=bounds, maxfev=3000)
    return sigmoid(x, *popt)

file_path_1 = 'Interaction_energy.csv'
csv1 = pd.read_csv(file_path_1, skiprows=612, header=None)
last_column = csv1.iloc[:, -1]
last_column.to_csv("last_column.csv", index=False, header=False)
file_path = 'last_column.csv'
data = pd.read_csv(file_path, header=None)

data = pd.to_numeric(data[0], errors='coerce').dropna()
original_indices = data.index
filtered_data = data[data < 0].reset_index(drop=True)
filtered_indices = pd.Series(original_indices[data < 0].values)

x = np.arange(1, len(filtered_data) + 1)
y = filtered_data.values

# Fit sigmoid curve
fitted_curve = fit_sigmoid(x, y)

threshold = 0.75

# Normalize the fitted curve
scaled_curve = (fitted_curve - np.min(fitted_curve)) / (np.max(fitted_curve) - np.min(fitted_curve))

# Loop to adjust the threshold
while True:
    frame_at_threshold = np.where(scaled_curve >= threshold)[0][0]
    
    # Divide data into lower and upper halves
    lower_half_mask = x < frame_at_threshold
    upper_half_mask = x >= frame_at_threshold

    if len(x[upper_half_mask]) >= 30:
        break  
    
    threshold -= 0.05  

lower_x, lower_y = x[lower_half_mask], y[lower_half_mask]
upper_x, upper_y = x[upper_half_mask], y[upper_half_mask]

fitted_lower_curve = scaled_curve[lower_half_mask]
fitted_upper_curve = scaled_curve[upper_half_mask]

# Find closest 10 points in lower half
lower_diff = np.abs(lower_y - fitted_curve[lower_half_mask])
n = len(lower_diff)
bin_size = n // 10
remainder = n % 10
bin_sizes = [bin_size + 1 if i < remainder else bin_size for i in range(10)]
lower_closest_indices = [np.argmin(lower_diff[start : start + size]) + start for start, size in zip(np.cumsum([0] + bin_sizes[:-1]), bin_sizes)]
marked_lower_x = lower_x[lower_closest_indices]
marked_lower_y = lower_y[lower_closest_indices]

# Find closest 30 points in upper half
upper_diff = np.abs(upper_y - fitted_curve[upper_half_mask])
n = len(upper_diff)
bin_size = n // 30
remainder = n % 30
bin_sizes = [bin_size + 1 if i < remainder else bin_size for i in range(30)]
upper_closest_indices = [np.argmin(upper_diff[start : start + size]) + start for start, size in zip(np.cumsum([0] + bin_sizes[:-1]), bin_sizes)]
marked_upper_x = upper_x[upper_closest_indices]
marked_upper_y = upper_y[upper_closest_indices]

# Retrieve original frame numbers for the selected points
marked_lower_original_indices = filtered_indices.iloc[lower_closest_indices]
marked_upper_original_indices = filtered_indices.iloc[upper_closest_indices]
marked_x = np.concatenate([marked_lower_x, marked_upper_x])
marked_y = np.concatenate([marked_lower_y, marked_upper_y])
marked_original_indices = np.concatenate([marked_lower_original_indices, marked_upper_original_indices])

plt.figure(figsize=(10, 8))
plt.scatter(x, y, color='dodgerblue', s=10, label='All data points')
plt.plot(x, fitted_curve, color='blue', linewidth=2, label='Sigmoid Curve fit')
plt.axvline(frame_at_threshold, color='magenta', linestyle='--', linewidth=1, label=f'x_{threshold:.2f}')
plt.scatter(marked_x, marked_y, color='red', s=50, label='Selected frames')
plt.xlabel("frames", fontsize=12)
plt.ylabel("BFE(kcal/mol)", fontsize=12)
plt.title("Sigmoid Curve Fit with Selected Frames", fontsize=14)
plt.legend()
plt.show()


def calculate_rmsd(traj, ref):
    return md.rmsd(traj, ref)

def plot_rmsd(time, rmsd_values, label, xlabel, ylabel):
    plt.plot(time, rmsd_values, label=label)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.legend()
    
traj = md.load('plum.dcd', top='sim1_last.pdb')  
ref = md.load('sim1_last.pdb')
protein_sel = traj.topology.select('protein')
protein_traj = traj.atom_slice(protein_sel)
ligand_sel = traj.topology.select('resname LIG')
ligand_traj = traj.atom_slice(ligand_sel)
system_sel = traj.topology.select('resid 1 to 126')
system_traj = traj.atom_slice(system_sel)

rmsd_protein = md.rmsd(protein_traj, protein_traj, 0)
rmsd_ligand_ligand = md.rmsd(ligand_traj, ligand_traj, 0)
#rmsd_system = md.rmsd(system_traj, system_traj, 0)

plot_rmsd(traj.time, rmsd_protein, 'Protein RMSD', 'Time (ps)', 'RMSD (nm)')
plot_rmsd(traj.time, rmsd_ligand_ligand, 'Ligand-Ligand RMSD', 'Time (ps)', 'RMSD (nm)')
#plot_rmsd(traj.time, rmsd_system, 'System RMSD', 'Time (ps)', 'RMSD (nm)')

plt.show()
