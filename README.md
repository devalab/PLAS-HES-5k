# PLAS-HES-5k Dataset

## Protein-Ligand Complex Trajectories for ML/DL in Data Driven Drug Discovery

### Overview

The **PLAS-HES-5k** dataset is a curated collection of 5,000 protein-ligand complexes (PLCs) designed specifically for machine learning applications in drug discovery. This comprehensive dataset includes both bound and unbound conformations, annotated with binding free energies derived from non-equilibrium molecular dynamics simulations and approximate free energy calculations. The structural and energetic diversity represented in PLAS-HES-5k makes it an ideal benchmark and training resource for predictive and generative ML/DL models aimed at improving drug-target binding predictions.

### Dataset Contents

The dataset contains:

- **5,000 protein-ligand complexes** in two conformational variants:
  - **PLAS**: Bound conformations starting from PDB database
  - **HES**: Synthetic unbound conformations derived from PLAS-20K, including both low and high energy states

- **For each complex**:
  - Atomic coordinates
  - Binding free energy data
  - Parameter files for simulation

### Structure and File Organization

Each protein-ligand complex in the dataset includes:

#### PLAS Data
- Bound conformations from PLAS-20K
- Complete atomic coordinates in PDB format
- Parameter files (.prmtop)
- Binding affinity calculations

#### HES Data
- Synthetic unbound conformations seeded from PLAS-20K
- Both low and high energy conformational states
- Complete atomic coordinates in PDB format
- Parameter files (.prmtop)
- Binding affinity calculations

### File Formats

- `.tar.gz` archives containing collections of PLCs
- `.txt` files listing which PLCs are present in each archive
- `.pdb` files containing structural information
- `.prmtop` files containing molecular dynamics parameters
- `.csv` files with binding affinity data

### Intended Use

The PLAS-HES-5k dataset is designed for:

1. **Training ML/DL models** for predicting protein-ligand binding affinities
2. **Evaluating generative models** in drug design
3. **Research on conformational sampling** and binding dynamics
4. **Benchmarking** novel ML/DL architectures for drug discovery applications

### Usage Instructions

1. **Environment Setup**
   - Environment configuration files for simulations are located in the `Env` folder
   - A YAML file is provided for setting up the conda environment with Plumed2 activated

2. **Steered Molecular Dynamics**
   - Simulation scripts are available in the `Steered_Molecular_Dynamics` folder
   - Submit jobs using the command:
     ```
     sh Complete_simulation_setup_hs.sh "Index_Number" "PDB_ID" "cpu-number" "partition"
     ```

3. **Trajectory Validation**
   - The `Trajectory_validation` folder contains scripts to validate trajectories through:
     - Sigmoidal curve fitting
     - RMSD analysis of protein and ligand
     - Center-of-mass distance separation measurements

4. **Dataset Access**
   - File structures are generated for each PDB ID
   - The complete dataset is publicly available on the India-Data website:
     [https://india-data.org/dataset-details/ef3a1c5b-6ff2-49f7-ae7a-a99f69003849](https://india-data.org/dataset-details/ef3a1c5b-6ff2-49f7-ae7a-a99f69003849)
   - Extract the `.tar.gz` archives to access individual PLC data

5. **Energy Component Analysis**
   - The `Distribution_Of_Energy_Components` folder contains scripts for:
     - Reproducing energy component distributions across all PLCs
     - Analyzing energy components for individual PLCs

6. **Training Machine Learning Models**
   - The dataset is suitable for various ML/DL approaches:
     - Graph Neural Networks
     - 3D Convolutional Networks
     - Equivariant Neural Networks
     - Attention-based models

7. **Benchmarking**
   - Use the binding affinity data to evaluate model performance
   - Compare model predictions across both PLAS (bound) and HES (unbound) conformations

### Citation

If you use this dataset in your research, please cite:
*[Citation information to be provided by dataset creators after the publication of the dataset]*

### Contact

If you have any query regarding the dataset, you can reach out to Prathit Chatterjee (prathit.chatterjee@ihub-data.iiit.ac.in).

## License

This work is licensed under a [Creative Commons Attribution-NonCommercial-NoDerivatives 4.0 International License](https://creativecommons.org/licenses/by-nc-nd/4.0/).

[![CC BY-NC-ND 4.0][cc-by-nc-nd-image]][cc-by-nc-nd]

[cc-by-nc-nd]: https://creativecommons.org/licenses/by-nc-nd/4.0/
[cc-by-nc-nd-image]: https://licensebuttons.net/l/by-nc-nd/4.0/88x31.png

### Version

1.0.0

