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

1. **Download the Dataset**
   - Environment for undertaking simulations
     The yml file corresponding to setting up the environments for undertaking simulations are present in the Env folder. Plumed2 is actiavted in conda environment.

   - Steered_Molecular_Dynamics Scripts
     The scripts for running the simulations are provided in this folder. 
     The command to submit jobs: sh Complete_simulation_setup_hs.sh "Index_Number" "PDB_ID" "cpu-number" "partition". 

   - Trajectory_validation
     The codes to validate the vtrajectory in terms of sigmoidal curve fitting, RMSD of protein and ligand, center-of-mass distance separation, are provided in this folder.

   - Dataset
     The file structure corresponding to each PDB ID, are generated with the corrsponding code. These will be displayed and made public from the India-Data website ([https://india-data.org/dataset-details/ef3a1c5b-6ff2-49f7-ae7a-a99f69003849%22]), as mentioned in the manuscript.
     
   - Extract the `.tar.gz` archives to access individual PLC data

2. **Training Models**
   - The dataset is suitable for various ML/DL approaches:
     - Graph Neural Networks
     - 3D Convolutional Networks
     - Equivariant Neural Networks
     - Attention-based models

3. **Benchmarking**
   - Use the binding affinity data to evaluate model performance
   - Compare model predictions across both PLAS (bound) and HES (unbound) conformations

### Citation

If you use this dataset in your research, please cite:
*[Citation information to be provided by dataset creators after the publication of the dataset]*

### Contact

If you have any query regarding the dataset, you can reach out to Prathit Chatterjee (prathit.chatterjee@ihub-data.iiit.ac.in).

### License

*[License information to be provided by dataset creators: this should evidently include citing the dataset for all academic purposes, and to contact the authors for any commercial purpose]*

### Version

PLAS-HES-5k_

