# HES5K
A dataset for the decoys of protein-ligand binding, useful for training machine learning models.
 
# Environment for undertaking simulations
The yml file corresponding to setting up the environments for undertaking simulations are present in the Env folder. Plumed2 is actiavted in conda environment.

# Steered_Molecular_Dynamics Scripts
The scripts for running the simulations are provided in this folder. 
The command to submit jobs: sh Complete_simulation_setup_hs.sh "Index_Number" "PDB_ID" "cpu-number" "partition". 

# Trajectory_validation
The codes to validate the vtrajectory in terms of sigmoidal curve fitting, RMSD of protein and ligand, center-of-mass distance separation, are provided in this folder.

# File_Structure
The file structure corresponding to each PDB ID, are generated with the corrsponding code. These will be displayed and ade public from the India-Data website, as mentioned in the manuscript.
