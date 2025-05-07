#!/bin/bash
#SBATCH -A account
#SBATCH --job-name=job_name
#SBATCH -p partition
#SBATCH --gres=gpu:1
#SBATCH -n ncores
##SBATCH --exclude=gnode[100,102]
#SBATCH --mem-per-cpu=3000
#SBATCH --time=4-00:00:00


module load u18/cuda/11.6
module add u18/cuda/11.6
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/cuda-11.6/lib64
eval "$(conda shell.bash hook)"
conda activate plumed2

script_dir="/path/to/smd_scripts"
input_dir="/scratch"
home_dir="/path/to_store/prepared_input"
store_dir="/path/to_store/outputs"

cp -r ${home_dir}/index.pdbid ${input_dir}/
cd ${input_dir}/index.pdbid/ 
python ${script_dir}/input.py index pdbid > input_output
mkdir ${input_dir}/index.pdbid/Results
cp ${input_dir}/index.pdbid/complex_solvated.prmtop ${input_dir}/index.pdbid/Results
cp ${input_dir}/index.pdbid/run0/sim1.dcd ${input_dir}/index.pdbid/Results
cp ${script_dir}/cpptraj_2.in ${input_dir}/index.pdbid/Results
cp ${script_dir}/cpptraj_200.in ${input_dir}/index.pdbid/Results
cp ${script_dir}/mmpbsa.in ${input_dir}/index.pdbid/Results
cp ${script_dir}/water_selection_no_skip.py ${input_dir}/index.pdbid/Results
cd ${input_dir}/index.pdbid/Results/
cpptraj -i cpptraj_2.in
ante-MMPBSA.py -p 2.complex_solvated.prmtop -c complex.prmtop -r protein.prmtop -l ligand.prmtop -s :0 -n :LIG
rsync -azp ${input_dir}/index.pdbid ada:$store_dir/
for NRUN in 0 1 2 3 4
do
  cd ${input_dir}/index.pdbid/run${NRUN}/
  python pull.py > simulation_output_${NRUN}
  result=$(grep -o "Simulation complete in [0-9.]* seconds" simulation_output_${NRUN})
  if [ -n "$result" ]; then
      echo -n "index pdbid run${NRUN} plumed_status True" >> ${input_dir}/index.pdbid/Result_status
  else
      echo -n "index pdbid run${NRUN} plumed_status False" >> ${input_dir}/index.pdbid/Result_status
  fi
  mkdir ${input_dir}/index.pdbid/Results/run${NRUN}
  cp ${input_dir}/index.pdbid/run${NRUN}/plum.dcd ${input_dir}/index.pdbid/Results/run${NRUN}
  cp ${input_dir}/index.pdbid/run${NRUN}/sim1_last.pdb ${input_dir}/index.pdbid/Results/run${NRUN}
  cd ${input_dir}/index.pdbid/Results/run${NRUN}/
  cpptraj -i ../cpptraj_200.in >> mmpbsa_output
  python  ../water_selection_no_skip.py -traj 200.dcd -top 200.complex_solvated.prmtop -nwater 2 -solvated_complex sim1_last.pdb >> mmpbsa_output
  mpirun -np ncores MMPBSA.py.MPI -O -i ../mmpbsa.in -rp ../protein.prmtop -lp ../ligand.prmtop -cp ../complex.prmtop -eo Interaction_energy.csv -y com*.dcd >> mmpbsa_output
  if [ -e "${input_dir}/index.pdbid/Results/run${NRUN}/FINAL_RESULTS_MMPBSA.dat" ]; then
      if grep -q 'DELTA TOTAL' "${input_dir}/index.pdbid/Results/run${NRUN}/FINAL_RESULTS_MMPBSA.dat"; then
          echo " mmpbsa_status TRUE" >> ${input_dir}/index.pdbid/Result_status
      else
          echo " mmpbsa_status FALSE" >> ${input_dir}/index.pdbid/Result_status
      fi
  else
      echo " FINAL_RESULTS_MMPBSA.dat not found" >> ${input_dir}/index.pdbid/Result_status
  fi
  rsync -azp ${input_dir}/index.pdbid/Result_status ada:${home_dir}/index.pdbid/
  rm ${store_dir}/index.pdbid/Results/run${NRUN}/plum.dcd
  rsync -azp ${input_dir}/index.pdbid/run${NRUN} ada:$store_dir/index.pdbid/
  rm ${store_dir}/index.pdbid/run${NRUN}/complex_solvated.prmtop
  rm ${store_dir}/index.pdbid/run${NRUN}/sim1.dcd
  rm ${store_dir}/index.pdbid/run${NRUN}/cpptraj_starting.in
  rm ${store_dir}/index.pdbid/run${NRUN}/pull.py
  #rm -r ${input_dir}/index.pdbid/run${NRUN}/
done
rsync -azp ${input_dir}/index.pdbid ada:$store_dir/
rsync -azp ${input_dir}/peptide_input.txt ada:$store_dir/
#rm -r ${input_dir}/index.pdbid
