#!/bin/bash

## sh Complete_simulation_setup_hes.sh 4018 2v2v 12 d4 ihub
index=$1
pdbid=$2
Ncores=$3
Account=$4
Partition=$5

input_dir="/home2/ananaya.jain/hes_input"
script_dir="/home2/ananaya.jain/smd_scripts_2"

mkdir -p $input_dir/${index}.${pdbid}
cp $script_dir/submit.sh $input_dir/${index}.${pdbid}/
sed -i "s/pdbid/${pdbid}/g" $input_dir/${index}.${pdbid}/submit.sh
sed -i "s/index/${index}/g" $input_dir/${index}.${pdbid}/submit.sh
sed -i "s/job_name/${index}_${pdbid}_${Account}/g" $input_dir/${index}.${pdbid}/submit.sh
sed -i "s/account/${Account}/g" $input_dir/${index}.${pdbid}/submit.sh
sed -i "s/partition/${Partition}/g" $input_dir/${index}.${pdbid}/submit.sh
sed -i "s/ncores/${Ncores}/g" $input_dir/${index}.${pdbid}/submit.sh

cd $input_dir/${index}.${pdbid}

sbatch submit.sh 
