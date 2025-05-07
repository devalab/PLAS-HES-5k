#!/usr/bin/env python
import numpy as np
import os
import argparse
import sys
import subprocess
import concurrent.futures
from os.path import exists
from colorama import Fore, Style


ada_input_location='/Path/to/hes_output/'
password='enter password' 
ada_user='Ada_username'
absolute_path='/path/to/Trajectory_Validation_folder'


BASH = "/bin/bash"
ZSH = "/bin/zsh"




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run scripts")
    parser.add_argument("index", type=str, nargs="?")
    parser.add_argument("pdbid", type=str, nargs="?")
    parser.add_argument("-r", "--read_input", action="store_true")

    if parser.parse_args().read_input:
        input_list = f"{absolute_path}/input.txt"
    elif parser.parse_args().index and parser.parse_args().pdbid:
        args = parser.parse_args()
        index = args.index
        pdbid = args.pdbid
        input_list = f"{absolute_path}/{index}.{pdbid}.txt"
        with open(input_list, "w", encoding="utf-8") as f:
            f.write(f"{index} {pdbid}")
    else:
        print("Please enter the index and pdbid, or provide an input file")
        sys.exit()


    with open(input_list, "r", encoding="utf-8") as id_list:
        for line in id_list:
            if not line.strip():
                continue
            index, pdbid = line.strip().split()[0],line.strip().split()[1]

            if not index or not pdbid:
                continue
            input_path = f"{absolute_path}/{index}.{pdbid}"

            print('###############################################')
            message = 'Processing index {} and pdbid {}'.format(index,pdbid)
            highlighted_message = "{}{}{}".format(Fore.GREEN, message, Style.RESET_ALL)
            print(highlighted_message)

            os.system('mkdir {}'.format(input_path))
            os.chdir(f"{input_path}")
            command = f"""sshpass -p {password} scp -r {ada_user}:{ada_input_location}/{index}.{pdbid}/Results/2.complex_solvated.prmtop {input_path}"""
            subprocess.run(command, shell=True, check=True, executable=BASH)

            try:
                command = f"""sshpass -p {password} scp -r {ada_user}:{ada_input_location}/{index}.{pdbid}/complex_solvated.pdb {input_path}"""
                subprocess.run(command, shell=True, check=True, executable=BASH)
            except subprocess.CalledProcessError:
                try:
                    command = f"""sshpass -p {password} scp -r {ada_user}:{ada_input_location}/{index}.{pdbid}/complex_solvated.prmtop {input_path}"""
                    subprocess.run(command, shell=True, check=True, executable=BASH) 
                except subprocess.CalledProcessError:
                    print(Fore.RED + "Initial Structure not found in the database" + Style.RESET_ALL)
                    continue  

            try:
                command = f"""sshpass -p {password} scp -r {ada_user}:{ada_input_location}/{index}.{pdbid}/complex_solvated.prmtop {input_path}"""
                subprocess.run(command, shell=True, check=True, executable=BASH)
            except subprocess.CalledProcessError:
                try:
                    command = f"""sshpass -p {password} scp -r {ada_user}:{ada_input_location}/{index}.{pdbid}/complex_solvated.prmtop {input_path}"""
                    subprocess.run(command, shell=True, check=True, executable=BASH) 
                except subprocess.CalledProcessError:
                    print(Fore.RED + "initial prmtop is not found in the database" + Style.RESET_ALL)
                    continue
                    
                    
            def download_files(i):
                run_folder = os.path.join(input_path, f"run{i}")
                os.makedirs(run_folder, exist_ok=True)
                message = f'Processing index {index} and pdbid {pdbid} ---------> run {i}'
                highlighted_message = "{}{}{}".format(Fore.YELLOW, message, Style.RESET_ALL)
                print(highlighted_message)
                scp_commands = [
                    f"sshpass -p {password} scp -r {ada_user}:{ada_input_location}/{index}.{pdbid}/run{i}/sim1_last.pdb {run_folder}",
                    f"sshpass -p {password} scp -r {ada_user}:{ada_input_location}/{index}.{pdbid}/run{i}/plum.dcd {run_folder}",
                    f"sshpass -p {password} scp -r {ada_user}:{ada_input_location}/{index}.{pdbid}/run{i}/plumed.dat {run_folder}",
                    f"sshpass -p {password} scp -r {ada_user}:{ada_input_location}/{index}.{pdbid}/Results/run{i}/Interaction_energy.csv {run_folder}"
                ]
                for command in scp_commands:
                    subprocess.run(command, shell=True, check=True, executable="/bin/bash")

            with concurrent.futures.ThreadPoolExecutor() as executor:
                runs = range(5)
                executor.map(download_files, runs)
                
            os.chdir(f"{absolute_path}")
