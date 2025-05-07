#!/usr/bin/env python
from prody import *
import numpy as np
import os
import argparse
import sys
import subprocess
import concurrent.futures
from os.path import exists
from colorama import Fore, Style
from concurrent.futures import ThreadPoolExecutor
from google.oauth2.service_account import Credentials
from googleapiclient.discovery import build



absolute_path='Path/to/Trajectory_Validation_folder/'

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

            
            os.chdir(f"{input_path}")
            
            message = 'Processing index {} and pdbid {} initial structure'.format(index,pdbid)
            highlighted_message = "{}{}{}".format(Fore.YELLOW, message, Style.RESET_ALL)
            print(highlighted_message)
#                #VMD Startup Command
#            command = f"""/Applications/vmd.app/Contents/MacOS/startup.command -e ../state_initial.vmd"""
#            subprocess.run(command,
#                            shell=True,
#                            check=True,
#                            executable=BASH,
#                        )              
                
            for i in range(5):
                
                message = 'Processing index {} and pdbid {} ---------> run {}'.format(index,pdbid,i)
                highlighted_message = "{}{}{}".format(Fore.GREEN, message, Style.RESET_ALL)
                print(highlighted_message)
                
                                
                os.chdir(f"{absolute_path}/{index}.{pdbid}/run{i}")
                                                #VMD Startup Command
                command = f"""/usr/local/bin/vmd -e ../../state_1.vmd"""
                subprocess.run(command, shell=True, check=True, executable=BASH)
                
                
                                                #VMD Startup Command
                command = f"""/usr/local/bin/vmd -e ../../state.vmd"""
                subprocess.run(command, shell=True, check=True, executable=BASH)
                command = f"""cp {absolute_path}/Distance_Sigmoid_RMSD.py {absolute_path}/{index}.{pdbid}/run{i}/"""
                subprocess.run(command, shell=True, check=True, executable=BASH)
                command = f"""cp {input_path}/complex_solvated.prmtop {absolute_path}/{index}.{pdbid}/run{i}/"""
                subprocess.run(command, shell=True, check=True, executable=BASH)
                command = f"""python Distance_Sigmoid_RMSD.py"""
                subprocess.run(command, shell=True, check=True, executable=BASH)
                

            os.chdir(f"{absolute_path}")
            os.system('rm -rf {}.{}'.format(index,pdbid))
