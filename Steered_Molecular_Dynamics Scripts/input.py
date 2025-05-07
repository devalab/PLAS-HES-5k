from prody import *
import parmed as pmd
import numpy as np
import os
import argparse
import subprocess
from os.path import exists
import sys
from colorama import Fore, Style
import csv


ada_input_location='/path/to/plas20k_outputs'
ada_user='username@ada'
password='password'
input_path = '/scratch'
scipt_dir = '/path/to/smd_scripts'

BASH = "/bin/bash" 
    
def parameter():
    comp = parsePDB('sim1_last.pdb')
    p = comp.select('protein')
    a = comp.select('protein within 5 of resname LIG')
    a.getSerials()
    p_serial=list(a.getSerials())
    l = comp.select('resname LIG')
    l_serial=list(l.getSerials())
    rg = calcGyradius(p)
    cent_p = calcCenter(p)
    cent_l = calcCenter(l)
    dc2c = calcDistance(cent_p,cent_l)
    dist = (calcDistance(calcCenter(a),calcCenter(l)).round(3))*0.1
    return p_serial,l_serial,dist,rg,dc2c
   
def string_replace(f,search,replace):
    fin = open(f, "rt")
    data = fin.read()
    data = data.replace(search, replace)
    fin.close()
    fin = open(f, "wt")
    fin.write(data)
    fin.close()
    
def ligand_status(index,pdbid):
    file_path = f'{input_path}/{index}.{pdbid}/ligand.pdb'
    file_exists = os.path.exists(file_path)
    return file_exists

def prepare_input():
    print("Preparing input files...")
    
    command = f"""sshpass -p {password} scp -r {ada_user}:{ada_input_location}/{index}.{pdbid}* {input_path}"""
    subprocess.run(command, shell=True, check=True, executable=BASH)
    print("Input files fetched")
    
    os.chdir(f"{input_path}")
    
    command = f"""
        if [ -f "{index}.{pdbid}.tar" ]
        then
            tar -xf {index}.{pdbid}.tar
            rm -rf {index}.{pdbid}.tar
        fi
    """
    subprocess.run(command, shell=True, check=True, executable=BASH)
    print("Input files untared")

    os.chdir(f"{input_path}/{index}.{pdbid}")
    os.system('rm *.sh cpptraj* heat* mmpbsa* minimize* min_* *.xml *.log')

    complex_pdb = parsePDB('complex_solvated.pdb')
    ligand = complex_pdb.select('resname LIG')
    if ligand:
        ligand_extract = "ligand.pdb"
        prody.writePDB(ligand_extract, ligand)
        print ("ligand extracted")
    else:
        try:
           f = open('{0}/peptide_input.txt'.format(input_path),'a')
           L='{} {}\n'.format(index,pdbid)
           f.write(L)
           f.close()
           os.system('rm -r {}/{}.{}'.format(input_path,index,pdbid))
           os.system('rm {}/{}.{}.txt'.format(input_path,index,pdbid))
        except FileNotFoundError:
           f = open('{0}/peptide_input.txt'.format(input_path),'w')
           L='{} {}\n'.format(index,pdbid)
           f.write(L)
           f.close()
           os.system('rm -r {}/{}.{}'.format(input_path,index,pdbid))
           os.system('rm {}/{}.{}.txt'.format(input_path,index,pdbid))
        
        print("ligand not found, peptide present")
        sys.exit()
    os.system(f'cp {input_path}/{index}.{pdbid}/2-water/2.complex_solvated.prmtop {input_path}/{index}.{pdbid}/')

if __name__ == "__main__":
   parser = argparse.ArgumentParser(description="Run scripts")
   parser.add_argument("index", type=str, nargs="?")
   parser.add_argument("pdbid", type=str, nargs="?")
   parser.add_argument("-r", "--read_input", action="store_true")

   if parser.parse_args().read_input:
       input_list = f"{input_path}/input.txt"
   elif parser.parse_args().index and parser.parse_args().pdbid:
       args = parser.parse_args()
       index = args.index
       pdbid = args.pdbid
       input_list = f"{input_path}/{index}.{pdbid}.txt"
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
           #input_path = f"{input_path}/{index}.{pdbid}"
            
           print('###############################################')
           message = 'Processing index {} and pdbid {}'.format(index,pdbid)
           highlighted_message = "{}{}{}".format(Fore.GREEN, message, Style.RESET_ALL)
           print(highlighted_message)

           prepare_input()
           
           if ligand_status(index,pdbid):
               directory_path = f"{input_path}/{index}.{pdbid}"
               os.system(f'rm {input_path}/{index}.{pdbid}.txt')
               bfe_allrun = []
               for i in range(5):
                   os.chdir(f'{directory_path}/Results/2-water/indi_4/run{i}/')
                   bfe = []
                   inside_row = False
                   with open('Interaction_energy.csv','r') as f:
                       csv_file = csv.reader(f, delimiter=',')
                       for row in csv_file:
                           if row and row[0] == "DELTA Energy Terms":
                              inside_row = True
                           elif inside_row and row:   
                                bfe.append(row[-1])
                   bfe = bfe[1:]
                   bfe = [float(x) for x in bfe]  
                   bfe_allrun.extend(bfe)
                   bfe_all = np.array(bfe_allrun)
                   #print(bfe_all)
    
               mean = np.mean(bfe_all)    
               std = np.std(bfe_all)
               fav_bfe = std + mean
               mean_output = f"{directory_path}/mean_output.txt"
               data = {"Mean:": mean, "Standard deviation:": std, "Favourable BFE:": fav_bfe}
               with open(mean_output, "w") as file:
                   for name, value in data.items():
                       file.write(f"{name}: {value}\n")
    
               print(mean,std,fav_bfe)
    
               fav_bfe_final = []
               frames = []
               for i in range(0, len(bfe_all), 40): 
                   runs = bfe_all[i:i+40]
                   index = np.argmin(np.abs(runs - fav_bfe))
                   fav_bfe_fin = runs[index]
                   fav_bfe_final = np.append(fav_bfe_final, fav_bfe_fin)
                   frame = index+1
                   frames = np.append(frames, frame)
               frames = frames.astype(int)
               print(fav_bfe_final)
               print(frames)
               data1 = {"Starting BFE values:": fav_bfe_final, "Frames:": frames}
               with open(mean_output, "a") as file:
                    for name, value in data1.items():
                        file.write(f"{name}: {value}\n")
               os.system(f'mv {directory_path}/Results/2-water/indi_4/ {directory_path}/')
               os.chdir(f'{directory_path}')
               os.system('rm -r 2-water Results')
    
               for run in range(len(frames)):
                   frm = frames[run]
                   os.chdir(f'{directory_path}/run{run}/')
                   os.system(f'rm *.xml *.log run*')
                   cur_dir = os.getcwd()
                   os.system(f'cp {directory_path}/complex_solvated.prmtop {cur_dir}/')
                   cpptraj_starting = f"""
                   parm complex_solvated.prmtop
                   trajin sim1.dcd {frm} {frm}
                   trajout sim1.rst7 rst7
                   trajout sim1_last.pdb
                   trajout sim1_last.dcd 
                   run
                   """  
                   cpptraj_file = f'cpptraj_starting.in'  
                   with open(cpptraj_file, 'w') as f:
                       f.write(cpptraj_starting)
                   command1 = f"""
                   cpptraj  -i cpptraj_starting.in
                   """
                   subprocess.run(command1, shell=True, check=True, executable=BASH)
        
                   os.system('cp {}/plumed_template.dat {}/plumed.dat'.format(scipt_dir,cur_dir))
                   os.system('cp {}/pull_template.py {}/pull.py'.format(scipt_dir,cur_dir))
            
                   protein_replace,ligand_replace,distance,rog,c2c=parameter()
                   print('RoG = ',rog)
                   print('P_binding_pocket-L distance = ',distance*10)
                   print('C2C P-L distance = ',c2c)
                   run_output = f"{directory_path}/run{run}/run_output.txt"
                   data2 = {"RoG:": rog, "P_binding_pocket-L distance:": distance*10, "C2C P-L distance:": c2c, "Frame:": frm}
                   with open(run_output, "w") as file:
                       for name, value in data2.items():
                           file.write(f"{name}: {value}\n")
                   
                   p_r = ""
                   for elem in protein_replace:
                       p_r = p_r + str(elem) + ','
                   p_r = p_r[0:len(p_r)-1]
        
                   l_r = ""
                   for elem in ligand_replace:
                       l_r = l_r + str(elem) + ','
                   l_r = l_r[0:len(l_r)-1]
 
                   string_replace('plumed.dat','protein_serial',p_r)
                   string_replace('plumed.dat','ligand_serial',l_r)
                   os.system('sed -i "s/\[//g" plumed.dat')
                   os.system('sed -i "s/\]//g" plumed.dat')
           
                   dist = [distance]
                   for i in range(0,2):
                       string_replace('plumed.dat','cal_d{}'.format(i),str(round(distance,4)))
                       distance=(abs((rog-c2c))*0.1+1.0)+distance   #10angstrom, 1 nm
                       dist.append(distance)    
                       #print(dist)
                   time = round(((dist[1] - dist[0])*10),1)
                   step = (time/2)*1000000
                   stride = step/200 
                   step = step.astype(int)
                   stride = stride.astype(int)
                   step = step.astype(str)
                   stride = stride.astype(str)
                   string_replace('plumed.dat','step_size',step)
                   string_replace('plumed.dat','stride_size',stride)
                   string_replace('pull.py','step_size',step)
                   string_replace('pull.py','stride_size',stride)
                   
           #os.system(f'cp {input_path}/{index}.{pdbid}.txt')
print("process complete")
