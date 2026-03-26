import os
import csv
import pickle
import itertools
import re
import cupy as cp
import numpy as np
import pandas as pd
import dask.dataframe as dd
from collections import OrderedDict
from multiprocessing import Pool
from tqdm import tqdm
from dask.diagnostics import ProgressBar

# ------------------------------
# Protein Parsing Class
# ------------------------------
class ParseProtein():
    def __init__(self, rec_fpath):
        # keep original line spacing (do NOT strip leading spaces)
        with open(rec_fpath) as f:
            # choose startswith to match fixed-column PDB records
            self.lines = [x.rstrip("\n") for x in f if x.startswith("ATOM") or x.startswith("HETATM")]
        self.defined_res = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'TYR', 'TRP', 'SER',
                            'THR', 'CYS', 'MET', 'ASN', 'GLN', 'ASP', 'GLU', 'LYS', 'ARG', 'HIS', 'OTH']
        self.all_res_xyz_list = []
        self.res_list = []
          
    def get_res(self, res):
        return res if res in self.defined_res else "OTH"
        
    def extract_letter(self, ele):
        pattern = re.compile(r'([A-Za-z]+)\d+[+-]')
        match = pattern.match(ele)
        return match.group(1) if match else ele
            
    def parse_receptor(self): 
        sym_pool = []
        num = -1      
        _temp_res_xyz = []
        for line in self.lines:
            # use fixed-column parsing (line already preserves spaces)
            ele = line.split()[-1]
            ele = self.extract_letter(ele)
            if ele == "H":
                continue
            num += 1
            res = line[17:20].strip()
            res = self.get_res(res)
            sym = line[17:27].strip()
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            if num == 0:
                self.res_list.append(res)
                sym_pool.append(sym)
                _temp_res_xyz.append([x, y, z])
            else:
                if sym == sym_pool[-1]:
                    _temp_res_xyz.append([x, y, z])
                else:
                    arr = np.array(_temp_res_xyz, dtype=float) * 0.1
                    # ensure (N,3) even if N==1
                    if arr.ndim == 1:
                        if arr.size == 3:
                            arr = arr.reshape(1, 3)
                        else:
                            raise ValueError(f"Unexpected coordinate shape from receptor: {arr.shape}")
                    self.all_res_xyz_list.append(arr)
                    _temp_res_xyz = [[x, y, z]]
                    sym_pool.append(sym)
                    self.res_list.append(res)
        # append last residue
        arr = np.array(_temp_res_xyz, dtype=float) * 0.1
        if arr.ndim == 1:
            if arr.size == 3:
                arr = arr.reshape(1, 3)
            else:
                raise ValueError(f"Unexpected coordinate shape from receptor: {arr.shape}")
        self.all_res_xyz_list.append(arr)
        return self

# ------------------------------
# Ligand Parsing Class
# ------------------------------
class ParseLigand():
    def __init__(self, lig_fpath):
        # preserve spacing for fixed-column parsing
        with open(lig_fpath) as f:
            self.lines = [x.rstrip("\n") for x in f if x.startswith("ATOM") or x.startswith("HETATM")]
        self.defined_eles = ['H', 'C', 'O', 'N', 'P', 'S', 'Hal', 'DU']
        self.hal_ele = ["F", "Cl", "Br", "I"]
        self.lig_ele_list = []
        self.lig_xyz_array = np.array([])
    
    def get_ele(self, ele):
        if ele in self.defined_eles:
            return ele
        elif ele in self.hal_ele:
            return "Hal"
        else:
            return "DU"
        
    def extract_letter(self, ele):
        pattern = re.compile(r'([A-Za-z]+)\d+[+-]')
        match = pattern.match(ele)
        return match.group(1) if match else ele
            
    def parse_ligand(self):
        lig_xyz = []
        for line in self.lines:
            x = float(line[30:38].strip())
            y = float(line[38:46].strip())
            z = float(line[46:54].strip())
            ele = line.split()[-1]
            ele = self.extract_letter(ele)
            ele = self.get_ele(ele)
            self.lig_ele_list.append(ele)
            lig_xyz.append([x, y, z])
        arr = np.array(lig_xyz, dtype=float) * 0.1
        # make sure ligand array is (N,3) even for single-atom ligands
        if arr.ndim == 1:
            if arr.size == 3:
                arr = arr.reshape(1, 3)
            else:
                raise ValueError(f"Unexpected ligand coords shape: {arr.shape}")
        self.lig_xyz_array = arr
        return self

# ------------------------------
# GPU-based Distance Computation
# ------------------------------
def gpu_cdist(x, y):
    x = cp.asarray(x)
    y = cp.asarray(y)
    diff = x[:, None, :] - y[None, :, :]
    return cp.linalg.norm(diff, axis=-1)

# ------------------------------
# Feature Extraction Class (GPU-Accelerated)
# ------------------------------
class GetFeatures():
    def __init__(self, rec, lig, shell, keys):
        self.rec = rec
        self.lig = lig
        self.N_shell = shell
        self.keys = keys

    def count_contacts(self):
        lig_xyz = cp.asarray(self.lig.lig_xyz_array)
        res_xyz_list = [cp.asarray(res) for res in self.rec.all_res_xyz_list]
        res_atom_pairs = []
        res_atom_dist = []
        
        for res, res_xyz in zip(self.rec.res_list, res_xyz_list):
            pairs = [f"{res}_{ele}" for ele in self.lig.lig_ele_list]
            dists = cp.min(gpu_cdist(lig_xyz, res_xyz), axis=1)
            res_atom_pairs.extend(pairs)
            res_atom_dist.extend(dists)
        
        res_atom_dist = cp.array(res_atom_dist)
        outermost = 0.05 * (self.N_shell + 1)
        ncutoffs = cp.linspace(0.1, outermost, self.N_shell)
        results = []
        prev_counts = cp.zeros_like(res_atom_dist, dtype=cp.int32)

        for cutoff in ncutoffs:
            contact_bool = (res_atom_dist <= cutoff).astype(cp.int32)
            onion_counts = contact_bool - prev_counts
            prev_counts = contact_bool
            d = OrderedDict.fromkeys(self.keys, 0)
            for pair, count in zip(res_atom_pairs, onion_counts.tolist()):
                if pair in d:
                    d[pair] += count
            results.append(cp.array(list(d.values())))
        
        return cp.concatenate(results).get()
    
# ------------------------------
# Optimized Standardization
# ------------------------------
def standardize_features_sequential(input_file, output_file):

    # Read CSV with a smaller blocksize to force multiple partitions
    df = dd.read_csv(input_file, blocksize="32MB")

    # Compute global mean & std with a progress bar
    with ProgressBar():
        mean = df.drop("pdb", axis=1).mean().compute()
        std = df.drop("pdb", axis=1).std().compute()

    # Function to standardize a partition
    def standardize_partition(partition):
        labels = partition["pdb"]
        features = (partition.drop("pdb", axis=1) - mean) / std
        features.insert(0, "pdb", labels)
        return features

    # Convert partitions to delayed objects for sequential processing
    partitions = df.to_delayed()

    # Process each partition one by one and write to CSV sequentially
    header_written = False
    with tqdm(total=len(partitions), desc="Standardizing Features", unit="partition") as pbar:
        for part in partitions:
            partition_df = part.compute()
            std_partition = standardize_partition(partition_df)
            
            # Determine file open mode: write (w) for first partition, append (a) for subsequent ones
            mode = "w" if not header_written else "a"
            # Write header only once
            std_partition.to_csv(output_file, mode=mode, index=False, header=not header_written, line_terminator="\n")
            header_written = True
            pbar.update(1)
    
    print("Standardization complete. Saved to:", output_file)

# ------------------------------
# Parallel Pickle Saving
# ------------------------------
def save_pickle(output_csv, out_pickle_dir, chunk_size=10000):
    os.makedirs(out_pickle_dir, exist_ok=True)
    df = dd.read_csv(output_csv)
    
    num_partitions = df.npartitions
    with tqdm(total=num_partitions, desc="Saving Pickle Files", unit="partition") as pbar:
        for i, partition in enumerate(df.to_delayed()):
            chunk = partition.compute()  # Compute one partition at a time
            with open(os.path.join(out_pickle_dir, f"part_{i}.pkl"), "wb") as f:
                pickle.dump(chunk, f)
            pbar.update(1)
    
    print("Pickle files saved to:", out_pickle_dir)

# ------------------------------
# Parallel Feature Extraction
# ------------------------------

def process_input_line_wrapper(args):
    return process_input_line(*args)

def process_input_line(line, shell, keys):
    pdb_id, rec_fpath, lig_fpath = line.split()
    try:
        rec = ParseProtein(rec_fpath).parse_receptor()
        lig = ParseLigand(lig_fpath).parse_ligand()
        feat = GetFeatures(rec, lig, shell, keys)
        result = feat.count_contacts()
        return [pdb_id] + result.tolist()
    except Exception as e:
        print("Error processing", line, e)
        return None

# ------------------------------
# Main Execution Function
# ------------------------------
def main():
    import argparse
    parser = argparse.ArgumentParser(description="Optimized Protein-Ligand Feature Extraction")
    parser.add_argument("-inp", type=str, required=True, help="Input file")
    parser.add_argument("-out", type=str, required=True, help="Output CSV file")
    parser.add_argument("-out_pickle", type=str, required=True, help="Output Pickle directory")
    parser.add_argument("-shells", type=int, default=62, help="Number of shells")
    args = parser.parse_args()
    
    lig_defined_ele = ['H', 'C', 'O', 'N', 'P', 'S', 'Hal', 'DU']
    rec_defined_res = ['GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'PRO', 'PHE', 'TYR', 'TRP', 'SER', 'THR', 'CYS', 'MET', 'ASN', 'GLN', 'ASP', 'GLU', 'LYS', 'ARG', 'HIS', 'OTH']
    keys = ["_".join(x) for x in itertools.product(rec_defined_res, lig_defined_ele)]
    header = ["pdb"] + [f"{n}_{i}" for i, n in enumerate(keys * args.shells)]
    
    with open(args.inp) as f:
        inputs = [line.strip() for line in f.readlines() if line.strip() and not line.startswith("#")]
    
    with open("/scratch/yashas.b/temp_expenergy_features.csv", "w", newline="") as fout:
        writer = csv.writer(fout)
        writer.writerow(header)
        with Pool(os.cpu_count()) as p:
            with tqdm(total=len(inputs), desc="Processing", unit="task") as pbar:
                for res in p.imap_unordered(process_input_line_wrapper, [(line, args.shells, keys) for line in inputs]):
                    if res:
                        writer.writerow(res)
                    pbar.update(1)
    
    print("Feature extraction complete. Standardizing...")
    standardize_features_sequential("/scratch/yashas.b/temp_expenergy_features.csv", args.out)
    save_pickle(args.out, args.out_pickle)

if __name__ == "__main__":
    main()
