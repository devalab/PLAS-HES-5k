import os
import subprocess
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit
import shutil
import concurrent.futures
import glob
import csv
from collections import defaultdict
import time
import logging
import datetime
import traceback
import fnmatch

BASH = "/bin/bash"

# ========== Configuration Section ==========
REMOTE_DIRS = [
    "path/to/files"
]

LOCAL_DIR = "/path/to/copyingHES"
TRANSFER_DIR = "/path/to/transferHES"
INPUT_TXT = os.path.join(LOCAL_DIR, "input.txt")
ADA_USER = "dAda_username"
PASSWORD = "password"
HES5K_DIR = "/path/to/store_final.tar.gz"        # For final .tar.gz files
HES_CSV_PATH = os.path.join(TRANSFER_DIR, "HES.csv")
LOG_FILE = os.path.join(LOCAL_DIR, "log.txt")
# ===========================================

# Set up logging
logging.basicConfig(
    filename=LOG_FILE,
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

# Add console output for logs
console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
console.setFormatter(formatter)
logging.getLogger('').addHandler(console)

def log_start():
    """Log script start with timestamp"""
    logging.info("=" * 80)
    logging.info(f"STARTING AUTO_TRANSFER SCRIPT - {datetime.datetime.now()}")
    logging.info(f"Input file: {INPUT_TXT}")
    logging.info(f"Local directory: {LOCAL_DIR}")
    logging.info(f"Transfer directory: {TRANSFER_DIR}")
    logging.info("=" * 80)

def sigmoid(x, a, b, c, d):
    return a / (1 + np.exp(-c * (x - d))) + b

def fit_sigmoid(x, y):
    initial_guess = [np.max(y)-np.min(y), np.min(y), 0.01, (np.max(x)+np.min(x))/2]
    bounds = ([0, -np.inf, 0, -np.inf], [np.inf, np.inf, np.inf, np.inf])
    popt, _ = curve_fit(sigmoid, x, y, p0=initial_guess, bounds=bounds, maxfev=3000)
    return sigmoid(x, *popt)

def process_csv(input_csv, final_frames, output_csv):
    """Process all four energy terms in CSV and save filtered version"""
    logging.info(f"Processing CSV: {input_csv} -> {output_csv}")
    
    try:
        with open(input_csv, 'r') as f:
            all_lines = f.readlines()

        # Parse file structure
        global_headers = []
        sections = []
        current_section = None
        current_headers = []
        current_data = []

        for line in all_lines:
            stripped = line.strip()
            if any(stripped.startswith(term) for term in [
                "Complex Energy Terms",
                "Receptor Energy Terms", 
                "Ligand Energy Terms",
                "DELTA Energy Terms"
            ]):
                if current_section:
                    sections.append({'headers': current_headers, 'data': current_data})
                current_section = stripped.split(',')[0]
                current_headers = [line]
                current_data = []
            else:
                if current_section:
                    if len(current_headers) < 2:
                        current_headers.append(line)
                    else:
                        if all(c == ',' for c in stripped.replace(',', '')):
                            sections.append({'headers': current_headers, 'data': current_data})
                            current_section = None
                            current_headers = []
                            current_data = []
                        else:
                            current_data.append(line)
                else:
                    global_headers.append(line)

        if current_section:
            sections.append({'headers': current_headers, 'data': current_data})

        # Process each section
        processed = []
        for section in sections:
            selected = []
            for new_frame, orig_frame in enumerate(final_frames):
                if orig_frame < len(section['data']):
                    parts = section['data'][orig_frame].strip().split(',')
                    parts[0] = str(new_frame)  # Renumber frame
                    selected.append(','.join(parts) + '\n')
            
            processed.append({
                'headers': section['headers'],
                'data': selected
            })

        # Write new CSV
        with open(output_csv, 'w') as f:
            f.writelines(global_headers)
            for i, section in enumerate(processed):
                f.writelines(section['headers'])
                f.writelines(section['data'])
                if i < len(processed)-1:
                    cols = len(section['headers'][1].split(',')) if len(section['headers'])>1 else 0
                    f.write(','*cols + '\n')
        
        logging.info(f"Successfully processed CSV with {len(final_frames)} frames")
        return True
    except Exception as e:
        logging.error(f"Error in process_csv: {str(e)}")
        logging.error(traceback.format_exc())
        return False

def ensure_remote_dir_exists():
    """Create the remote tmp_HES5K directory if it doesn't exist"""
    logging.info(f"Ensuring remote directory exists: {TMP_HES_DIR}")
    mkdir_cmd = f"sshpass -p {PASSWORD} ssh {ADA_USER} 'mkdir -p {TMP_HES_DIR}'"
    try:
        subprocess.run(mkdir_cmd, shell=True, check=True)
        logging.info(f"Successfully created or confirmed remote directory: {TMP_HES_DIR}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to create remote directory: {str(e)}")
        return False

def process_delta_energy(pdbid, index_dir):
    """Process DELTA energy terms from all runs and create CSVs"""
    logging.info(f"Processing DELTA energy for {pdbid}")
    
    # Ensure remote directory exists
    if not ensure_remote_dir_exists():
        logging.error("Failed to ensure remote directory exists")
        return False
    
    delta_columns = ['DELTA_TOTAL', 'VDWAALS', 'EEL', 'EPB', 'ENPOLAR']
    run_means = defaultdict(list)  # Stores mean of each run
    intermediate_data = []  # For intermediate CSV

    # Process all 5 runs
    for run_num in range(5):
        csv_path = os.path.join(index_dir, "Results", f"run{run_num}",
                              f"Interaction_energy_40_{pdbid}_run{run_num}.csv")
        if not os.path.exists(csv_path):
            logging.warning(f"CSV not found for run{run_num}: {csv_path}")
            continue

        try:
            # Find DELTA section dynamically
            with open(csv_path, 'r') as f:
                reader = csv.reader(f)
                delta_found = False
                header_found = False
                delta_data = []
                
                for row in reader:
                    # Find DELTA Energy Terms section
                    if not delta_found and any("DELTA Energy Terms" in cell for cell in row):
                        delta_found = True
                        continue
                    
                    if delta_found and not header_found:
                        # Look for column headers
                        if row and row[0] == "Frame #":
                            header_found = True
                            continue
                    
                    # Collect all data rows for DELTA section
                    if delta_found and header_found:
                        # Stop if we reach a blank row or a new section
                        if not row or not row[0] or row[0].strip() == '':
                            break
                        
                        try:
                            # Check if this is a data row (should start with a frame number)
                            frame_num = int(row[0])
                            delta_data.append(row)
                        except (ValueError, IndexError):
                            # Not a frame row, might be end of section
                            break
                
                if not delta_data:
                    raise ValueError("Could not find DELTA Energy Terms data")
                
                # Check if exactly 40 frames are present
                if len(delta_data) != 40:
                    logging.warning(f"Run {run_num} has {len(delta_data)} frames, expected 40. Skipping.")
                    continue
                
                logging.info(f"Found 40 frames in DELTA section for run{run_num}")

            # Define correct column indices
            column_indices = {
                'VDWAALS': 7,  # Column H
                'EEL': 8,      # Column I
                'EPB': 11,     # Column L
                'ENPOLAR': 12, # Column M
                'DELTA_TOTAL': 16  # Column Q
            }

            # Calculate mean values for each column (40 frames)
            run_avgs = {}
            for col_name, col_idx in column_indices.items():
                col_values = []
                for row in delta_data:
                    if col_idx < len(row) and row[col_idx]:
                        try:
                            col_values.append(float(row[col_idx]))
                        except ValueError:
                            logging.warning(f"Invalid value in {csv_path} for {col_name}: {row[col_idx]}")
                            col_values.append(np.nan)
                
                # Validate 40 values per column
                if len(col_values) != 40:
                    logging.warning(f"Column {col_name} in run{run_num} has {len(col_values)} values. Skipping run.")
                    run_avgs = None
                    break
                
                run_avgs[col_name] = np.nanmean(col_values)
            
            if run_avgs is None:
                continue  # Skip invalid run

            # Populate data structures
            run_entry = {'Run': run_num}
            for col in delta_columns:
                run_entry[col] = run_avgs[col]
                run_means[col].append(run_avgs[col])
            
            intermediate_data.append(run_entry)
            logging.info(f"Calculated mean values for run{run_num}: {run_avgs}")

        except Exception as e:
            logging.error(f"Error processing {csv_path}: {str(e)}")
            logging.error(traceback.format_exc())
            continue

    # Proceed only if there are valid run means
    if not all(len(run_means[col]) > 0 for col in delta_columns):
        logging.error(f"No valid run data found for {pdbid}")
        return False

    # Create local temp directory
    local_temp_dir = os.path.join(LOCAL_DIR, "temp_delta_files")
    try:
        os.makedirs(local_temp_dir, exist_ok=True)
        logging.info(f"Created local temp directory: {local_temp_dir}")
    except Exception as e:
        logging.error(f"Failed to create local temp directory: {str(e)}")
        return False
    
    local_intermediate_path = os.path.join(local_temp_dir, f"DELTA_Energy_terms_{pdbid}.csv")
    
    try:
        # Create intermediate CSV with run data
        with open(local_intermediate_path, 'w') as f:
            writer = csv.DictWriter(f, fieldnames=['Run'] + delta_columns)
            writer.writeheader()
            writer.writerows(intermediate_data)
            
            # Add MEAN row using run_means
            mean_row = {'Run': 'MEAN'}
            for col in delta_columns:
                if run_means[col]:  # Check if data exists
                    mean_row[col] = np.nanmean(run_means[col])
                else:
                    mean_row[col] = np.nan
            writer.writerow(mean_row)
        
        logging.info(f"Created intermediate CSV: {local_intermediate_path}")

        # Copy to remote directory
        remote_intermediate_path = f"{TMP_HES_DIR}/DELTA_Energy_terms_{pdbid}.csv"
        copy_cmd = (
            f"sshpass -p {PASSWORD} scp "
            f"{local_intermediate_path} "
            f"{ADA_USER}:{remote_intermediate_path}"
        )
        subprocess.run(copy_cmd, shell=True, check=True)
        logging.info(f"Copied intermediate CSV to {remote_intermediate_path}")

    except Exception as e:
        logging.error(f"Error creating/transferring intermediate file: {str(e)}")
        return False

    try:
        # Create/append to HES.csv
        os.makedirs(TRANSFER_DIR, exist_ok=True)
        
        experimental = get_experimental_value(pdbid)
        hes_data = {
            'PDB_ID': pdbid,
            'Experimental': experimental,
            **{col: np.nanmean(run_means[col]) for col in delta_columns}
        }
        
        header = not os.path.exists(HES_CSV_PATH)
        with open(HES_CSV_PATH, 'a') as f:
            writer = csv.DictWriter(f, fieldnames=['PDB_ID', 'Experimental'] + delta_columns)
            if header:
                writer.writeheader()
            writer.writerow(hes_data)
        
        logging.info(f"Updated HES.csv at {HES_CSV_PATH}")
        return True
    except Exception as e:
        logging.error(f"Error creating HES.csv: {str(e)}")
        logging.error(traceback.format_exc())
        return False
        
def get_experimental_value(pdbid):
    """Get experimental value from extended_PLAS20K.csv"""
    plas_csv = os.path.join(LOCAL_DIR, "extended_PLAS20K.csv")
    try:
        df = pd.read_csv(plas_csv)
        row = df[df['PDB_ID'] == pdbid.upper()]
        if not row.empty:
            val = row['Experimental'].values[0]
            logging.info(f"Found experimental value for {pdbid}: {val}")
            return val
        else:
            logging.warning(f"No experimental value found for {pdbid}")
            return "N/A"
    except Exception as e:
        logging.error(f"Error reading experimental value: {str(e)}")
        logging.error(traceback.format_exc())
    return "N/A"

def setup_transfer_structure(index, pdbid):
    """Enhanced file copying with path validation for HES5K structure"""
    logging.info(f"Setting up transfer structure for {index}.{pdbid}")
    
    source_dir = os.path.join(LOCAL_DIR, f"{index}.{pdbid}")
    target_dir = os.path.join(TRANSFER_DIR, f"{index}.{pdbid}")
    
    try:
        # Create HES directory structure
        hes_dir = os.path.join(target_dir, 'HES')
        os.makedirs(hes_dir, exist_ok=True)
        
        # Create HES 2.complex_solvated.pdb
        create_complex_solvated_pdb(source_dir, pdbid)
        
        # Copy required files to HES - top level
        for fname in ['2.complex_solvated.pdb', '2.complex_solvated.prmtop', 'complex_solvated.prmtop', 'complex_solvated.pdb']:
            src = os.path.join(source_dir, fname)
            dst = os.path.join(hes_dir, fname)
            if os.path.exists(src):
                shutil.copy2(src, dst)
                logging.info(f"Copied {fname} to HES directory")
            else:
                logging.warning(f"File not found: {src}")
        
        # Setup run directories
        for run_num in range(5):
            src_run_dir = os.path.join(source_dir, "Results", f"run{run_num}")
            dest_run_dir = os.path.join(hes_dir, f"run{run_num}")
            os.makedirs(dest_run_dir, exist_ok=True)
            
            # Copy CSV and rename to Interaction_energy.csv
            src_csv = os.path.join(src_run_dir, f"Interaction_energy_40_{pdbid}_run{run_num}.csv")
            dst_csv = os.path.join(dest_run_dir, "Interaction_energy.csv")
            
            if os.path.exists(src_csv):
                shutil.copy2(src_csv, dst_csv)
                logging.info(f"Copied and renamed CSV for run{run_num}")
            else:
                logging.warning(f"CSV file not found: {src_csv}")
            
            # Generate PDB instead of DCD
            src_dcd = os.path.join(src_run_dir, f"40_com_wat2_{pdbid}_run{run_num}.dcd")
            dst_pdb = os.path.join(dest_run_dir, f"40_com_wat2_{pdbid}_run{run_num}.pdb")
            
            if os.path.exists(src_dcd):
                # Convert DCD to PDB
                prmtop_path = os.path.join(source_dir, "2.complex_solvated.prmtop")
                if os.path.exists(prmtop_path):
                    cpptraj_input = os.path.join(src_run_dir, f"dcd2pdb_{run_num}.in")
                    with open(cpptraj_input, 'w') as f:
                        f.write(f"""parm {prmtop_path}
trajin {src_dcd}
trajout {dst_pdb} pdb
run
quit
""")
                    subprocess.run(f"cpptraj -i {cpptraj_input}", shell=True, executable=BASH)
                    logging.info(f"Converted DCD to PDB for run{run_num}")
                    
                    # Clean up
                    if os.path.exists(cpptraj_input):
                        os.remove(cpptraj_input)
                else:
                    logging.warning(f"Prmtop file not found: {prmtop_path}")
            else:
                logging.warning(f"DCD file not found: {src_dcd}")
        
        # Create PLAS directory structure
        plas_dir = os.path.join(target_dir, 'PLAS')
        os.makedirs(plas_dir, exist_ok=True)
        
        # Copy variation2's runs to PLAS root
        for run_num in range(5):
            run_dir = os.path.join(plas_dir, f"run{run_num}")
            os.makedirs(run_dir, exist_ok=True)
            
            # Copy Interaction_energy.csv from remote
            remote_ie_source = f"/zfs/d4/plas20k_interaction_outputs/{index}.{pdbid}/run{run_num}/Interaction_energy.csv"
            local_ie_dest = os.path.join(run_dir, "Interaction_energy.csv")
            copy_ie_cmd = f"sshpass -p {PASSWORD} scp {ADA_USER}:{remote_ie_source} {local_ie_dest}"
            try:
                subprocess.run(copy_ie_cmd, shell=True, check=True)
                logging.info(f"Copied Interaction_energy.csv for run{run_num} to {local_ie_dest}")
            except subprocess.CalledProcessError as e:
                logging.error(f"Error copying Interaction_energy.csv for run{run_num}: {e}")
            
            # Copy expli_2wat_{run_num}.pdb from remote
            remote_pdb_source = f"/zfs/d4/PLAS20k-Cleanedup/chunks/{pdbid}/expli_2wat_{run_num}.pdb"
            local_pdb_dest = os.path.join(run_dir, f"expli_2wat_{run_num}.pdb")
            copy_pdb_cmd = f"sshpass -p {PASSWORD} scp {ADA_USER}:{remote_pdb_source} {local_pdb_dest}"
            try:
                subprocess.run(copy_pdb_cmd, shell=True, check=True)
                logging.info(f"Copied expli_2wat_{run_num}.pdb to {local_pdb_dest}")
            except subprocess.CalledProcessError as e:
                logging.error(f"Error copying expli_2wat_{run_num}.pdb: {e}")
        
        logging.info(f"Transfer structure setup complete for {index}.{pdbid}")
        return True
    except Exception as e:
        logging.error(f"Error setting up transfer structure: {str(e)}")
        logging.error(traceback.format_exc())
        return False
        
def process_index(index_pdbid, base_dir):
    """Process all 5 runs for a single index.pdbid"""
    try:
        index, pdbid = index_pdbid.split()
        logging.info(f"Starting to process {index}.{pdbid}")
        
        index_dir = os.path.join(base_dir, f"{index}.{pdbid}")
        if not os.path.exists(index_dir):
            logging.error(f"Index directory not found: {index_dir}")
            return False
        
        success_count = 0
        for run_num in range(5):  # Process run0-run4
            run_dir = os.path.join(index_dir, "Results", f"run{run_num}")
            cpptraj_file = None
            temp_csv = None
            
            if not os.path.exists(run_dir):
                logging.warning(f"Run directory not found: {run_dir}")
                continue

            logging.info(f"\n🔨 Processing {index}.{pdbid} - run{run_num}...")
            
            try:
                # File paths
                input_csv = os.path.join(run_dir, "Interaction_energy.csv")
                if not os.path.exists(input_csv):
                    logging.error(f"Interaction_energy.csv not found in {run_dir}")
                    continue
                    
                output_dcd = os.path.join(run_dir, f"40_com_wat2_{pdbid}_run{run_num}.dcd")
                output_csv = os.path.join(run_dir, f"Interaction_energy_40_{pdbid}_run{run_num}.csv")
                temp_csv = os.path.join(run_dir, "last_column.csv")

                # Frame selection logic
                logging.info(f"Reading CSV: {input_csv}")
                csv1 = pd.read_csv(input_csv, skiprows=612, header=None)
                logging.info(f"CSV shape: {csv1.shape}")
                
                if csv1.shape[1] == 0:
                    logging.error(f"CSV file {input_csv} is empty or improperly formatted")
                    continue
                
                last_column = csv1.iloc[:, -1]
                last_column.to_csv(temp_csv, index=False, header=False)
                
                data = pd.read_csv(temp_csv, header=None)
                data = pd.to_numeric(data[0], errors='coerce').dropna()
                filtered_data = data[data < 0].reset_index(drop=True)
                filtered_indices = data.index[data < 0].values.astype(int)
                
                logging.info(f"Found {len(filtered_data)} negative values for curve fitting")
                if len(filtered_data) < 40:
                    logging.error(f"Not enough negative values ({len(filtered_data)}) for curve fitting")
                    continue

                x = np.arange(1, len(filtered_data)+1)
                y = filtered_data.values

                logging.info("Fitting sigmoid curve")
                fitted_curve = fit_sigmoid(x, y)
                scaled_curve = (fitted_curve - np.min(fitted_curve))/(np.max(fitted_curve)-np.min(fitted_curve))
                threshold = 0.75

                while True:
                    frame_at_threshold = np.where(scaled_curve >= threshold)[0][0]
                    upper_half_mask = x >= frame_at_threshold
                    if len(x[upper_half_mask]) >= 30: break
                    threshold -= 0.05
                    if threshold < 0.05:
                        logging.warning(f"Threshold dropped below 0.05, using 0.05")
                        frame_at_threshold = np.where(scaled_curve >= 0.05)[0][0]
                        break

                # Frame selection calculations
                lower_half_mask = x < frame_at_threshold
                upper_half_mask = x >= frame_at_threshold
                
                lower_y = y[lower_half_mask]
                upper_y = y[upper_half_mask]
                
                logging.info(f"Frame separation at threshold {threshold}: Lower={len(lower_y)}, Upper={len(upper_y)}")
                
                # Lower binning
                lower_diff = np.abs(lower_y - fitted_curve[lower_half_mask])
                n = len(lower_diff)
                if n < 10:
                    logging.warning(f"Not enough frames in lower half ({n}), adjusting bin count")
                    bin_count = max(1, n)
                else:
                    bin_count = 10
                    
                bin_size = n // bin_count
                remainder = n % bin_count
                bin_sizes = [bin_size + 1 if i < remainder else bin_size for i in range(bin_count)]
                lower_closest_indices = [
                    np.argmin(lower_diff[start:start+size]) + start 
                    for start, size in zip(np.cumsum([0]+bin_sizes[:-1]), bin_sizes)
                ]
                
                # Upper binning
                upper_diff = np.abs(upper_y - fitted_curve[upper_half_mask])
                n = len(upper_diff)
                if n < 30:
                    logging.warning(f"Not enough frames in upper half ({n}), adjusting bin count")
                    bin_count = max(1, n)
                else:
                    bin_count = 30
                    
                bin_size = n // bin_count
                remainder = n % bin_count
                bin_sizes = [bin_size + 1 if i < remainder else bin_size for i in range(bin_count)]
                upper_closest_indices = np.array([
                    np.argmin(upper_diff[start:start+size]) + start 
                    for start, size in zip(np.cumsum([0]+bin_sizes[:-1]), bin_sizes)
                ])
                
                # Get original indices
                marked_lower = filtered_indices[lower_closest_indices]
                marked_upper = filtered_indices[upper_closest_indices + len(lower_y)]
                final_frames = np.concatenate([marked_lower, marked_upper]).astype(int)
                
                logging.info(f"Selected {len(final_frames)} frames: {len(marked_lower)} lower + {len(marked_upper)} upper")

                # Process CSV with all energy terms
                if not process_csv(input_csv, final_frames, output_csv):
                    logging.error(f"Failed to process CSV for {pdbid} run{run_num}")
                    continue

                # Generate DCD
                cpptraj_content = f"""parm {os.path.join(index_dir, "Results", "2.complex_solvated.prmtop")}
trajin {os.path.join(run_dir, "com_wat2.dcd")}
trajout {output_dcd} onlyframes {','.join(map(str, final_frames))}
run
"""
                cpptraj_file = os.path.join(run_dir, f"cpptraj_40_{pdbid}_run{run_num}.in")
                with open(cpptraj_file, 'w') as f:
                    f.write(cpptraj_content)
                
                logging.info(f"Running cpptraj for {pdbid} run{run_num}")
                try:
                    result = subprocess.run(f"cpptraj -i {cpptraj_file}", 
                                          shell=True, 
                                          executable=BASH, 
                                          capture_output=True, 
                                          text=True)
                    if result.returncode != 0:
                        logging.error(f"cpptraj error: {result.stderr}")
                    else:
                        logging.info(f"cpptraj completed successfully")
                        success_count += 1
                except Exception as e:
                    logging.error(f"cpptraj execution error: {str(e)}")
                
            except Exception as e:
                logging.error(f"Error processing run{run_num}: {str(e)}")
                continue

        # After all runs are processed
        if success_count >= 1:
            # Process DELTA energy terms and update HES.csv
            logging.info(f"Processing DELTA energy terms for {pdbid}")
            delta_success = process_delta_energy(pdbid, index_dir)
            if delta_success:
                logging.info(f"Successfully processed DELTA energy terms for {pdbid}")
            else:
                logging.warning(f"Failed to process DELTA energy terms for {pdbid}")
            return delta_success
        else:
            logging.error(f"No successful runs for {pdbid}, skipping DELTA energy processing")
            return False
    
    except Exception as e:
        logging.error(f"Fatal error in process_index: {str(e)}")
        logging.error(traceback.format_exc())
        return False
        
def run_command(command):
    """Run a shell command and return success/failure and output"""
    try:
        result = subprocess.run(command, shell=True, executable=BASH, 
                               capture_output=True, text=True)
        return result.returncode == 0, result.stdout + result.stderr
    except Exception as e:
        return False, str(e)

def create_complex_solvated_pdb(work_dir, pdbid):
    """Create 2.complex_solvated.pdb file from existing pdb or prmtop"""
    try:
        # Check for existing files
        existing_pdb = os.path.join(work_dir, 'complex_solvated.pdb')
        output_pdb = os.path.join(work_dir, '2.complex_solvated.pdb')
        
        # If complex_solvated.pdb exists, use that directly
        if os.path.exists(existing_pdb):
            logging.info(f"Found complex_solvated.pdb, using it to create 2.complex_solvated.pdb")
            
            # Create cpptraj input file to process the PDB
            cpptraj_input = os.path.join(work_dir, "create_pdb.in")
            with open(cpptraj_input, 'w') as f:
                f.write(f"""parm {existing_pdb}
trajin {existing_pdb}
strip :Na+,Cl-
closest 2 :LIG first
trajout {output_pdb} pdb :PROT,:LIG,CLosest
run
quit
""")
            
            # Run cpptraj
            logging.info(f"Running cpptraj to create 2.complex_solvated.pdb")
            cpptraj_cmd = f"cpptraj -i {cpptraj_input}"
            success, output = run_command(cpptraj_cmd)
            
            # If that fails, try alternative approach with just copying the file
            if not success or not os.path.exists(output_pdb) or os.path.getsize(output_pdb) == 0:
                logging.warning(f"cpptraj processing failed: {output}. Trying direct copy.")
                shutil.copy(existing_pdb, output_pdb)
                logging.info(f"Copied {existing_pdb} to {output_pdb}")
                success = True
            
            # Clean up input file
            if os.path.exists(cpptraj_input):
                os.remove(cpptraj_input)
            
            return success
        else:
            logging.warning(f"No complex_solvated.pdb found in {work_dir}")
            return False
    except Exception as e:
        logging.error(f"Error creating 2.complex_solvated.pdb: {str(e)}")
        logging.error(traceback.format_exc())
        return False

def process_final_structure(index, pdbid):
    """Post-processing to ensure the final structure is correct"""
    logging.info(f"Post-processing final structure for {index}.{pdbid}")
    
    target_dir = os.path.join(TRANSFER_DIR, f"{index}.{pdbid}")
    
    try:
        # Process HES directory
        hes_dir = os.path.join(target_dir, 'HES')
        if os.path.exists(hes_dir):
            # Remove any DCD files
            for root, dirs, files in os.walk(hes_dir):
                for file in files:
                    if file.endswith('.dcd'):
                        os.remove(os.path.join(root, file))
                        logging.info(f"Removed DCD file: {os.path.join(root, file)}")
        
        # Process PLAS directory
        plas_dir = os.path.join(target_dir, 'PLAS')
        if os.path.exists(plas_dir):
            # Check for variation2 directory and flatten if needed
            var2_dir = os.path.join(plas_dir, 'variation2')
            if os.path.exists(var2_dir):
                for item in os.listdir(var2_dir):
                    item_path = os.path.join(var2_dir, item)
                    dest_path = os.path.join(plas_dir, item)
                    
                    if os.path.exists(dest_path):
                        if os.path.isdir(dest_path):
                            shutil.rmtree(dest_path)
                        else:
                            os.remove(dest_path)
                    
                    shutil.move(item_path, dest_path)
                    logging.info(f"Moved {item} from variation2 to PLAS root")
                
                # Remove empty variation2 directory
                shutil.rmtree(var2_dir)
                logging.info("Removed variation2 directory")
            
            # Remove variation1 directory
            var1_dir = os.path.join(plas_dir, 'variation1')
            if os.path.exists(var1_dir):
                shutil.rmtree(var1_dir)
                logging.info("Removed variation1 directory")
        
        logging.info(f"Post-processing complete for {index}.{pdbid}")
        return True
    except Exception as e:
        logging.error(f"Error in post-processing: {str(e)}")
        logging.error(traceback.format_exc())
        return False
        
def post_process_structure(target_dir):
    """Clean up and finalize structure to match requirements"""
    logging.info(f"Post-processing directory structure: {target_dir}")
    try:
        # Process HES directory - remove all DCD files
        hes_dir = os.path.join(target_dir, 'HES')
        if os.path.exists(hes_dir):
            for root, dirs, files in os.walk(hes_dir):
                for file in files:
                    if file.endswith('.dcd'):
                        file_path = os.path.join(root, file)
                        os.remove(file_path)
                        logging.info(f"Removed DCD file: {file_path}")
        
        # Ensure all HES run folders have correct CSV naming
        for run_num in range(5):
            run_dir = os.path.join(hes_dir, f"run{run_num}")
            if os.path.exists(run_dir):
                for file in os.listdir(run_dir):
                    if file.startswith("Interaction_energy_40_") and file.endswith(".csv"):
                        src = os.path.join(run_dir, file)
                        dst = os.path.join(run_dir, "Interaction_energy.csv")
                        if src != dst:
                            shutil.move(src, dst)
                            logging.info(f"Renamed {file} to Interaction_energy.csv")
        
        return True
    except Exception as e:
        logging.error(f"Error in post-processing: {str(e)}")
        logging.error(traceback.format_exc())
        return False
                
def find_tar_gz(index, pdbid):
    """Search for {index}.{pdbid}.tar.gz in remote directories"""
    for remote_dir in REMOTE_DIRS:
        # Use wildcard pattern to match any valid path
        check_command = f"sshpass -p {PASSWORD} ssh {ADA_USER} 'ls {remote_dir}/{index}.{pdbid}.tar.gz 2>/dev/null'"
        result = subprocess.run(check_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        if result.returncode == 0:
            return os.path.join(remote_dir, f"{index}.{pdbid}.tar.gz")
    return None

def copy_tar_gz(index, pdbid, tar_gz_path):
    """Copy .tar.gz from remote to local"""
    local_path = os.path.join(LOCAL_DIR, f"{index}.{pdbid}.tar.gz")
    subprocess.run(f"sshpass -p {PASSWORD} scp {ADA_USER}:{tar_gz_path} {local_path}", 
                   shell=True, check=True)
    print(f"Copied {tar_gz_path} to {local_path}")

def untar_file(index, pdbid):
    """Extract to copyingHES for processing with support for specific unusual nested structures"""
    tar_gz_path = os.path.join(LOCAL_DIR, f"{index}.{pdbid}.tar.gz")
    final_dir = os.path.join(LOCAL_DIR, f"{index}.{pdbid}")
    temp_extract_dir = os.path.join(LOCAL_DIR, f"temp_extract_{index}.{pdbid}")
    
    # Clean up existing directories
    for dir_path in [final_dir, temp_extract_dir]:
        if os.path.exists(dir_path):
            shutil.rmtree(dir_path)
    
    # Create extraction directories
    os.makedirs(final_dir, exist_ok=True)
    os.makedirs(temp_extract_dir, exist_ok=True)
    
    # Extract to temporary directory first (without stripping components)
    logging.info(f"Extracting {tar_gz_path} to temporary directory")
    subprocess.run(
        f"tar -xf {tar_gz_path} -C {temp_extract_dir}", 
        shell=True, 
        check=True
    )
    
    # Check for specific known structures
    known_paths = [
        # Structure 1: For archita's files
        os.path.join(temp_extract_dir, "share6", "archita.b", "HigherEnergy_OutputFiles", f"{index}.{pdbid}"),
        # Structure 2: For tanmay's files
        os.path.join(temp_extract_dir, "share6", "tanmay.goel", "higher_energy", f"{index}.{pdbid}"),
        # Structure 3: For sabarno's files
        os.path.join(temp_extract_dir, "share6", "sabarno.baral", "HigherEnergy_OutputFiles", f"{index}.{pdbid}"),
        # Structure 4: For khabeer's files
        os.path.join(temp_extract_dir, "share6" ,"khabeer.azhar", "higher_energy_output", f"{index}.{pdbid}"),
    ]
    
    content_dir = None
    
    # Check if any of the known structures exists
    for path in known_paths:
        if os.path.exists(path) and os.path.isdir(path):
            content_dir = path
            logging.info(f"Found known directory structure: {content_dir}")
            break
    
    # If none of the specific paths matched, try recursive search
    if not content_dir:
        logging.info("Known structures not found, attempting recursive search...")
        
        def find_content_dir(start_dir, target_files=None):
            """Recursively find directory containing important files/subdirectories"""
            if target_files is None:
                # Look for these key files or directories
                target_files = ["Results", "2.complex_solvated.prmtop", "complex_solvated.prmtop"]
            
            # Check if any key files exist in current directory
            if any(os.path.exists(os.path.join(start_dir, target)) for target in target_files):
                return start_dir
            
            # Check if this directory matches the index.pdbid pattern
            if os.path.basename(start_dir) == f"{index}.{pdbid}":
                if any(os.path.exists(os.path.join(start_dir, target)) for target in target_files):
                    return start_dir
            
            # Check subdirectories
            for item in os.listdir(start_dir):
                item_path = os.path.join(start_dir, item)
                if os.path.isdir(item_path):
                    # Search recursively in this subdirectory
                    result = find_content_dir(item_path, target_files)
                    if result:
                        return result
            
            return None
        
        content_dir = find_content_dir(temp_extract_dir)
    
    if not content_dir:
        logging.error(f"Could not locate content directory in extracted files for {index}.{pdbid}")
        # Clean up
        shutil.rmtree(temp_extract_dir)
        return False
    
    logging.info(f"Found content in: {content_dir}")
    
    # Move all content from the identified directory to the final directory
    for item in os.listdir(content_dir):
        source_path = os.path.join(content_dir, item)
        dest_path = os.path.join(final_dir, item)
        
        # Handle existing files/directories
        if os.path.exists(dest_path):
            if os.path.isdir(dest_path):
                shutil.rmtree(dest_path)
            else:
                os.remove(dest_path)
        
        # Move the file/directory
        shutil.move(source_path, dest_path)
    
    # Clean up temporary extraction directory
    shutil.rmtree(temp_extract_dir)
    
    # Verify key directories/files exist
    expected_paths = [
        os.path.join(final_dir, "Results"),
        os.path.join(final_dir, "Results", "2.complex_solvated.prmtop"),
        os.path.join(final_dir, "complex_solvated.prmtop")
    ]
    
    # At least one of these paths should exist
    if not any(os.path.exists(path) for path in expected_paths):
        logging.warning(f"Expected directory structure not found after extraction for {index}.{pdbid}")
        # Continue anyway since we might have a different structure
    
    logging.info(f"Successfully extracted {tar_gz_path} to {final_dir}")
    return True
    
def ensure_hes5k_dir_exists():
    """Create the HES5K directory if it doesn't exist"""
    logging.info(f"Ensuring HES5K directory exists: {HES5K_DIR}")
    mkdir_cmd = f"sshpass -p {PASSWORD} ssh {ADA_USER} 'mkdir -p {HES5K_DIR}'"
    try:
        subprocess.run(mkdir_cmd, shell=True, check=True)
        logging.info(f"Confirmed HES5K directory: {HES5K_DIR}")
        return True
    except subprocess.CalledProcessError as e:
        logging.error(f"Failed to create HES5K directory: {str(e)}")
        return False


def create_and_transfer_tar(index, pdbid):
    """Package and transfer processed results to remote server"""
    try:
        # Define paths
        target_dir = os.path.join(TRANSFER_DIR, f"{index}.{pdbid}")
        tar_path = os.path.join(TRANSFER_DIR, f"{index}.{pdbid}.tar.gz")
        
        # Check if target directory exists
        if not os.path.exists(target_dir):
            logging.error(f"Target directory not found: {target_dir}")
            return False
        
        # Post-process directory structure
        process_final_structure(index, pdbid)
            
        # Create tar archive
        logging.info(f"Creating tar archive: {tar_path}")
        subprocess.run(f"tar -czf {tar_path} -C {TRANSFER_DIR} {index}.{pdbid}", 
                     shell=True, check=True)
        logging.info(f"Successfully created tar archive: {tar_path}")
        
        # Ensure HES5K directory exists
        if not ensure_hes5k_dir_exists():
            logging.error("Failed to ensure HES5K directory exists")
            return False
        
        # Transfer to HES5K directory
        remote_path = f"{HES5K_DIR}"
        logging.info(f"Transferring tar archive to {ADA_USER}:{remote_path}/")
        subprocess.run(
            f"sshpass -p {PASSWORD} scp {tar_path} {ADA_USER}:{remote_path}/",
            shell=True, 
            check=True
        )
        logging.info(f"Successfully transferred {tar_path} to remote server")
        
        # Clean up local files
        os.remove(tar_path)
        logging.info(f"Cleaned up local tar file: {tar_path}")
        
        return True
    except Exception as e:
        logging.error(f"Error in create_and_transfer_tar: {str(e)}")
        logging.error(traceback.format_exc())
        return False
        
def process_index_pdbid(index, pdbid):
    """Full processing pipeline for a single index.pdbid"""
    try:
        logging.info(f"Starting full processing pipeline for {index}.{pdbid}")
        
        # Paths for cleanup
        local_dir = os.path.join(LOCAL_DIR, f"{index}.{pdbid}")
        local_tar = os.path.join(LOCAL_DIR, f"{index}.{pdbid}.tar.gz")
        transfer_dir = os.path.join(TRANSFER_DIR, f"{index}.{pdbid}")
        
        # 1. Find and download tar.gz from remote
        tar_path = find_tar_gz(index, pdbid)
        if not tar_path:
            logging.error(f"Could not find tar.gz for {index}.{pdbid} in any remote directory")
            return False
        
        logging.info(f"Found tar.gz at: {tar_path}")
        
        # 2. Copy tar.gz to local
        try:
            copy_tar_gz(index, pdbid, tar_path)
            logging.info(f"Successfully copied tar.gz for {index}.{pdbid}")
        except Exception as e:
            logging.error(f"Failed to copy tar.gz: {str(e)}")
            # Clean up if download failed
            if os.path.exists(local_tar):
                os.remove(local_tar)
                logging.info(f"Cleaned up incomplete tar file: {local_tar}")
            return False
        
        # 3. Extract tar.gz
        try:
            untar_file(index, pdbid)
            logging.info(f"Successfully extracted tar.gz for {index}.{pdbid}")
        except Exception as e:
            logging.error(f"Failed to extract tar.gz: {str(e)}")
            # Clean up downloaded tar and any partial extraction
            if os.path.exists(local_tar):
                os.remove(local_tar)
                logging.info(f"Cleaned up tar file after extraction failure: {local_tar}")
            if os.path.exists(local_dir):
                shutil.rmtree(local_dir)
                logging.info(f"Cleaned up partially extracted directory: {local_dir}")
            return False
        
        # 4. Process index directory
        if not process_index(f"{index} {pdbid}", LOCAL_DIR):
            logging.error(f"Failed to process {index}.{pdbid}")
            # Clean up downloaded and extracted files if processing fails
            cleanup_files(index, pdbid)
            return False
        logging.info(f"Successfully processed {index}.{pdbid}")
        
        # 5. Setup transfer structure
        if not setup_transfer_structure(index, pdbid):
            logging.error(f"Failed to setup transfer structure for {index}.{pdbid}")
            cleanup_files(index, pdbid)
            return False
        logging.info(f"Successfully setup transfer structure for {index}.{pdbid}")
        
        # 6. Create and transfer final tar
        if not create_and_transfer_tar(index, pdbid):
            logging.error(f"Failed to create and transfer tar for {index}.{pdbid}")
            cleanup_files(index, pdbid)
            return False
        logging.info(f"Successfully created and transferred tar for {index}.{pdbid}")
        
        # 7. Clean up local directories after successful processing
        cleanup_files(index, pdbid)
        return True
    except Exception as e:
        logging.error(f"Error in process_index_pdbid for {index}.{pdbid}: {str(e)}")
        logging.error(traceback.format_exc())
        # Clean up on any unexpected exception
        cleanup_files(index, pdbid)
        return False
        
def cleanup_files(index, pdbid):
    """Clean up all local files for a given index.pdbid - used after success or failure"""
    try:
        local_dir = os.path.join(LOCAL_DIR, f"{index}.{pdbid}")
        local_tar = os.path.join(LOCAL_DIR, f"{index}.{pdbid}.tar.gz")
        transfer_dir = os.path.join(TRANSFER_DIR, f"{index}.{pdbid}")
        
        if os.path.exists(local_dir):
            shutil.rmtree(local_dir)
            logging.info(f"Cleaned up extracted directory: {local_dir}")
        
        if os.path.exists(local_tar):
            os.remove(local_tar)
            logging.info(f"Cleaned up tar file: {local_tar}")
        
        if os.path.exists(transfer_dir):
            shutil.rmtree(transfer_dir)
            logging.info(f"Cleaned up transfer directory: {transfer_dir}")
            
        # Also clean temp directories if they exist
        temp_delta_dir = os.path.join(LOCAL_DIR, "temp_delta_files")
        if os.path.exists(temp_delta_dir):
            shutil.rmtree(temp_delta_dir)
            logging.info(f"Cleaned up temp delta directory: {temp_delta_dir}")
            
        logging.info(f"Completed cleanup for {index}.{pdbid}")
    except Exception as e:
        logging.warning(f"Warning: Error during cleanup for {index}.{pdbid}: {str(e)}")

def main():
    """Main entry point with concurrent processing"""
    log_start()
    
    # Create directories if they don't exist
    os.makedirs(LOCAL_DIR, exist_ok=True)
    os.makedirs(TRANSFER_DIR, exist_ok=True)
    
    # Check if input file exists
    if not os.path.exists(INPUT_TXT):
        logging.error(f"Input file not found: {INPUT_TXT}")
        return
    
    # Read and process input file
    with open(INPUT_TXT, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]
    
    logging.info(f"Found {len(lines)} entries to process in {INPUT_TXT}")
    
    # Process entries in parallel with thread pool
    with concurrent.futures.ThreadPoolExecutor(max_workers=5) as executor:
        futures = []
        results = {}
        
        for line in lines:
            try:
                parts = line.split()
                if len(parts) != 2:
                    logging.error(f"Invalid line format: {line}")
                    continue
                    
                index, pdbid = parts
                logging.info(f"Submitting {index}.{pdbid} for processing")
                future = executor.submit(process_index_pdbid, index, pdbid)
                futures.append(future)
                results[future] = f"{index}.{pdbid}"
            except Exception as e:
                logging.error(f"Error submitting job for line '{line}': {str(e)}")
        
        # Track progress and results
        success_count = 0
        fail_count = 0
        
        for future in concurrent.futures.as_completed(futures):
            pdb_id = results[future]
            try:
                success = future.result()
                if success:
                    logging.info(f"Successfully processed {pdb_id}")
                    success_count += 1
                else:
                    logging.error(f"Failed to process {pdb_id}")
                    fail_count += 1
            except Exception as e:
                logging.error(f"Processing failed for {pdb_id}: {str(e)}")
                logging.error(traceback.format_exc())
                fail_count += 1
    
    logging.info("=" * 80)
    logging.info(f"PROCESSING COMPLETE - {datetime.datetime.now()}")
    logging.info(f"Successful: {success_count}, Failed: {fail_count}, Total: {len(lines)}")
    logging.info("=" * 80)

if __name__ == "__main__":
    main()
