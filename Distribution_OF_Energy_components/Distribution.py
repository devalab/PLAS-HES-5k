import os
import tarfile
import pandas as pd
import numpy as np
import logging
import tempfile
import re
import glob
from tqdm import tqdm

# Setup logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[
        logging.FileHandler("processing.log"),
        logging.StreamHandler()
    ]
)

def process_interaction_energy(file_path):
    """Process a single Interaction_energy.csv file and return the data needed for the CSV"""
    try:
        # Read the CSV file
        # First try to determine the structure
        with open(file_path, 'r') as f:
            lines = f.readlines()
        
        # Look for "DELTA Energy Terms" section - it's the last section of the file
        delta_start_idx = None
        for i, line in enumerate(lines):
            if "DELTA Energy Terms" in line:
                delta_start_idx = i
                break
        
        if delta_start_idx is None:
            logging.error(f"Could not find 'DELTA Energy Terms' section in {file_path}")
            return [], [], file_path, False
        
        # Find the header line (contains Frame #, BOND, ANGLE, etc.)
        header_idx = delta_start_idx + 1
        if header_idx >= len(lines):
            logging.error(f"File structure issue in {file_path}")
            return [], [], file_path, False
        
        # Find where the data starts (after the header)
        data_start_idx = header_idx + 1
        if data_start_idx >= len(lines):
            logging.error(f"File structure issue in {file_path}")
            return [], [], file_path, False
        
        # Read the data into a DataFrame
        df = pd.read_csv(file_path, skiprows=132)
        
        # Check if DataFrame has at least 40 rows (for 10 lower + 30 upper)
        if len(df) < 40:
            logging.warning(f"Not enough data points in {file_path}")
            return [], [], file_path, False
        
        # Extract columns - Frame # is col 0, VDWAALS is col 7, EEL is col 8, 
        # EPB is col 11, ENPOLAR is col 12, DELTA TOTAL is the last column
        frame_column = df['Frame #']
        col7 = df['VDWAALS']
        col8 = df['EEL']
        col11 = df['EPB']
        col12 = df['ENPOLAR']
        last_column = df['DELTA TOTAL']
        
        # Convert to numeric, drop NaN values
        frames = pd.to_numeric(frame_column, errors='coerce').dropna()
        data = pd.to_numeric(last_column, errors='coerce').dropna()
        col7_data = pd.to_numeric(col7, errors='coerce').dropna()
        col8_data = pd.to_numeric(col8, errors='coerce').dropna()
        col11_data = pd.to_numeric(col11, errors='coerce').dropna()
        col12_data = pd.to_numeric(col12, errors='coerce').dropna()
        
        # Make sure we have matching length data
        min_length = min(len(frames), len(data), len(col7_data), len(col8_data), len(col11_data), len(col12_data))
        frames = frames[:min_length]
        data = data[:min_length]
        col7_data = col7_data[:min_length]
        col8_data = col8_data[:min_length]
        col11_data = col11_data[:min_length]
        col12_data = col12_data[:min_length]
        
        # Create DataFrame with all needed columns
        full_df = pd.DataFrame({
            'frame': frames,
            'col7': col7_data,
            'col8': col8_data,
            'col11': col11_data,
            'col12': col12_data,
            'value': data
        })
        
        # Remove rows with NaN values and filter out values greater than 10
        full_df = full_df.dropna()
        full_df = full_df[full_df['value'] < 10]
        
        # Check if we have enough data after filtering
        if len(full_df) < 40:
            logging.warning(f"Not enough data points after filtering in {file_path}")
            return [], [], file_path, False
        
        # Sort by frame number
        full_df = full_df.sort_values(by='frame').reset_index(drop=True)
        
        # Take first 10 frames as lower threshold
        lower_frames = full_df.iloc[:10].to_dict('records')
        for frame in lower_frames:
            frame['segment'] = 'lower'
        
        # Take next 30 frames as upper threshold
        upper_frames = full_df.iloc[10:40].to_dict('records')
        for frame in upper_frames:
            frame['segment'] = 'upper'
        
        return lower_frames, upper_frames, file_path, True
        
    except Exception as e:
        logging.error(f"Error processing {file_path}: {e}")
        return [], [], file_path, False

def process_tar_file(tar_path, experimental_pdbs):
    """Process a single tar.gz file and extract the required data"""
    file_basename = os.path.basename(tar_path)
    # Extract PDB ID from the filename (format: INDEX.PDB_ID.tar.gz)
    match = re.search(r'(\d+)\.(\w+)\.tar\.gz', file_basename)
    
    if match:
        index = match.group(1)
        pdb_id = match.group(2).lower()
        pdb_name = f"{index}.{pdb_id}"
    else:
        pdb_name = file_basename.replace('.tar.gz', '')
        pdb_id = pdb_name.lower()
    
    logging.info(f"Processing {pdb_name} (PDB ID: {pdb_id})")
    
    # Determine if this is an experimental PDB
    is_experimental = pdb_id in experimental_pdbs
    
    # Create temporary directory for extraction
    with tempfile.TemporaryDirectory() as temp_dir:
        try:
            # Extract tar.gz file
            with tarfile.open(tar_path) as tar:
                tar.extractall(path=temp_dir)
            
            # Find HES directory - first check if it's in a subdirectory with the name pattern
            base_dir_name = pdb_name  # e.g., "4935.2h4k"
            possible_parent_dir = os.path.join(temp_dir, base_dir_name)
            
            if os.path.exists(possible_parent_dir) and os.path.isdir(possible_parent_dir):
                # The tar file extracted to INDEX.PDB_ID/HES/
                hes_dir = os.path.join(possible_parent_dir, 'HES')
            else:
                # Try the original expected path
                hes_dir = os.path.join(temp_dir, 'HES')
                
            if not os.path.exists(hes_dir):
                logging.warning(f"HES directory not found in {tar_path}")
                return None
            
            all_results = []
            
            # Process each run directory
            for run_idx in range(5):  # run0 to run4
                run_dir = os.path.join(hes_dir, f'run{run_idx}')
                if not os.path.exists(run_dir):
                    logging.warning(f"Run directory {run_dir} not found in {tar_path}")
                    continue
                
                # Find and process Interaction_energy.csv
                interaction_file = os.path.join(run_dir, 'Interaction_energy.csv')
                if not os.path.exists(interaction_file):
                    logging.warning(f"Interaction_energy.csv not found in {run_dir}")
                    continue
                
                lower_frames, upper_frames, file_path, success = process_interaction_energy(interaction_file)
                
                if not success:
                    logging.warning(f"Failed to process {interaction_file}")
                    continue
                
                # Add results to our data collection
                results_data = []
                
                # Process lower frames
                for info in lower_frames:
                    results_data.append({
                        'file_name': pdb_name,
                        'run': f'run{run_idx}',
                        'frame': int(info['frame']),
                        'col7_value': info['col7'],
                        'col8_value': info['col8'],
                        'col11_value': info['col11'],
                        'col12_value': info['col12'],
                        'bfe_value': info['value'],
                        'threshold': info['segment']
                    })
                
                # Process upper frames
                for info in upper_frames:
                    results_data.append({
                        'file_name': pdb_name,
                        'run': f'run{run_idx}',
                        'frame': int(info['frame']),
                        'col7_value': info['col7'],
                        'col8_value': info['col8'],
                        'col11_value': info['col11'],
                        'col12_value': info['col12'],
                        'bfe_value': info['value'],
                        'threshold': info['segment']
                    })
                
                all_results.extend(results_data)
            
            # Return the results along with classification
            return {
                'results': all_results,
                'is_experimental': is_experimental
            }
            
        except Exception as e:
            logging.error(f"Error processing {tar_path}: {e}")
            return None

def load_experimental_pdbs(csv_path):
    """Load experimental PDSs from CSV file"""
    experimental_pdbs = set()
    try:
        if os.path.exists(csv_path):
            exp_df = pd.read_csv(csv_path)
            # Extract PDB_ID column and convert to lowercase for case-insensitive matching
            if 'PDB_ID' in exp_df.columns:
                experimental_pdbs = {pdb_id.strip().lower() for pdb_id in exp_df['PDB_ID']}
            logging.info(f"Loaded {len(experimental_pdbs)} experimental PDBs from {csv_path}")
        else:
            # If the file doesn't exist, create it with the example data
            example_data = {
                'INDEX': [5150, 7404],
                'PDB_ID': ['1j96', '1ttj']
            }
            pd.DataFrame(example_data).to_csv(csv_path, index=False)
            experimental_pdbs = {'1j96', '1ttj'}
            logging.info(f"Created experimental PDBs file with example data: {csv_path}")
    except Exception as e:
        logging.error(f"Error reading experimental PDSs list: {e}")
    
    return experimental_pdbs

def main():
    base_dir = "path_to_csv_directory"
    experimental_list_path = "experimental_pdb.csv"
    
    # Create output directories
    output_dirs = {
        'Experimental': os.path.join(os.getcwd(), 'Experimental'),
        'Non-Experimental': os.path.join(os.getcwd(), 'Non-Experimental'),
        'All5000': os.path.join(os.getcwd(), 'All5000')
    }
    
    for directory in output_dirs.values():
        os.makedirs(directory, exist_ok=True)
    
    # Initialize CSV files
    csv_files = {}
    for category in output_dirs:
        csv_files[f"{category}_lower"] = os.path.join(output_dirs[category], "lowerthreshold_energy.csv")
        csv_files[f"{category}_upper"] = os.path.join(output_dirs[category], "upperthreshold_energy.csv")
        
        # Create CSV files with headers
        header = ['file_name', 'run', 'frame', 'col7_value', 'col8_value', 'col11_value', 'col12_value', 'bfe_value', 'threshold']
        for csv_file in csv_files.values():
            if not os.path.exists(csv_file):
                pd.DataFrame(columns=header).to_csv(csv_file, index=False)
    
    # Read experimental PDSs list from CSV file
    experimental_pdbs = load_experimental_pdbs(experimental_list_path)
    
    # Get list of all tar.gz files - use specific pattern matching
    tar_files_pattern = os.path.join(base_dir, "*.tar.gz")
    tar_files = glob.glob(tar_files_pattern)
    
    # Check if base_dir exists
    if not os.path.exists(base_dir):
        logging.warning(f"Base directory {base_dir} does not exist")
        
        # Check the current directory as fallback
        current_dir_files = glob.glob("*.tar.gz")
        if current_dir_files:
            tar_files = current_dir_files
            logging.info(f"Found {len(tar_files)} tar.gz files in current directory")
        else:
            logging.error("No tar.gz files found. Please check file paths.")
            # Use the example files as a last resort
            example_files = [
                "5150.1j96.tar.gz",
                "6269.2c5b.tar.gz",
                "7404.1ttj.tar.gz",
                "8523.4w4w.tar.gz"
            ]
            tar_files = [file for file in example_files if os.path.exists(file)]
            if tar_files:
                logging.info(f"Using {len(tar_files)} example files")
            else:
                logging.error("No example files found either. Exiting.")
                return
    
    logging.info(f"Found {len(tar_files)} tar.gz files")
    
    # Process each tar file
    for tar_file in tqdm(tar_files, desc="Processing tar files"):
        result = process_tar_file(tar_file, experimental_pdbs)
        
        if not result or not result['results']:
            continue
        
        # Separate lower and upper threshold data
        lower_data = [item for item in result['results'] if item['threshold'] == 'lower']
        upper_data = [item for item in result['results'] if item['threshold'] == 'upper']
        
        # Determine which categories to write to
        categories = ['All5000']
        if result['is_experimental']:
            categories.append('Experimental')
        else:
            categories.append('Non-Experimental')
        
        # Write to appropriate CSV files
        for category in categories:
            if lower_data:
                lower_df = pd.DataFrame(lower_data)
                lower_df.to_csv(csv_files[f"{category}_lower"], mode='a', header=False, index=False)
            
            if upper_data:
                upper_df = pd.DataFrame(upper_data)
                upper_df.to_csv(csv_files[f"{category}_upper"], mode='a', header=False, index=False)
    
    logging.info("Processing complete")

if __name__ == "__main__":
    main()
