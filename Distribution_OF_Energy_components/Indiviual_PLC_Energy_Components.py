import os
import pandas as pd
import matplotlib.pyplot as plt
import logging
import numpy as np
import matplotlib as mpl
from matplotlib import font_manager

# Set up logging
logging.basicConfig(
    filename='plot_generation.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)

def plot_average_energy_components(csv_file):
    """Create the average energy components plot across all runs with all components properly displayed"""
    plt.ioff()
    
    # ========================================
    # FONT CONFIGURATION (Arial Regular + Bold)
    # ========================================
    try:
        arial_regular_path = "/home/user/miniconda3/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf"
        arial_bold_path = "/home/user/miniconda3/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arialbd.ttf"

        if os.path.exists(arial_regular_path) and os.path.exists(arial_bold_path):
            font_regular = font_manager.FontProperties(fname=arial_regular_path)
            font_bold = font_manager.FontProperties(fname=arial_bold_path)
        else:
            raise FileNotFoundError("Arial fonts not found at specified path.")
    except Exception as e:
        print(f"Font warning: {str(e)}. Using default font.")
        font_regular = None
        font_bold = None

    try:
        # ===== DATA PROCESSING AND PLOTTING START HERE =====

        df = pd.read_csv(csv_file)
        print(f"Initial data shape: {df.shape}")
        print(f"Initial columns: {df.columns}")
        
        # Clean column names - important to preserve special characters
        df.columns = [col.strip() for col in df.columns]
        print(f"\nCleaned columns: {df.columns.tolist()}")
        
        # Clean data
        df = df.dropna(how='all')
        
        # Convert Frame to numeric
        if 'Frame' in df.columns:
            df['Frame'] = pd.to_numeric(df['Frame'], errors='coerce')
            df = df[df['Frame'].notna()]
            df['Frame'] = df['Frame'].astype(int)
        else:
            raise ValueError("Frame column not found after cleaning")

        # Explicitly identify expected column names (with correct special characters)
        expected_columns = ['∆Evdw', '∆Eele', '∆Gpol', '∆Gnp', '∆GMM–PBSA']
        
        # Check for columns and print proper error message
        found_columns = []
        missing_columns = []
        for col in expected_columns:
            if col in df.columns:
                found_columns.append(col)
            else:
                missing_columns.append(col)
        
        print(f"Found energy columns: {found_columns}")
        if missing_columns:
            print(f"Warning: Could not find these columns: {missing_columns}")
            print("Available columns:", df.columns.tolist())
            
        energy_columns = found_columns
        
        if not energy_columns:
            raise ValueError("No energy component columns found. Check column names.")

        # Convert energy columns to numeric
        for col in energy_columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')
            if df[col].isna().any():
                print(f"Warning: {col} contains non-numeric values that were coerced to NaN")

        # Add Run column if missing
        if 'Run' not in df.columns:
            df['Run'] = 0

        # ========================================
        # PLOTTING SECTION
        # ========================================
        plt.figure(figsize=(12, 8), dpi=300)
        
        # Color mapping - explicitly use the actual column names from the dataframe
        colors = {
            '∆Evdw': '#000000',    # Black
            '∆Eele': '#1f77b4',    # Blue
            '∆Gpol': '#ff7f0e',    # Orange
            '∆Gnp': '#228B22',     # Green
            '∆GMM–PBSA': '#FF0000' # Red - note the en-dash character
        }
        
        # Plot each component
        for component in energy_columns:
            # Skip empty columns
            if df[component].isna().all():
                print(f"Skipping {component} - all values NaN")
                continue
            
            # Get average data
            avg_data = df.groupby('Frame')[component].mean().reset_index()
            
            # Create label with subscript
            clean_component = component.replace('∆', '').replace('Delta', '')
            if clean_component.startswith(('E', 'G')):
                if clean_component == 'GMM–PBSA':  # Handle the total energy differently
                    label = r'$\Delta G_{\mathrm{MM-PBSA}}$'
                else:
                    base = f"Δ{clean_component[0]}"
                    suffix = clean_component[1:]
                    label = f"${base}_{{{suffix}}}$"
            else:
                label = f"${component}$"
            
            # Use a thicker line for the total energy
            lw = 2.0 if component == '∆GMM–PBSA' else 1.5
            
            plt.plot(avg_data['Frame'], avg_data[component],
                     marker='.', 
                     markersize=6,
                     linewidth=lw,
                     color=colors.get(component, '#9467bd'), 
                     label=label)

        # Style adjustments
        plt.tick_params(axis='both', which='both', direction='in', length=8, width=1.0, labelsize=20)
        
        # Set labels with appropriate font
        plt.xlabel("Frame", fontsize=22, fontproperties=font_bold)
        plt.ylabel("Energy (kcal/mol)", fontsize=22, fontproperties=font_bold)
        
        # Save output
        output_file = os.path.splitext(csv_file)[0] + "_all_components.png"
        plt.savefig(output_file, bbox_inches='tight')
        
        # Also save as PDF for publication quality
        pdf_output = os.path.splitext(csv_file)[0] + "_all_components.pdf"
        plt.savefig(pdf_output, bbox_inches='tight')
        
        plt.close()
        print(f"\nSaved plots to {output_file} and {pdf_output}")
        return output_file

    except Exception as e:
        logging.error(f"Error: {str(e)}")
        print(f"Critical error: {str(e)}")
        import traceback
        traceback.print_exc()
        return None

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Generate energy components plot")
    parser.add_argument("csv_file", help="Path to the CSV file")
    args = parser.parse_args()
    plot_average_energy_components(args.csv_file)
