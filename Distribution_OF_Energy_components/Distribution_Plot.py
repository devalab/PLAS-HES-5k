import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os
from matplotlib.font_manager import FontProperties
plt.rcParams.update({'font.size': 22})
font_path ="/home/user/miniconda3/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arial.ttf"
arial_font = FontProperties(fname=font_path)
bold_path= "/home/user/miniconda3/lib/python3.12/site-packages/matplotlib/mpl-data/fonts/ttf/arialbd.ttf"
bold= FontProperties(fname=bold_path)
categories = ["Experimental", "Non-Experimental", "All5000"]
category_titles = ["Experimental", "Non-Experimental", "Entire Datapoints"]

columns_to_plot = ["bfe_value", "col7_value", "col8_value", "col11_value", "col12_value"]
column_labels = {
    "bfe_value": "a)",
    "col7_value": "b)",
    "col8_value": "c)",
    "col11_value": "d)",
    "col12_value": "e)"
}

# Define specific x and y ranges for each column
x_ranges = {
    "bfe_value": (-100, 0),        # Range of 100
    "col7_value": (-100, 0),       # Range of 100
    "col8_value": (-50, 50),       # Range of 100
    "col11_value": (-50, 50),      # Range of 100
    "col12_value": (-80, 20)       # Range of 100
}

y_ranges = {
    "bfe_value": (0, 0.065),       # 0 to 0.065
    "col7_value": (0, 0.065),      # 0 to 0.065
    "col8_value": (0, 0.15),       # 0 to 0.15
    "col11_value": (0, 0.15),      # 0 to 0.15
    "col12_value": (0, 0.3)        # 0 to 0.3
}

def load_csvs(base_dir):
    data = {}
    for category in categories:
        lower_path = os.path.join(base_dir, category, "lowerthreshold_energy.csv")
        upper_path = os.path.join(base_dir, category, "upperthreshold_energy.csv")
        try:
            df_lower = pd.read_csv(lower_path)
            df_lower["threshold"] = "lower"
            df_upper = pd.read_csv(upper_path)
            df_upper["threshold"] = "upper"
            data[category] = pd.concat([df_lower, df_upper])
        except Exception as e:
            print(f"Error loading {category}: {e}")
            data[category] = pd.DataFrame()
    return data

# Simple custom KDE function using numpy instead of scipy
def simple_kde(data, x_grid, bandwidth=None):
    """
    Compute a simple kernel density estimate using Gaussian kernels
    """
    data = np.asarray(data)
    n = len(data)
    
    if n == 0:
        return np.zeros_like(x_grid)
    
    # Use Scott's rule for bandwidth if not specified
    if bandwidth is None:
        bandwidth = 1.06 * np.std(data) * n**(-1/5)
        # Ensure minimum bandwidth
        bandwidth = max(bandwidth, 0.5)
    
    # Compute the KDE
    kde = np.zeros_like(x_grid, dtype=float)
    for x in data:
        # Add a Gaussian kernel for each data point
        kernel = np.exp(-0.5 * ((x_grid - x) / bandwidth)**2) / (bandwidth * np.sqrt(2 * np.pi))
        kde += kernel
    
    # Normalize
    kde /= n
    
    return kde

def plot_matrix(data_dict, output_file="combined_threshold_distribution_scaled.png"):
    # Set the default to Arial globally
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.sans-serif'] = ['Arial']
    
    fig, axes = plt.subplots(nrows=len(columns_to_plot), ncols=len(categories), figsize=(19, 16))
    # Changed threshold_colors to match the legend labels in your example image
    threshold_colors = {"lower": "skyblue", "upper": "lightcoral"}
    edge_colors = {"lower": "blue", "upper": "red"}
    line_colors = {"lower": "blue", "upper": "red"}
    bins = 50  # Fixed number of bins for all histograms
    
    for col_idx, column in enumerate(columns_to_plot):
        for cat_idx, category in enumerate(categories):
            ax = axes[col_idx, cat_idx]
            df = data_dict.get(category, pd.DataFrame())
            if df.empty or column not in df.columns:
                ax.set_visible(False)
                continue

            # Set predefined x and y ranges instead of calculating from data
            x_min, x_max = x_ranges[column]
            y_min, y_max = y_ranges[column]
            
            # Calculate bin width based on the range and number of bins
            bin_width = (x_max - x_min) / bins
            
            # Set x and y limits
            ax.set_xlim(x_min, x_max)
            ax.set_ylim(y_min, y_max)
            
            # Set nicer round-number ticks based on the column type
            if column in ["bfe_value", "col7_value"]:  # Range -100 to 0
                ax.set_xticks([-100, -75, -50, -25, 0])
            elif column in ["col8_value", "col11_value"]:  # Range -50 to 50
                ax.set_xticks([-50, -25, 0, 25, 50])
            elif column == "col12_value":  # Range -80 to 20
                ax.set_xticks([-80, -60, -40, -20, 0, 20])
            
            # Set nice round y-ticks based on the y-range
            if y_max <= 0.1:
                ax.set_yticks([ 0.02, 0.04, 0.06])
            elif y_max <= 0.2:
                ax.set_yticks([ 0.05, 0.10, 0.15])
            else:  # For larger ranges (col12_value)
                ax.set_yticks([ 0.1, 0.2, 0.3])
            
            # Format the tick labels to avoid scientific notation and reduce decimals
            ax.xaxis.set_major_formatter(plt.FormatStrFormatter('%d'))
            ax.yaxis.set_major_formatter(plt.FormatStrFormatter('%.2f'))
            
            # Remove grid
            ax.grid(False)

            # Create common x grid for KDE
            x_kde = np.linspace(x_min, x_max, 1000)
            
            for threshold in ["lower", "upper"]:
                subset = pd.to_numeric(df[df["threshold"] == threshold][column], errors="coerce").dropna()
                if len(subset) == 0:
                    continue
                
                # Use fixed range for histogram bins with density=True
                ax.hist(subset, bins=bins, alpha=0.7, density=True,
                        color=threshold_colors[threshold], 
                        edgecolor=edge_colors[threshold],
                        linewidth=0.8,
                        range=(x_min, x_max))
                        
                # Add KDE curve using our simple function
                if len(subset) > 1:  # Need at least 2 points for KDE
                    # Adjust bandwidth based on bin width
                    bandwidth = bin_width * 2
                    
                    y_kde = simple_kde(subset, x_kde, bandwidth)
                    ax.plot(x_kde, y_kde, color=line_colors[threshold], linewidth=2)
            
            # Only set title for top row
            if col_idx == 0:
                ax.set_title(category_titles[cat_idx], fontsize=7, fontproperties=arial_font)
                
            # Add a), b), c), d), e) labels to the left column plots
            if cat_idx == 0:
                # Position the label closer to the y-axis
                ax.text(-0.3, 0.5, column_labels[column], 
                        fontsize=27, fontweight='bold', ha='right', va='center',
                        transform=ax.transAxes, fontproperties=bold)  # Use axis coordinates
            
            # Set tick parameters for both axes - increased length to 8
            ax.tick_params(axis='both', direction='in', which='both', length=8, pad=10)
            
            # Only show y-axis labels and tick values on the leftmost column
            if cat_idx > 0:
                # Hide y-axis ticks and labels for middle and right columns
                ax.set_yticklabels([])
            
            # Make tick labels smaller for better appearance
            ax.tick_params(axis='both', labelsize=25)
            for label in ax.get_xticklabels() + ax.get_yticklabels():
                label.set_fontproperties(arial_font)
            
            # Only set x-axis label for the last row, middle column
            if col_idx == len(columns_to_plot) - 1 and cat_idx == 1:
                # Increased padding for x-axis label (moved down from ticks)
                ax.xaxis.set_label_coords(0.5, -0.3)  # Adjust y-coordinate to move label further down
                ax.set_xlabel("Energy (kcal/mol)", fontsize=22, fontweight='bold', fontproperties=bold)
            else:
                ax.set_xlabel("")

    # Use increased spacing between subplots
    plt.subplots_adjust(hspace=0.6, wspace=0.4, bottom=0.13)  # Increased bottom margin
    
    # Use tight_layout with adjusted parameters for more space between plots
    plt.tight_layout(rect=[0, 0.05, 1, 0.98], pad=1, h_pad=0.9, w_pad=0.7)
    
    fig.text(0.03, 0.5, "Probability Density", va='center', rotation='vertical',
         fontsize=22, fontweight='bold', fontproperties=bold)
    
    # Save with high dpi for publication quality
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()

if __name__ == "__main__":
    base_dir = "."  # Set path to where your data folders are
    data = load_csvs(base_dir)
    plot_matrix(data)
