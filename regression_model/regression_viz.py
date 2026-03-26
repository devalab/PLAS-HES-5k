import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from scipy import stats
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
import gc

def smart_sample(y_true, y_pred, max_points=50000):
    """Intelligently sample data if it's too large"""
    if len(y_true) <= max_points:
        return y_true, y_pred
        
    # Stratified sampling to preserve distribution
    n_bins = 100
    bins = np.linspace(y_true.min(), y_true.max(), n_bins + 1)
    bin_indices = np.digitize(y_true, bins) - 1
    bin_indices = np.clip(bin_indices, 0, n_bins - 1)
    
    sampled_idx = []
    points_per_bin = max_points // n_bins
    
    for i in range(n_bins):
        mask = bin_indices == i
        bin_size = mask.sum()
        if bin_size > 0:
            n_samples = min(points_per_bin, bin_size)
            bin_idx = np.where(mask)[0]
            sampled_idx.extend(np.random.choice(bin_idx, n_samples, replace=False))
    
    sampled_idx = np.array(sampled_idx)
    return y_true[sampled_idx], y_pred[sampled_idx]


def calculate_metrics(y_true, y_pred):
    """Calculate all metrics at once"""
    r2 = r2_score(y_true, y_pred)
    mae = mean_absolute_error(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    pearson_r = np.corrcoef(y_true, y_pred)[0, 1]
    return r2, mae, rmse, pearson_r


def plot_hexbin_regression(y_true, y_pred, title="Hexbin Regression Plot", save_path=None, max_points=50000):
    """Memory-efficient hexbin plot"""
    # Sample if necessary
    y_true_sampled, y_pred_sampled = smart_sample(y_true, y_pred, max_points)
    
    # Calculate metrics on full dataset
    r2, mae, rmse, pearson_r = calculate_metrics(y_true, y_pred)
    
    plt.figure(figsize=(10, 8))
    
    # Create hexbin plot with reasonable gridsize
    hexbin = plt.hexbin(y_true_sampled, y_pred_sampled, gridsize=40, cmap='Blues', 
                        mincnt=1, reduce_C_function=np.sum)
    
    # Add perfect prediction line
    min_val = min(y_true_sampled.min(), y_pred_sampled.min())
    max_val = max(y_true_sampled.max(), y_pred_sampled.max())
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', lw=2)
    
    # Add text box with metrics
    textstr = f'R² = {r2:.3f}\nPearson r = {pearson_r:.3f}\nMAE = {mae:.3f}\nRMSE = {rmse:.3f}'
    plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, 
             fontsize=10, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.colorbar(hexbin, label='Count')
    plt.xlabel('True Values', fontsize=12)
    plt.ylabel('Predicted Values', fontsize=12)
    plt.title(title, fontsize=14)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    gc.collect()
    print(f"Saved: {save_path}")


def plot_residual_analysis(y_true, y_pred, title="Residual Analysis", save_path=None, max_points=50000):
    """Memory-efficient residual analysis"""
    residuals = y_pred - y_true
    
    # Sample for scatter plots but use full data for histogram and Q-Q
    y_true_sampled, y_pred_sampled = smart_sample(y_true, y_pred, max_points)
    residuals_sampled = y_pred_sampled - y_true_sampled
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # 1. Residual vs Predicted
    axes[0, 0].scatter(y_pred_sampled, residuals_sampled, alpha=0.5, s=5)
    axes[0, 0].axhline(y=0, color='r', linestyle='--', lw=2)
    axes[0, 0].set_xlabel('Predicted Values', fontsize=10)
    axes[0, 0].set_ylabel('Residuals', fontsize=10)
    axes[0, 0].set_title('Residuals vs Predicted', fontsize=11)
    axes[0, 0].grid(alpha=0.3)
    
    # 2. Residual vs True
    axes[0, 1].scatter(y_true_sampled, residuals_sampled, alpha=0.5, s=5)
    axes[0, 1].axhline(y=0, color='r', linestyle='--', lw=2)
    axes[0, 1].set_xlabel('True Values', fontsize=10)
    axes[0, 1].set_ylabel('Residuals', fontsize=10)
    axes[0, 1].set_title('Residuals vs True', fontsize=11)
    axes[0, 1].grid(alpha=0.3)
    
    # 3. Histogram of residuals (use full data if possible)
    sample_for_hist = min(100000, len(residuals))
    if len(residuals) > sample_for_hist:
        residuals_hist = np.random.choice(residuals, sample_for_hist, replace=False)
    else:
        residuals_hist = residuals
    
    axes[1, 0].hist(residuals_hist, bins=50, edgecolor='black', alpha=0.7)
    axes[1, 0].axvline(x=0, color='r', linestyle='--', lw=2)
    axes[1, 0].set_xlabel('Residuals', fontsize=10)
    axes[1, 0].set_ylabel('Frequency', fontsize=10)
    axes[1, 0].set_title(f'Residual Distribution (μ={residuals.mean():.3f}, σ={residuals.std():.3f})', fontsize=11)
    axes[1, 0].grid(alpha=0.3)
    
    # 4. Q-Q plot (sample for memory efficiency)
    sample_for_qq = min(10000, len(residuals))
    if len(residuals) > sample_for_qq:
        residuals_qq = np.random.choice(residuals, sample_for_qq, replace=False)
    else:
        residuals_qq = residuals
    
    stats.probplot(residuals_qq, dist="norm", plot=axes[1, 1])
    axes[1, 1].set_title('Q-Q Plot', fontsize=11)
    axes[1, 1].grid(alpha=0.3)
    
    plt.suptitle(title, fontsize=14, y=1.00)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    gc.collect()
    print(f"Saved: {save_path}")


def plot_error_bands(y_true, y_pred, title="Prediction with Error Bands", save_path=None, max_points=10000):
    """Memory-efficient error band plot"""
    # Sample for this visualization
    y_true_sampled, y_pred_sampled = smart_sample(y_true, y_pred, max_points)
    
    plt.figure(figsize=(12, 8))
    
    # Sort by true values
    sorted_idx = np.argsort(y_true_sampled)
    y_true_sorted = y_true_sampled[sorted_idx]
    y_pred_sorted = y_pred_sampled[sorted_idx]
    
    # Calculate rolling statistics
    window = max(50, len(y_true_sorted) // 100)
    rolling_mean = np.convolve(y_pred_sorted, np.ones(window)/window, mode='same')
    residuals = y_pred_sorted - y_true_sorted
    rolling_std = np.array([np.std(residuals[max(0, i-window//2):min(len(residuals), i+window//2)]) 
                            for i in range(len(residuals))])
    
    # Plot
    plt.scatter(y_true_sorted, y_pred_sorted, alpha=0.3, s=3, label='Predictions')
    plt.plot(y_true_sorted, y_true_sorted, 'r--', lw=2)
    plt.plot(y_true_sorted, rolling_mean, 'g-', lw=2, label='Rolling Mean')
    plt.fill_between(y_true_sorted, rolling_mean - 2*rolling_std, rolling_mean + 2*rolling_std, 
                     alpha=0.2, color='green', label='±2σ band')
    
    plt.xlabel('True Values', fontsize=12)
    plt.ylabel('Predicted Values', fontsize=12)
    plt.title(title, fontsize=14)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    gc.collect()
    print(f"Saved: {save_path}")


def plot_binned_performance(y_true, y_pred, n_bins=20, title="Binned Performance Analysis", save_path=None):
    """Analyze performance across different value ranges - memory efficient"""
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    # Create bins
    bins = np.linspace(y_true.min(), y_true.max(), n_bins + 1)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    bin_indices = np.digitize(y_true, bins) - 1
    bin_indices = np.clip(bin_indices, 0, n_bins - 1)
    
    # Calculate metrics per bin
    mae_per_bin = []
    rmse_per_bin = []
    r2_per_bin = []
    counts_per_bin = []
    
    for i in range(n_bins):
        mask = bin_indices == i
        if mask.sum() > 1:
            mae_per_bin.append(mean_absolute_error(y_true[mask], y_pred[mask]))
            rmse_per_bin.append(np.sqrt(mean_squared_error(y_true[mask], y_pred[mask])))
            r2_per_bin.append(r2_score(y_true[mask], y_pred[mask]))
            counts_per_bin.append(mask.sum())
        else:
            mae_per_bin.append(np.nan)
            rmse_per_bin.append(np.nan)
            r2_per_bin.append(np.nan)
            counts_per_bin.append(0)
    
    # Plot 1: MAE per bin
    axes[0, 0].bar(bin_centers, mae_per_bin, width=(bins[1]-bins[0])*0.8, alpha=0.7, edgecolor='black')
    axes[0, 0].set_xlabel('True Value Range', fontsize=10)
    axes[0, 0].set_ylabel('MAE', fontsize=10)
    axes[0, 0].set_title('Mean Absolute Error per Bin', fontsize=11)
    axes[0, 0].grid(alpha=0.3)
    
    # Plot 2: RMSE per bin
    axes[0, 1].bar(bin_centers, rmse_per_bin, width=(bins[1]-bins[0])*0.8, alpha=0.7, 
                   color='orange', edgecolor='black')
    axes[0, 1].set_xlabel('True Value Range', fontsize=10)
    axes[0, 1].set_ylabel('RMSE', fontsize=10)
    axes[0, 1].set_title('Root Mean Squared Error per Bin', fontsize=11)
    axes[0, 1].grid(alpha=0.3)
    
    # Plot 3: R² per bin
    axes[1, 0].bar(bin_centers, r2_per_bin, width=(bins[1]-bins[0])*0.8, alpha=0.7, 
                   color='green', edgecolor='black')
    axes[1, 0].set_xlabel('True Value Range', fontsize=10)
    axes[1, 0].set_ylabel('R²', fontsize=10)
    axes[1, 0].set_title('R² Score per Bin', fontsize=11)
    axes[1, 0].axhline(y=0, color='red', linestyle='--', alpha=0.5)
    axes[1, 0].grid(alpha=0.3)
    
    # Plot 4: Sample count per bin
    axes[1, 1].bar(bin_centers, counts_per_bin, width=(bins[1]-bins[0])*0.8, alpha=0.7, 
                   color='purple', edgecolor='black')
    axes[1, 1].set_xlabel('True Value Range', fontsize=10)
    axes[1, 1].set_ylabel('Count', fontsize=10)
    axes[1, 1].set_title('Sample Distribution', fontsize=11)
    axes[1, 1].grid(alpha=0.3)
    
    plt.suptitle(title, fontsize=14, y=1.00)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    gc.collect()
    print(f"Saved: {save_path}")


def plot_regression_results(y_true, y_pred, title="Regression Results", save_path=None):
    """Simple scatter plot for regression results"""
    plt.figure(figsize=(8, 6))
    plt.scatter(y_true, y_pred, alpha=0.5, color="blue", s=1)
    
    # Add perfect prediction line
    min_val, max_val = y_true.min(), y_true.max()
    plt.plot([min_val, max_val], [min_val, max_val], "r--", lw=2)
    
    plt.title(title)
    plt.xlabel("True Values")
    plt.ylabel("Predicted Values")
    plt.grid(alpha=0.3)
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    else:
        # Fallback to title-based filename if no save_path provided (matches train.py behavior)
        plt.savefig(f"{title}.png", dpi=300, bbox_inches='tight')
    plt.close()


def plot_reg_metrics(history, metrics=["rmse", "mse"], save_path="metrics_plot.png"):
    """Plot training and validation metrics over epochs"""
    plt.figure(figsize=(10, 6))
    for metric in metrics:
        if metric in history.history:
            plt.plot(history.history[metric], label=f"Train {metric.upper()}")
        val_metric = f"val_{metric}"
        if val_metric in history.history:
            plt.plot(history.history[val_metric], label=f"Validation {metric.upper()}")
    plt.title("Training and Validation Metrics")
    plt.xlabel("Epochs")
    plt.ylabel("Metric Value")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_prediction_intervals(y_true, y_pred, title="Prediction Intervals", save_path=None):
    """Plot showing individual prediction intervals (sampled for clarity)"""
    plt.figure(figsize=(14, 8))
    
    # Sort by true values
    sorted_idx = np.argsort(y_true)
    y_true_sorted = y_true[sorted_idx]
    y_pred_sorted = y_pred[sorted_idx]
    
    # Sample points for clarity (if too many)
    if len(y_true_sorted) > 500:
        step = len(y_true_sorted) // 500
        y_true_sorted = y_true_sorted[::step]
        y_pred_sorted = y_pred_sorted[::step]
    
    x_pos = np.arange(len(y_true_sorted))
    errors = np.abs(y_pred_sorted - y_true_sorted)
    
    plt.errorbar(x_pos, y_pred_sorted, yerr=errors, fmt='o', markersize=3, 
                 alpha=0.6, elinewidth=0.5, capsize=2, label='Predictions with Errors')
    plt.plot(x_pos, y_true_sorted, 'r-', lw=2, alpha=0.7, label='True Values')
    
    plt.xlabel('Sample Index (sorted by true value)', fontsize=12)
    plt.ylabel('Values', fontsize=12)
    plt.title(title, fontsize=14)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()


def plot_simple_scatter(y_true, y_pred, title="Regression Results", save_path=None, max_points=100000):
    """Simple but effective scatter plot - most memory efficient"""
    print(f"Creating simple scatter plot for {len(y_true)} points...")
    
    # Sample if necessary
    y_true_sampled, y_pred_sampled = smart_sample(y_true, y_pred, max_points)
    
    # Calculate metrics on full dataset
    r2, mae, rmse, pearson_r = calculate_metrics(y_true, y_pred)
    
    plt.figure(figsize=(10, 8))
    
    # Scatter plot with transparency
    plt.scatter(y_true_sampled, y_pred_sampled, alpha=0.3, s=3, color='blue')
    
    # Add perfect prediction line
    min_val = min(y_true_sampled.min(), y_pred_sampled.min())
    max_val = max(y_true_sampled.max(), y_pred_sampled.max())
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', lw=2)
    
    # Add text box with metrics
    textstr = f'R² = {r2:.3f}\nPearson r = {pearson_r:.3f}\nMAE = {mae:.3f}\nRMSE = {rmse:.3f}'
    plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, 
             fontsize=11, verticalalignment='top', bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
    
    plt.xlabel('True Values', fontsize=12)
    plt.ylabel('Predicted Values', fontsize=12)
    plt.title(title, fontsize=14)
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.close()
    gc.collect()
    print(f"Saved: {save_path}")


def plot_all_visualizations(y_train_true, y_train_pred, y_valid_true, y_valid_pred, prefix="model"):
    """Generate all visualizations - memory efficient version"""
        
    # Simple scatter plots (most reliable)
    plot_simple_scatter(y_train_true, y_train_pred, 
                       title=f"Training Data - Scatter Plot",
                       save_path=f"{prefix}_train_scatter.png")
    plot_simple_scatter(y_valid_true, y_valid_pred, 
                       title=f"Validation Data - Scatter Plot",
                       save_path=f"{prefix}_valid_scatter.png")
    
    # Hexbin plots
    plot_hexbin_regression(y_train_true, y_train_pred, 
                          title=f"Training Data - Hexbin Plot",
                          save_path=f"{prefix}_train_hexbin.png")
    plot_hexbin_regression(y_valid_true, y_valid_pred, 
                          title=f"Validation Data - Hexbin Plot",
                          save_path=f"{prefix}_valid_hexbin.png")
    
    # Residual analysis
    plot_residual_analysis(y_train_true, y_train_pred, 
                          title=f"Training Data - Residual Analysis",
                          save_path=f"{prefix}_train_residuals.png")
    plot_residual_analysis(y_valid_true, y_valid_pred, 
                          title=f"Validation Data - Residual Analysis",
                          save_path=f"{prefix}_valid_residuals.png")
    
    # Binned performance (uses full data)
    print("\n--- BINNED PERFORMANCE ---")
    plot_binned_performance(y_train_true, y_train_pred, 
                           title=f"Training Data - Binned Performance",
                           save_path=f"{prefix}_train_binned.png")
    plot_binned_performance(y_valid_true, y_valid_pred, 
                           title=f"Validation Data - Binned Performance",
                           save_path=f"{prefix}_valid_binned.png")
    
    # Error bands
    plot_error_bands(y_train_true, y_train_pred, 
                    title=f"Training Data - Error Bands",
                    save_path=f"{prefix}_train_error_bands.png")
    plot_error_bands(y_valid_true, y_valid_pred, 
                    title=f"Validation Data - Error Bands",
                    save_path=f"{prefix}_valid_error_bands.png")
        

if __name__ == "__main__":
    import argparse
    import tensorflow as tf
    
    parser = argparse.ArgumentParser(description="Generate regression visualizations from saved model")
    parser.add_argument("-model", type=str, required=True, help="Path to saved model (.h5)")
    parser.add_argument("-test_data", type=str, help="Path to test data pickle/directory")
    parser.add_argument("-labels", type=str, help="Path to labels CSV")
    parser.add_argument("-prefix", type=str, default="visualization", help="Prefix for output files")
    
    args = parser.parse_args()
    