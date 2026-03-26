import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import gaussian_kde
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
import gc
import matplotlib as mpl

def smart_sample(y_true, y_pred, max_points=200000):
    y_true = np.asarray(y_true)
    y_pred = np.asarray(y_pred)
    n = len(y_true)
    if n <= max_points:
        return y_true, y_pred
    n_bins = 100
    bins = np.linspace(y_true.min(), y_true.max(), n_bins + 1)
    bin_idx = np.digitize(y_true, bins) - 1
    bin_idx = np.clip(bin_idx, 0, n_bins - 1)
    points_per_bin = max(1, max_points // n_bins)
    keep = []
    for i in range(n_bins):
        inds = np.where(bin_idx == i)[0]
        if len(inds) == 0:
            continue
        k = min(points_per_bin, len(inds))
        keep.extend(np.random.choice(inds, k, replace=False))
    keep = np.unique(keep)
    return y_true[keep], y_pred[keep]


def calculate_metrics(y_true, y_pred):
    r2 = r2_score(y_true, y_pred)
    mae = mean_absolute_error(y_true, y_pred)
    rmse = np.sqrt(mean_squared_error(y_true, y_pred))
    pearson_r = np.corrcoef(y_true, y_pred)[0, 1]
    return r2, mae, rmse, pearson_r


def set_plot_style():
    sns.set_style("white")
    plt.rcParams.update({
        "figure.dpi": 150,
        "font.size": 11,
        "axes.titlesize": 14,
        "axes.labelsize": 12,
        "axes.labelweight": "bold",
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "lines.linewidth": 1.6,
        "axes.linewidth": 0.9,
        "axes.edgecolor": "#333333",
        "xtick.direction": "out",
        "ytick.direction": "out",
    })


def _infer_true_pred_cols(df):
    lower = {c.lower(): c for c in df.columns}
    true_cands = ['y_true', 'true', 'target', 'label', 'ground_truth', 'gt']
    pred_cands = ['y_pred', 'pred', 'prediction', 'yhat', 'y_hat', 'y_predicted']
    tcol = next((lower[c] for c in true_cands if c in lower), None)
    pcol = next((lower[c] for c in pred_cands if c in lower), None)
    if tcol is None or pcol is None:
        numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
        if len(numeric_cols) >= 2:
            tcol, pcol = numeric_cols[0], numeric_cols[1]
        else:
            raise ValueError("Couldn't infer true/pred columns.")
    return tcol, pcol



def plot_density_grid_from_csvs(csv_paths, labels=None, ncols=2, max_points=200000,
                                figsize=(12, 10), save_path=None, cmap='viridis_r'):
    """
    Grid of density-colored scatter plots (one subplot per CSV).
    Uses a single shared colorbar on the far right to avoid overlapping with axes.
    Colorbar is on a log scale (LogNorm). Safeguards added for zero/non-positive KDE outputs.
    """
    # keep your style (bold labels etc.)
    set_plot_style()

    if labels is None:
        labels = [p.split('/')[-1] for p in csv_paths]

    n = len(csv_paths)
    ncols = max(1, ncols)
    nrows = int(np.ceil(n / ncols))

    # load data + compute KDE densities for sampled points, store for plotting
    datasets = []
    global_min = np.inf
    global_max = -np.inf
    global_zmin = np.inf
    global_zmax = -np.inf
    global_zposmin = np.inf  # track smallest positive density for LogNorm

    for p in csv_paths:
        df = pd.read_csv(p)
        tcol, pcol = _infer_true_pred_cols(df)
        x_all = df[tcol].to_numpy(dtype=float)
        y_all = df[pcol].to_numpy(dtype=float)
        # sample for plotting (use your smart_sample if present)
        try:
            xs, ys = smart_sample(x_all, y_all, max_points=max_points)
        except Exception:
            # fallback uniform sample
            npts = len(x_all)
            if npts > max_points:
                idx = np.random.choice(npts, max_points, replace=False)
                xs, ys = x_all[idx], y_all[idx]
            else:
                xs, ys = x_all, y_all

        # compute KDE density for sampled points (fallback to ones)
        xy = np.vstack([xs, ys])
        try:
            z = gaussian_kde(xy)(xy)
        except Exception:
            z = np.ones(xs.shape[0])

        # update global ranges for diagonal scaling and density normalization
        global_min = min(global_min, np.nanmin(xs), np.nanmin(ys))
        global_max = max(global_max, np.nanmax(xs), np.nanmax(ys))
        if np.isfinite(z).any():
            # for normalization display (may include zeros)
            finite_z = z[np.isfinite(z)]
            if finite_z.size:
                global_zmin = min(global_zmin, np.nanmin(finite_z))
                global_zmax = max(global_zmax, np.nanmax(finite_z))
            # consider positive values only for LogNorm vmin
            pos = z[(np.isfinite(z)) & (z > 0)]
            if pos.size:
                global_zposmin = min(global_zposmin, pos.min())

        datasets.append({
            "x_all": x_all,
            "y_all": y_all,
            "xs": xs,
            "ys": ys,
            "z": z
        })

    # safety defaults if KDE failed everywhere
    if not np.isfinite(global_zmin) or not np.isfinite(global_zmax) or global_zmin == global_zmax:
        global_zmin, global_zmax = 0.0, 1.0

    # Determine vmin for LogNorm: must be > 0
    if np.isfinite(global_zposmin) and global_zposmin > 0:
        vmin_log = float(global_zposmin)
    else:
        # If no positive densities found, pick a small positive fallback.
        # Use a fraction of global_zmax if available, otherwise 1e-3.
        if np.isfinite(global_zmax) and global_zmax > 0:
            vmin_log = float(max(global_zmax * 1e-6, np.nextafter(0.0, 1.0)))
        else:
            vmin_log = 1e-3

    # If somehow vmin is equal to vmax, widen the interval
    vmax_log = float(global_zmax) if np.isfinite(global_zmax) and global_zmax > 0 else vmin_log * 100.0
    if vmin_log >= vmax_log:
        vmax_log = vmin_log * 100.0

    pad = 0.05 * (global_max - global_min) if global_max > global_min else 1.0
    xmin, xmax = global_min - pad, global_max + pad

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
    axes = axes.flatten()

    # Use a logarithmic normalization for color mapping
    norm = mpl.colors.LogNorm(vmin=vmin_log, vmax=vmax_log)
    cmap_obj = plt.get_cmap(cmap)

    for i, data in enumerate(datasets):
        ax = axes[i]
        xs, ys, z = data["xs"], data["ys"], data["z"]

        # To avoid passing non-positive values to LogNorm, clip z for plotting
        z_plot = np.array(z, copy=True)
        # replace non-finite or non-positive with small positive number for plotting (they will all map near vmin)
        invalid_mask = ~np.isfinite(z_plot) | (z_plot <= 0)
        if invalid_mask.any():
            z_plot[invalid_mask] = vmin_log

        # sort by density so dense points draw on top
        order = np.argsort(z_plot)
        xs_s, ys_s, z_s = xs[order], ys[order], z_plot[order]

        sc = ax.scatter(xs_s, ys_s, c=z_s, s=18, cmap=cmap_obj, norm=norm,
                        alpha=0.85, edgecolors='none')

        # diagonal perfect line
        ax.plot([xmin, xmax], [xmin, xmax], ls='--', lw=2.2, color='firebrick', alpha=0.9)

        # subplot labels / title
        ax.set_title(labels[i], fontsize=12, fontweight='bold')
        ax.set_xlabel("True Values", fontsize=11, fontweight='bold')
        ax.set_ylabel("Predicted Values", fontsize=11, fontweight='bold')

        ax.set_xlim(xmin, xmax)
        ax.set_ylim(xmin, xmax)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.grid(False)

    # hide any empty subplots
    for j in range(len(datasets), len(axes)):
        axes[j].set_visible(False)

    # place a single shared colorbar on the far right (avoids overlap with y-axes)
    mappable = mpl.cm.ScalarMappable(norm=norm, cmap=cmap_obj)
    mappable.set_array([])  # required for colorbar

    plt.tight_layout()
    plt.subplots_adjust(right=0.88, hspace=0.28, wspace=0.28)  # leave space for colorbar
    cbar = fig.colorbar(mappable, ax=axes.tolist(), fraction=0.025, pad=0.02)
    cbar.set_label('Density (log scale)', fontsize=10, fontweight='bold')
    # Set log ticks and formatter
    try:
        cbar.ax.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0))
        cbar.ax.yaxis.set_major_formatter(mpl.ticker.LogFormatterSciNotation())
        cbar.update_ticks()
    except Exception:
        # fallback - if locator fails (e.g. tiny dynamic range), let matplotlib handle it
        pass
    cbar.ax.tick_params(labelsize=9)

    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")

    plt.show()
    plt.close(fig)
    gc.collect()


def plot_density_regression_from_csv(csv_path, title=None, max_points=200000, save_path=None, cmap='viridis_r'):
    set_plot_style()

    df = pd.read_csv(csv_path)
    tcol, pcol = _infer_true_pred_cols(df)
    y_true = df[tcol].to_numpy(dtype=float)
    y_pred = df[pcol].to_numpy(dtype=float)

    y_true_s, y_pred_s = smart_sample(y_true, y_pred, max_points=max_points)
    xy = np.vstack([y_true_s, y_pred_s])

    try:
        z = gaussian_kde(xy)(xy)
    except Exception:
        z = np.ones_like(y_true_s)

    # prepare z for log-display: find smallest positive
    pos = z[(np.isfinite(z)) & (z > 0)]
    if pos.size:
        vmin_log = float(pos.min())
    else:
        if np.isfinite(z).any():
            vmax_guess = float(np.nanmax(z[np.isfinite(z)]))
            vmin_log = float(max(vmax_guess * 1e-6, np.nextafter(0.0, 1.0))) if vmax_guess > 0 else 1e-3
        else:
            vmin_log = 1e-3
    vmax_log = float(np.nanmax(z[np.isfinite(z)]) if np.isfinite(z).any() and np.nanmax(z[np.isfinite(z)]) > 0 else (vmin_log * 100.0))
    if vmin_log >= vmax_log:
        vmax_log = vmin_log * 100.0

    idx = z.argsort()
    y_true_s, y_pred_s, z = y_true_s[idx], y_pred_s[idx], z[idx]

    # global bounds for diagonal
    pad = 0.05 * (max(y_true_s.max(), y_pred_s.max()) - min(y_true_s.min(), y_pred_s.min()))
    xmin = min(y_true_s.min(), y_pred_s.min()) - pad
    xmax = max(y_true_s.max(), y_pred_s.max()) + pad

    fig, ax = plt.subplots(figsize=(6, 5))

    # ensure all z values we pass to the scatter are finite & > 0
    z_plot = np.array(z, copy=True)
    invalid_mask = ~np.isfinite(z_plot) | (z_plot <= 0)
    if invalid_mask.any():
        z_plot[invalid_mask] = vmin_log

    norm = mpl.colors.LogNorm(vmin=vmin_log, vmax=vmax_log)
    sc = ax.scatter(y_true_s, y_pred_s, c=z_plot, s=20, cmap=cmap, norm=norm, alpha=0.85, edgecolors='none')
    ax.plot([xmin, xmax], [xmin, xmax], ls='--', lw=2.2, color='firebrick', alpha=0.9)

    ax.set_xlabel("True Values", fontsize=11, fontweight='bold')
    ax.set_ylabel("Predicted Values", fontsize=11, fontweight='bold')
    if title:
        ax.set_title(title, fontsize=13, fontweight='bold')

    ax.set_xlim(xmin, xmax)
    ax.set_ylim(xmin, xmax)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(False)

    # add colorbar on the right (log scale)
    cbar = fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label('Density (log scale)', fontsize=10, fontweight='bold')
    try:
        cbar.ax.yaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0))
        cbar.ax.yaxis.set_major_formatter(mpl.ticker.LogFormatterSciNotation())
        cbar.update_ticks()
    except Exception:
        pass
    cbar.ax.tick_params(labelsize=9)

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Saved: {save_path}")

    plt.show()
    plt.close(fig)
    gc.collect()


# Example usage:
# plot_density_regression_from_csv("test_predictions.csv", title="Test Predictions", save_path="test_density.png")

if __name__ == "__main__":
    # Example usage
    csv_files = [
        "test_predictions1.csv",
        "test_predictions2.csv",
        "test_predictions4.csv",
        "test_predictions5.csv",
    ]
    labels = [
        "Fold 1",
        "Fold 2",
        "Fold 4",
        "Fold 5"
    ]
    plot_density_grid_from_csvs(csv_files, labels=labels, save_path="overlay_density_plot.png")    
    plot_density_regression_from_csv("test_predictions.csv", title="Test Predictions", save_path="test_density.png")