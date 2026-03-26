import os
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"

import tensorflow as tf
import pandas as pd
import numpy as np
import joblib
import glob
import dask.dataframe as dd
import dask
from argparse import ArgumentParser
import sys

def PCC_RMSE(y_true, y_pred):
    """Custom loss function - must match training definition"""
    alpha = 0.7  # Default value, adjust if you used different alpha
    eps = 1e-7

    y_true = tf.cast(y_true, tf.float32)
    y_pred = tf.cast(y_pred, tf.float32)
    
    y_true = tf.where(tf.math.is_finite(y_true), y_true, tf.zeros_like(y_true))
    y_pred = tf.where(tf.math.is_finite(y_pred), y_pred, tf.zeros_like(y_pred))
    
    y_true_mean = tf.reduce_mean(y_true)
    y_pred_mean = tf.reduce_mean(y_pred)
    
    fsp = y_pred - y_pred_mean
    fst = y_true - y_true_mean
    
    cov = tf.reduce_mean(fsp * fst)
    devP = tf.math.reduce_std(y_pred) + eps
    devT = tf.math.reduce_std(y_true) + eps
    denom = devP * devT
    
    corr = cov / denom
    corr = tf.clip_by_value(corr, -1.0, 1.0)
    corr = tf.where(tf.math.is_finite(corr), corr, 0.0)
    
    pcc_loss_term = 1.0 - corr
    squared_diff = tf.square(y_pred - y_true)
    squared_diff = tf.clip_by_value(squared_diff, 0.0, 1e6)
    rmse = tf.sqrt(tf.reduce_mean(squared_diff) + eps)
    
    loss = alpha * pcc_loss_term + (1.0 - alpha) * rmse
    loss = tf.where(tf.math.is_finite(loss), loss, tf.constant(1.0, dtype=tf.float32))
    return loss


def pearson_correlation(y_true, y_pred):
    """Pearson correlation metric"""
    eps = 1e-7
    
    y_true = tf.where(tf.math.is_finite(y_true), y_true, tf.zeros_like(y_true))
    y_pred = tf.where(tf.math.is_finite(y_pred), y_pred, tf.zeros_like(y_pred))
    
    y_true_mean = tf.reduce_mean(y_true)
    y_pred_mean = tf.reduce_mean(y_pred)
    
    y_true_centered = y_true - y_true_mean
    y_pred_centered = y_pred - y_pred_mean
    
    covariance = tf.reduce_sum(y_true_centered * y_pred_centered)
    y_true_std = tf.sqrt(tf.reduce_sum(tf.square(y_true_centered)) + eps)
    y_pred_std = tf.sqrt(tf.reduce_sum(tf.square(y_pred_centered)) + eps)
    
    corr = covariance / (y_true_std * y_pred_std)
    corr = tf.clip_by_value(corr, -1.0, 1.0)
    corr = tf.where(tf.math.is_finite(corr), corr, 0.0)
    
    return corr


def build_dask_dataframe(pickle_path):
    """Build Dask DataFrame from pickle file(s)"""
    if os.path.isdir(pickle_path):
        pickle_files = sorted(glob.glob(os.path.join(pickle_path, "*.pkl")))
        if len(pickle_files) == 0:
            raise ValueError(f"No pickle files found in directory: {pickle_path}")
        print(f"Found {len(pickle_files)} pickle files in {pickle_path}")
        delayed_dfs = [dask.delayed(pd.read_pickle)(pf) for pf in pickle_files]
        ddf = dd.from_delayed(delayed_dfs)
    else:
        print(f"Loading single pickle file: {pickle_path}")
        ddf = dd.from_delayed([dask.delayed(pd.read_pickle)(pickle_path)])
    
    return ddf


def data_generator(ddf, n_features, scaler, shape, label_series=None):
    """Generator that yields (features, labels) or just features if no labels"""
    for partition in ddf.to_delayed():
        df_part = partition.compute()
        if df_part.empty:
            continue
        
        keys = df_part["pdb"]
        X_part = df_part.iloc[:, 1: n_features+1].fillna(0).values.astype(np.float32)
        
        if np.isnan(X_part).any() or np.isinf(X_part).any():
            X_part = np.nan_to_num(X_part, nan=0.0, posinf=0.0, neginf=0.0)
        
        try:
            X_std = scaler.transform(X_part)
            if np.isnan(X_std).any() or np.isinf(X_std).any():
                X_std = np.nan_to_num(X_std, nan=0.0, posinf=0.0, neginf=0.0)
        except Exception as e:
            print(f"Error during scaling: {e}")
            continue
        
        for key, row in zip(keys, X_std):
            try:
                # Clip extreme values
                row = np.clip(row, -10, 10)
                reshaped_row = row.reshape(shape)
                
                if np.isnan(reshaped_row).any() or np.isinf(reshaped_row).any():
                    continue
                
                if label_series is not None:
                    # If labels are provided, yield (features, label)
                    try:
                        label = label_series.loc[key]
                        if np.isnan(label) or np.isinf(label):
                            continue
                        yield reshaped_row, float(label)
                    except KeyError:
                        print(f"Warning: Key {key} not found in labels")
                        continue
                else:
                    # If no labels, yield just features and key
                    yield key, reshaped_row
                    
            except Exception as e:
                print(f"Error processing key {key}: {e}")
                continue


def load_model_and_predict(model_path, data_path, scaler_path, label_path=None, 
                          n_features=10416, shape=(84, 124, 1), batch_size=64, 
                          output_csv="predictions.csv"):
    """
    Load saved model and generate predictions
    
    Args:
        model_path: Path to saved .h5 model
        data_path: Path to pickle file or directory containing pickle files
        scaler_path: Path to saved scaler (.scaler file)
        label_path: Optional path to CSV with true labels (for evaluation)
        n_features: Number of features (default 10416)
        shape: Shape to reshape features (default (84, 124, 1))
        batch_size: Batch size for prediction (default 64)
        output_csv: Output CSV file for predictions
    """
        
    # Load the model with custom objects
    print(f"\nLoading model from: {model_path}")
    custom_objects = {
        'PCC_RMSE': PCC_RMSE,
        'pearson_correlation': pearson_correlation
    }
    
    try:
        # Try loading directly
        model = tf.keras.models.load_model(model_path, custom_objects=custom_objects)
        print("Model loaded successfully!")
    except (IOError, OSError) as e:
        print(f"Initial load failed: {e}")
        print("Attempting alternative loading method...")
        
        # Try loading without compilation, then recompile
        try:
            model = tf.keras.models.load_model(model_path, custom_objects=custom_objects, compile=False)
            print("Model loaded without compilation. Recompiling...")
            
            model.compile(
                optimizer=tf.keras.optimizers.SGD(learning_rate=0.001, momentum=0.9, clipnorm=1.0),
                loss=PCC_RMSE,
                metrics=[
                    tf.keras.metrics.MeanSquaredError(name="mse"),
                    tf.keras.metrics.MeanAbsoluteError(name="mae"),
                    tf.keras.metrics.RootMeanSquaredError(name="rmse"),
                    pearson_correlation
                ]
            )
            print("Model recompiled successfully!")
        except Exception as e2:
            print(f"ERROR: Could not load model: {e2}")
            print("\nTroubleshooting tips:")
            print("1. Check if model path is correct")
            print("2. If it's a directory, it should contain 'saved_model.pb'")
            print("3. If it's an HDF5 file, it should be a single .h5 file")
            raise
    
    
    # Load the scaler
    scaler = joblib.load(scaler_path)
    
    # Build Dask DataFrame
    print(f"\nLoading data from: {data_path}")
    ddf = build_dask_dataframe(data_path)
    
    # Load labels if provided
    label_series = None
    if label_path:
        label_df = pd.read_csv(label_path, index_col="name")
        label_series = label_df["binding_aff"]
    
    # Create dataset    
    if label_series is not None:
        # With labels - for evaluation
        def gen():
            return data_generator(ddf, n_features, scaler, shape, label_series)
        
        output_signature = (
            tf.TensorSpec(shape=shape, dtype=tf.float32),
            tf.TensorSpec(shape=(), dtype=tf.float32)
        )
        
        dataset = tf.data.Dataset.from_generator(gen, output_signature=output_signature)
        dataset = dataset.batch(batch_size).prefetch(tf.data.AUTOTUNE)
        
        # Make predictions and collect true labels
        predictions = []
        true_labels = []
        
        for i, (batch_x, batch_y) in enumerate(dataset):
            pred = model.predict(batch_x, verbose=0)
            predictions.append(pred)
            true_labels.append(batch_y.numpy())
            
            if (i + 1) % 10 == 0:
        
        predictions = np.concatenate(predictions).flatten()
        true_labels = np.concatenate(true_labels)
        
        # Calculate metrics
        from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
        
        r2 = r2_score(true_labels, predictions)
        mae = mean_absolute_error(true_labels, predictions)
        rmse = np.sqrt(mean_squared_error(true_labels, predictions))
        pearson_r = np.corrcoef(true_labels, predictions)[0, 1]
        
        print("\n" + "="*60)
        print("EVALUATION METRICS")
        print("="*60)
        print(f"R² Score:           {r2:.4f}")
        print(f"Pearson Correlation: {pearson_r:.4f}")
        print(f"MAE:                {mae:.4f}")
        print(f"RMSE:               {rmse:.4f}")
        print(f"Total samples:      {len(predictions)}")
        print("="*60)
        
        # Save predictions
        results_df = pd.DataFrame({
            'true_value': true_labels,
            'predicted_value': predictions,
            'absolute_error': np.abs(predictions - true_labels)
        })
        
    else:
        # Without labels - just predictions
        def gen():
            return data_generator(ddf, n_features, scaler, shape, None)
        
        output_signature = (
            tf.TensorSpec(shape=(), dtype=tf.string),  # key
            tf.TensorSpec(shape=shape, dtype=tf.float32)  # features
        )
        
        keys_list = []
        predictions = []
        
        # Process without batching to preserve keys
        for key, features in gen():
            pred = model.predict(np.expand_dims(features, axis=0), verbose=0)
            keys_list.append(key)
            predictions.append(pred[0, 0])
            
            if len(predictions) % 100 == 0:
        
        predictions = np.array(predictions)
        
        
        # Save predictions
        results_df = pd.DataFrame({
            'pdb': keys_list,
            'predicted_value': predictions
        })
    
    # Save to CSV
    results_df.to_csv(output_csv, index=False)
    print(f"\nResults saved to: {output_csv}")
    
    return results_df, predictions, (true_labels if label_series is not None else None)


def main():
    parser = ArgumentParser(description="Load trained model and make predictions on new data")
    
    parser.add_argument("-model", type=str, required=True,
                       help="Path to saved model (.h5 file)")
    parser.add_argument("-data", type=str, required=True,
                       help="Path to pickle file or directory with pickle files")
    parser.add_argument("-scaler", type=str, required=True,
                       help="Path to saved scaler (.scaler file)")
    parser.add_argument("-labels", type=str, default=None,
                       help="Optional: CSV file with true labels (for evaluation)")
    parser.add_argument("-n_features", type=int, default=10416,
                       help="Number of features (default: 10416)")
    parser.add_argument("-shape", type=str, default="84,124,1",
                       help="Shape to reshape features (default: '84,124,1')")
    parser.add_argument("-batch_size", type=int, default=64,
                       help="Batch size for prediction (default: 64)")
    parser.add_argument("-output", type=str, default="predictions.csv",
                       help="Output CSV file (default: predictions.csv)")
    parser.add_argument("-visualize", action="store_true",
                       help="Generate visualization plots (requires labels)")
    parser.add_argument("-plot_prefix", type=str, default="inference",
                       help="Prefix for plot files (default: inference)")
    
    args = parser.parse_args()
    
    # Parse shape
    shape = tuple(int(x) for x in args.shape.split(","))
    
    # Make predictions
    results_df, predictions, true_labels = load_model_and_predict(
        model_path=args.model,
        data_path=args.data,
        scaler_path=args.scaler,
        label_path=args.labels,
        n_features=args.n_features,
        shape=shape,
        batch_size=args.batch_size,
        output_csv=args.output
    )
    
    # Generate visualizations if requested and labels are available
    if args.visualize and true_labels is not None:
        try:
            from regression_viz import plot_all_visualizations
            
            # For inference, we only have one dataset
            plot_all_visualizations(
                true_labels, predictions,
                true_labels, predictions,  # Use same data for both
                prefix=args.plot_prefix
            )

if __name__ == "__main__":
    main()