import os
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"  

import tensorflow as tf
import pandas as pd
import numpy as np
from argparse import RawDescriptionHelpFormatter
import argparse
import joblib
from sklearn import preprocessing
import matplotlib.pyplot as plt
import dask.dataframe as dd
import dask
import glob
import random
import seaborn as sns
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
from scipy import stats
from regression_viz import plot_all_visualizations, plot_reg_metrics, plot_regression_results


SEED = 42
os.environ['PYTHONHASHSEED'] = str(SEED)
random.seed(SEED)
np.random.seed(SEED)
tf.random.set_seed(SEED)

def PCC_RMSE(y_true, y_pred):
    alpha = args.alpha
    eps = 1e-7

    # ensure float32
    y_true = tf.cast(y_true, tf.float32)
    y_pred = tf.cast(y_pred, tf.float32)


    # center
    y_true_mean = tf.reduce_mean(y_true)
    y_pred_mean = tf.reduce_mean(y_pred)
    
    fsp = y_pred - y_pred_mean
    fst = y_true - y_true_mean

    # covariance
    cov = tf.reduce_mean(fsp * fst)

    # standard deviations with epsilon guard
    devP = tf.math.reduce_std(y_pred) + eps
    devT = tf.math.reduce_std(y_true) + eps
    denom = devP * devT

    # pearson correlation (clamp to [-1, 1])
    corr = cov / denom
    corr = tf.clip_by_value(corr, -1.0, 1.0)
    

    pcc_loss_term = 1.0 - corr

    # RMSE with clipping to prevent extreme values
    squared_diff = tf.square(y_pred - y_true)
    rmse = tf.sqrt(tf.reduce_mean(squared_diff) + eps)

    loss = alpha * pcc_loss_term + (1.0 - alpha) * rmse
    
    return loss


def RMSE(y_true, y_pred):
    return tf.keras.backend.sqrt(tf.keras.backend.mean(tf.keras.backend.square(y_pred - y_true), axis=-1))

def PCC(y_true, y_pred):
    eps = 1e-7
    fsp = y_pred - tf.keras.backend.mean(y_pred)
    fst = y_true - tf.keras.backend.mean(y_true)
    devP = tf.keras.backend.std(y_pred) + eps
    devT = tf.keras.backend.std(y_true) + eps
    pcc = tf.keras.backend.mean(fsp * fst) / (devP * devT)
    return pcc

def pearson_correlation(y_true, y_pred):
    eps = 1e-7
    
    
    # Compute mean
    y_true_mean = tf.reduce_mean(y_true)
    y_pred_mean = tf.reduce_mean(y_pred)
    
    # Center
    y_true_centered = y_true - y_true_mean
    y_pred_centered = y_pred - y_pred_mean
    
    # Covariance and standard deviations
    covariance = tf.reduce_sum(y_true_centered * y_pred_centered)
    y_true_std = tf.sqrt(tf.reduce_sum(tf.square(y_true_centered)) + eps)
    y_pred_std = tf.sqrt(tf.reduce_sum(tf.square(y_pred_centered)) + eps)
    
    corr = covariance / (y_true_std * y_pred_std)
    corr = tf.clip_by_value(corr, -1.0, 1.0)
    
    
    return corr


def create_model(input_size, lr=0.001):
    model = tf.keras.Sequential()
    model.add(tf.keras.layers.Conv2D(32, 4, 1, input_shape=input_size))
    model.add(tf.keras.layers.Activation("relu"))
    model.add(tf.keras.layers.BatchNormalization())  # Moved before dropout
    model.add(tf.keras.layers.Dropout(rate=args.rate))
    
    model.add(tf.keras.layers.Conv2D(64, 4, 1))
    model.add(tf.keras.layers.Activation("relu"))
    model.add(tf.keras.layers.BatchNormalization())  # Moved before dropout
    model.add(tf.keras.layers.Dropout(rate=args.rate))
    
    model.add(tf.keras.layers.Conv2D(128, 4, 1))
    model.add(tf.keras.layers.Activation("relu"))
    model.add(tf.keras.layers.BatchNormalization())  # Moved before dropout
    model.add(tf.keras.layers.Dropout(rate=args.rate))
    
    model.add(tf.keras.layers.Flatten())
    
    model.add(tf.keras.layers.Dense(100, kernel_regularizer=tf.keras.regularizers.l2(0.01)))
    model.add(tf.keras.layers.Activation("relu"))
    model.add(tf.keras.layers.BatchNormalization())  # Moved before dropout
    model.add(tf.keras.layers.Dropout(rate=args.rate))
    
    model.add(tf.keras.layers.Dense(50, kernel_regularizer=tf.keras.regularizers.l2(0.01)))
    model.add(tf.keras.layers.Activation("relu"))
    model.add(tf.keras.layers.BatchNormalization())  # Moved before dropout
    model.add(tf.keras.layers.Dropout(rate=args.rate))
    
    model.add(tf.keras.layers.Dense(1, kernel_regularizer=tf.keras.regularizers.l2(0.01)))
    model.add(tf.keras.layers.Activation("linear"))
    
    # Use gradient clipping in optimizer
    sgd = tf.keras.optimizers.SGD(
        learning_rate=lr, 
        momentum=0.9, 
        clipnorm=1.0  # Use clipnorm for per-gradient clipping
    )
    
    model.compile(
        optimizer=sgd, 
        loss=PCC_RMSE, 
        metrics=[
            tf.keras.metrics.MeanSquaredError(name="mse"),
            tf.keras.metrics.MeanAbsoluteError(name="mae"),
            tf.keras.metrics.RootMeanSquaredError(name="rmse"),
            pearson_correlation
        ]
    )
    return model

def build_dask_dataframe(path_or_pattern):
    """
    Given a directory path or glob pattern containing pickle files,
    build and return a Dask DataFrame.
    
    Examples:
        - "/scratch/data/set1" -> loads all .pkl files in set1/
        - "/scratch/data/set*[123]*" -> loads all .pkl files in set1/, set2/, set3/
        - "/scratch/data/*.pkl" -> loads all .pkl files directly
    """
    pickle_files = []
    
    # Check if it's a plain directory
    if os.path.isdir(path_or_pattern):
        pickle_files = sorted(glob.glob(os.path.join(path_or_pattern, "*.pkl")))
    else:
        # It's a glob pattern - expand it
        expanded_paths = glob.glob(path_or_pattern)
        
        for path in expanded_paths:
            if os.path.isdir(path):
                # If pattern expanded to directories, get .pkl files from each
                pickle_files.extend(sorted(glob.glob(os.path.join(path, "*.pkl"))))
            elif path.endswith('.pkl'):
                # If pattern matched .pkl files directly
                pickle_files.append(path)
        
        pickle_files = sorted(pickle_files)
    
    
    print(f"Found {len(pickle_files)} pickle files to load")
    delayed_dfs = [dask.delayed(pd.read_pickle)(pf) for pf in pickle_files]
    ddf = dd.from_delayed(delayed_dfs)
    return ddf

def compute_scaler(ddf, n_features, batch_size):
    scaler = preprocessing.StandardScaler()
    for part in ddf.to_delayed():
        df = part.compute()
        # skip if fewer than 2 rows
        if df.shape[0] < 2:
            continue

        X = df.iloc[:, 1 : n_features+1].fillna(0).values.astype(np.float32)
        if X.size == 0:
            continue

        scaler.partial_fit(X)

    # Ensure all variances are positive
    scaler.var_ = np.maximum(scaler.var_, 1e-10)
    scaler.scale_ = np.sqrt(scaler.var_)
    
    return scaler


def data_generator(ddf, n_features, scaler, shape, label_series):
    """
    Generator that yields standardized (features, label) pairs with robust NaN handling.
    """
    for partition in ddf.to_delayed():
        df_part = partition.compute()
        if df_part.empty:
            continue
        # Get keys from the "pdb" column
        keys = df_part["pdb"]
        # Extract features (columns 1 to n_features+1) and fill NaNs
        X_part = df_part.iloc[:, 1: n_features+1].fillna(0).values.astype(np.float32)
        
        X_std = scaler.transform(X_part)
            
        for key, row in zip(keys, X_std):
                label = label_series.loc[key]
                reshaped_row = row.reshape(shape)
                yield reshaped_row, float(label)

if __name__ == "__main__":
    print("Start training ...")
    
    d = """
    This architecture is built with TensorFlow 2.x.
    
    Train the CNN model using features extracted from complexes.
    Training and validation features are loaded via Dask from directories containing pickle files.
    
    Usage:
      python train.py -train_fpath train_features_dir -valid_fpath valid_features_dir \
         -label_fpath_train train_labels.csv -label_fpath_validate valid_labels.csv
    """
    
    parser = argparse.ArgumentParser(description=d, formatter_class=RawDescriptionHelpFormatter)
    
    parser.add_argument("-train_fpath", type=str, default="train_features.pkl",
                        help="Input file or directory containing training features (pickle files).")
    parser.add_argument("-valid_fpath", type=str, default="valid_features.pkl",
                        help="Input file or directory containing validation features (pickle files).")
    parser.add_argument("-label_fpath_train", type=str, default="PDBbindv2019_pKd-label.csv",
                        help="CSV file containing training labels. Must have a 'name' column as index.")
    parser.add_argument("-label_fpath_validate", type=str, default="PDBbindv2019_pKd-label.csv",
                        help="CSV file containing validation labels. Must have a 'name' column as index.")
    parser.add_argument("-shape", type=str, default="84,124,1",
                        help="Reshape the features; e.g., '84,124,1'.")
    parser.add_argument("-lr", type=float, default=0.001,
                        help="Learning rate.")
    parser.add_argument("-batch_size", type=int, default=64,
                        help="Batch size for training.")
    parser.add_argument("-rate", type=float, default=0.0,
                        help="Dropout rate.")
    parser.add_argument("-alpha", type=float, default=0.7,
                        help="Alpha value in loss function.")
    parser.add_argument("-clipvalue", type=float, default=0.01,
                        help="Threshold for gradient clipping.")
    parser.add_argument("-n_features", type=int, default=10416,
                        help="Number of features for each complex. For shells N=62, n_features=21*8*62.")
    parser.add_argument("-epochs", type=int, default=300,
                        help="Number of training epochs.")
    parser.add_argument("-patience", type=int, default=30,
                        help="Epochs with no improvement after which training stops.")
    parser.add_argument("-out_model", type=str, default="bestmodel.h5",
                        help="File path to save the best model.")
    parser.add_argument("-plotas", type=str, default="example", help="Prefix for plot file names.")
    
    args = parser.parse_args()
    shape = [int(x) for x in args.shape.split(",")]
    n_features = args.n_features
    batch_size = args.batch_size
    
    ddf_train = build_dask_dataframe(args.train_fpath)
    ddf_valid = build_dask_dataframe(args.valid_fpath)

    # Load label CSV files; assume the index column is named "name"
    label_df_train = pd.read_csv(args.label_fpath_train, index_col="name")
    label_df_validate = pd.read_csv(args.label_fpath_validate, index_col="name")
    
    # Compute the StandardScaler over the training set (batchwise)
    scaler = compute_scaler(ddf_train, n_features, batch_size)

    # replace any NaN variance with 0 (so scale = sqrt(0) = 0)
    scaler.var_  = np.nan_to_num(scaler.var_,  nan=0.0, posinf=0.0, neginf=0.0)
    scaler.scale_ = np.sqrt(scaler.var_)
    # guard against scale==0
    scaler.scale_[scaler.scale_ == 0] = 1.0

    joblib.dump(scaler, 'train_scaler_3.scaler')
    
    # Create tf.data.Datasets using a generator.
    def train_generator():
        return data_generator(ddf_train, n_features, scaler, shape, label_df_train["binding_aff"])
    
    def valid_generator():
        return data_generator(ddf_valid, n_features, scaler, shape, label_df_validate["binding_aff"])
    
    output_signature = (tf.TensorSpec(shape=shape, dtype=tf.float32),
                        tf.TensorSpec(shape=(), dtype=tf.float32))
    
    train_dataset = tf.data.Dataset.from_generator(train_generator, output_signature=output_signature)
    valid_dataset = tf.data.Dataset.from_generator(valid_generator, output_signature=output_signature)
    

    
    train_dataset = train_dataset.batch(batch_size).prefetch(tf.data.AUTOTUNE)
    valid_dataset = valid_dataset.batch(batch_size).prefetch(tf.data.AUTOTUNE)

    # Create model
    model = create_model(shape, args.lr)
    
    
    # Add TerminateOnNaN for stability
    class TerminateOnNaN(tf.keras.callbacks.Callback):
        def on_batch_end(self, batch, logs=None):
            if np.isnan(logs.get('loss')) or np.isinf(logs.get('loss')):
                print(f'\nBatch {batch}: Invalid loss, terminating training')
                self.model.stop_training = True
    
    callbacks = [
        TerminateOnNaN(),
        tf.keras.callbacks.EarlyStopping(
            monitor='val_loss',
            patience=args.patience,
            restore_best_weights=True,
            verbose=1
        ),
        tf.keras.callbacks.ModelCheckpoint(
            args.out_model,
            monitor='val_loss',
            save_best_only=True,
            verbose=1
        ),
        tf.keras.callbacks.ReduceLROnPlateau(
            monitor='val_loss',
            factor=0.5,
            patience=10,
            min_lr=1e-7,
            verbose=1
        ),
    ]
    
    history = model.fit(
        train_dataset,
        validation_data=valid_dataset,
        epochs=args.epochs,
        verbose=1,
        callbacks=callbacks
    )
    
    # Gather predictions for plotting.
    print("Generating predictions...")
    y_train_pred = np.concatenate([model.predict(batch[0]) for batch in train_dataset])
    y_valid_pred = np.concatenate([model.predict(batch[0]) for batch in valid_dataset])
    
    y_train_true = np.concatenate([batch[1].numpy() for batch in train_dataset])
    y_valid_true = np.concatenate([batch[1].numpy() for batch in valid_dataset])
    
    print("Plotting results...")
    plot_regression_results(y_train_true, y_train_pred, title=f"Regression of Training Data!wholesim PCCl2+relu2 - {args.plotas}")
    plot_regression_results(y_valid_true, y_valid_pred, title=f"Regression of Validation Data!wholesim PCCl2+relu2- {args.plotas}")
    plot_reg_metrics(history, save_path=f"metrics_over_epochs_{args.plotas}.png")
    plot_all_visualizations(y_train_true, y_train_pred, y_valid_true, y_valid_pred, prefix=args.plotas)
