import os
from natsort import natsorted
import matplotlib.pyplot as  plt
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import GridSearchCV, cross_val_score, KFold, RandomizedSearchCV
from sklearn.metrics import accuracy_score, classification_report, mean_squared_error, \
     f1_score,mean_absolute_error, r2_score
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVR
from sklearn.model_selection import train_test_split
from sklearn.metrics import make_scorer
from sklearn.dummy import DummyRegressor
from sklearn.model_selection import LeaveOneOut
from sklearn.exceptions import ConvergenceWarning
from sklearn.utils._testing import ignore_warnings
from sklearn.linear_model import ElasticNet
from sklearn.pipeline import Pipeline
import glob
import logging
from joblib import dump, load
import streamlit as st
from datetime import datetime
from pathlib import Path




config = None
def init_config(cfg):
    global config
    config = cfg
    return config
ml_logger = logging.getLogger("ml")


def gaussian_scoring(y_true, y_pred):

    y_true = np.array(y_true)
    y_pred = np.array(y_pred)
    errors = np.abs(y_pred - y_true)

    scores = np.zeros_like(errors)

    # Piecewise for errors: 0 <= err <= 30
    mask1 = (errors >= 0) & (errors <= 30)
    scores[mask1] = 1 + 0.047 * (1 - np.exp(0.1 * errors[mask1]))

    # Piecewise for errors: 30 < err <= 100
    mask2 = (errors > 30) & (errors <= 100)
    mu = 16.08
    sigma = 4.21
    scores[mask2] = 24.36 * np.exp(-((errors[mask2] - mu) ** 2) / (2 * sigma ** 2)) + 0.000001
    return np.mean(scores)

def new_adjusted_metrics(y_true,y_pred):
    errors = np.abs(y_pred - y_true)
    scores = gaussian_scoring(y_true, y_pred)
    weighted_scores = [error*(1-weight) for error,weight in zip(errors,scores)]
    weighted_scores = np.array(weighted_scores)
    dict_metrics = {
        "NEW MSE": np.mean(weighted_scores**2),
        "NEW RMSE": np.sqrt(np.mean(weighted_scores**2)),
        "NEW MAE": np.mean(np.abs(weighted_scores)),
    }
    return dict_metrics


def calc_metrics(Y_test, y_pred):
    mse = mean_squared_error(Y_test, y_pred)
    dict_metrics = {
        "CUSTOM SCORE FUNCTION": gaussian_scoring(Y_test, y_pred).mean(),
        "MSE": mse,
        "RMSE": np.sqrt(mse),
        "MAE": mean_absolute_error(Y_test, y_pred),
        "R2": r2_score(Y_test, y_pred)
    }
    return dict_metrics


def plot_results(true_labels_for_testing, y_pred, model,
                 out_dir: str = "plots") -> str:
    this_file = Path(__file__).resolve()
    out_dir = os.path.join(this_file.parent, out_dir)
    """Build the scatter/line plot, save it, and return the file path."""
    # ---------- build the figure ----------
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(true_labels_for_testing, "o", label="Actual",  markersize=6)
    ax.plot(y_pred,               "x", label="Predicted", markersize=6)

    for i, (t, p) in enumerate(zip(true_labels_for_testing, y_pred)):
        ax.plot([i, i], [t, p], "k--", alpha=0.6)  # vertical dashed line

    all_vals = np.concatenate([true_labels_for_testing, y_pred])
    buffer = 0.3  # adjust the margin as needed
    y_min = all_vals.min() - buffer
    y_max = all_vals.max() + buffer
    ax.set_ylim(y_min, y_max)

    ax.set(title=f"{model} - Peptides",
           xlabel="Peptide index",
           ylabel="Target value")
    ax.legend()
    ax.grid(True)
    fig.tight_layout()

    # ---------- save ----------
    os.makedirs(out_dir, exist_ok=True)
    timestamp   = datetime.now().strftime("%Y%m%d-%H%M%S")
    file_name   = f"{model.replace(' ', '_')}_{timestamp}.png"
    file_path   = os.path.join(out_dir, file_name)
    fig.savefig(file_path, dpi=300, bbox_inches="tight", transparent=False)
    plt.close(fig)   

def filter_names(main_dictionary) -> dict:
    best_entries = {}

    for name, (smiles, percent) in main_dictionary.items():
        if smiles not in best_entries or percent > best_entries[smiles][1]:
            best_entries[smiles] = (name, percent)

    # Step 2: rebuild dictionary with the chosen names
    filtered_dict = {}
    for smiles, (name, percent) in best_entries.items():
        filtered_dict[name] = (smiles, percent)

    filtered_dict = dict(natsorted(filtered_dict.items()))
    return filtered_dict


def peptide_csv_to_array(name, feature, main_dir) -> np.ndarray:
    folder = os.path.join(main_dir, f"Peptide_{name}")
    pattern = os.path.join(folder, f"*{feature}.csv")
    matches = glob.glob(pattern)

    feature_file = matches[0]

    data = pd.read_csv(feature_file, header=None, index_col=None)
    if feature == "BWDihedralNormalized":
        #data = data.drop(columns=[2,3])
        data = np.array(data)
        data = data[:, :-1]
        data[:, -1] = np.round(data[:, -1])
        data = pd.DataFrame(data)

    return data.values.flatten()


def create_model_data(names,features, main_dir,target_value):
    X_all = []
    Y = []

    for folder in natsorted(os.listdir(main_dir)):
        if folder.startswith("Peptide_"):
            name = folder.split("_")[1]
            if name in names:
                working_dir = os.path.join(main_dir, folder)
                os.chdir(working_dir)
                # Load target
                target_file = os.path.join(working_dir, f"{name}_{target_value}.txt")
                with open(target_file) as f:
                    for line in f:
                        Y.append(float(line.split()[0]))

                X_temp = []
                print(name)
                for feature in features:
                    x_feat = peptide_csv_to_array(name, feature, main_dir)
                    X_temp.append(x_feat)

                # Stack features into single row for this peptide
                X_all.append(np.concatenate(X_temp))


            else:
                pass
    X = np.array(X_all)
    Y = np.array(Y)
    return X, Y


def true_errors(Y_test, y_pred):
    ranges = {
        "Excellent": 0,
        "Good": 0,
        "Fair": 0,
        "Poor": 0,
    }
    for Y, y in zip(Y_test, y_pred):
        if abs(Y-y) <= 10:
            ranges["Excellent"] += 1
        elif abs(Y-y) <= 20:
            ranges["Good"] += 1
        elif abs(Y-y) <= 30:
            ranges["Fair"] += 1
        else:
            ranges["Poor"] += 1
    return ranges


def dummy(X,Y):

   # X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size= test_size,random_state=42)
    # Dummy model that always predicts the mean of y_train
    dummy = DummyRegressor(strategy='mean')
    dummy.fit(X,Y)
    y_pred = dummy.predict(X)
    y_true = Y

    #scores = cross_val_score(dummy_model, X, y, cv=kf, scoring="accuracy")  # Example: accuracy for classification problems

    # Predict and evaluate
  
    print("dummy r2:", r2_score(y_true, y_pred))
    print("dummy mse:", mean_squared_error(y_true, y_pred))
    print("dummy rmse:", np.sqrt(mean_squared_error(y_true, y_pred)))
    print("dummy mae:", mean_absolute_error(y_true, y_pred))
    print("dummy gaussian score", gaussian_scoring(Y, y_pred))
    #print("DUMMY NEW METRICS:")
    #print(new_adjusted_metrics(y_true, y_pred))



def plot_Y_distribution(Y):
    plt.hist(Y, bins=50, edgecolor='k')
    plt.title('Distribution of Y values')
    plt.xlabel('Y')
    plt.ylabel('Frequency')
    plt.show()


def run_RFR(X, Y):
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)
    kf = KFold(n_splits=5, shuffle=True, random_state=42)
    pipeline = Pipeline([
        ('model', RandomForestRegressor(random_state=42))
    ])
    loo = LeaveOneOut()

    param_grid = {
        'model__n_estimators': [ 400,300,350,450],
        'model__max_depth': [None,15,20,22,23,19    ],
        'model__min_samples_split': [2, 5, 3,4],
        'model__min_samples_leaf': [1, 2, 4, 10, 20],
        'model__max_features': ['sqrt', 'log2', 0.2, 0.4, 0.6, 0.8, 1.0],
        'model__bootstrap': [True],
        'model__max_leaf_nodes': [None, 10, 20, 50, 100],
        'model__criterion': ['squared_error', 'absolute_error'],
    }
    search = RandomizedSearchCV(
        estimator=pipeline,  # Your pipeline with 'model' step as SVR
        param_distributions=param_grid,
        scoring='neg_mean_absolute_error',
        cv=kf,
        n_iter=config["machine_learning"]["n_iter"],
        n_jobs=-1,
        verbose=1,
    )
    # Fit model and search best params
    search.fit(X_train, Y_train)
    best = search.best_estimator_

    y_true = Y_test
    y_pred = best.predict(X_test)
    #y_pred = np.clip(y_pred, 0,100)
    metrics = calc_metrics(y_true, y_pred)
    
    params_str = ",".join(f"{key}={value}" for key, value in search.best_params_.items())
    line = f"{params_str},best_score={search.best_score_:.4f}"
    ml_logger.info(line)
    ml_logger.info("RFR Results:")
    metric_line = ",".join(f"{key}={value:.4f}" for key, value in metrics.items())
    ml_logger.info(metric_line)

    plot_results(y_true, y_pred, f'RFR on {config["machine_learning"]["features_to_train_on"]}')
    os.chdir(config["data_generation"]["main_dir"])
    if config["machine_learning"]["save_model"]:
        dump(best, 'RFR_model.joblib')
        np.savetxt("X.csv", X, delimiter=",")
        np.savetxt("y.csv", Y, delimiter=",")



def run_SVR(X, Y):
    custom_scorer = make_scorer(gaussian_scoring, greater_is_better=True)
    loo = LeaveOneOut()
    pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('model', SVR())
    ])
    kf = KFold(n_splits=10, shuffle=True, random_state=42)
    param_grid = {
        'model__kernel': ['poly', 'rbf', 'sigmoid'],
        'model__degree': [2, 3],  # Added degree=2 for comparison
        'model__C': [1, 5, 10, 15, 20, 30],  # Added smaller and larger C values for regularization strength
        'model__gamma': [0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01],  # Broader gamma search
        'model__epsilon': [0.01, 0.05, 0.1, 0.15, 0.2],  # Added smaller epsilon for sensitivity
        'model__coef0': [0.0, 0.3, 0.5, 0.7, 1.0],  # Include zero for no independent term influence
    }

    search = RandomizedSearchCV(
        estimator=pipeline,  # Your pipeline with 'model' step as SVR
        param_distributions=param_grid,
        scoring='r2',
        cv=kf,
        n_iter=config["machine_learning"]["n_iter"],
        n_jobs=-1,
        verbose=1,
    )
    # Fit model and search best params
    search.fit(X, Y)
    best = search.best_estimator_

    loo = LeaveOneOut()

    y_pred = []
    y_true = []
    for train_index, test_index in kf.split(X):
        x_train, x_test = X[train_index], X[test_index]
        y_train, y_test = Y[train_index], Y[test_index]
        best.fit(x_train, y_train)
        y_pred.extend(best.predict(x_test))
        y_true.extend(y_test)
    #y_pred = np.clip(y_pred, 0,100)
    metrics = calc_metrics(y_true, y_pred)
    print(metrics)
    params_str = ",".join(f"{key}={value}" for key, value in search.best_params_.items())
    line = f"{params_str},best_score={search.best_score_:.4f}"
    ml_logger.info(line)
    ml_logger.info("SVR Results:")
    metric_line = ",".join(f"{key}={value:.4f}" for key, value in metrics.items())
    ml_logger.info(metric_line)

    plot_results(y_true, y_pred, f'SVR on {config["machine_learning"]["features_to_train_on"]}')
    os.chdir(config["data_generation"]["main_dir"])
    if config["machine_learning"]["save_model"]:
        dump(best, 'SVR_model.joblib')
        np.savetxt("X.csv", X, delimiter=",")
        np.savetxt("y.csv", Y, delimiter=",")



@ignore_warnings(category=ConvergenceWarning)
def run_ElasticNet(X, Y):
    custom_scorer = make_scorer(gaussian_scoring, greater_is_better=True)
    kf = KFold(n_splits=10, shuffle=True, random_state=42)
    from sklearn.preprocessing import StandardScaler, PolynomialFeatures
    from sklearn.decomposition import PCA

    pipeline = Pipeline([
       # ('poly', PolynomialFeatures(degree=3, include_bias=False)),
       # ('pca', PCA(n_components=0.9)),  # Keep 95% of variance
        ('model', ElasticNet(max_iter=10000))
    ])
    loo = LeaveOneOut()
    param_grid = {
        # Polynomial expansion (non-linear interactions)
        # ElasticNet regularization strength
        'model__alpha': [
            5e-5, 1e-4, 5e-4, 1e-3, 5e-3,
            0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.21,0.19, 0.18,0.17, 0.25, 0.3, 0.4, 0.5,
        ],

        # ElasticNet L1/L2 mixing - more granular values
        'model__l1_ratio': [
            0.0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        ],

        # Convergence tolerance - expanded with smaller values
        'model__tol': [
            0.0001,0.0005,0.0003
        ],

        # Max iterations - significantly increased to address convergence
        'model__max_iter': [
            1000, 2000, 5000, 10000, 20000, 50000, 100000
        ],

        # Optimization strategy
        'model__selection': ['cyclic', 'random'],

        # Whether to fit the intercept term
        'model__fit_intercept': [True, False]

    }

    # GridSearchCV setup

    search = RandomizedSearchCV(
        estimator=pipeline,  # Your pipeline with 'model' step as ElasticNet
        param_distributions=param_grid,
        scoring='r2',
        cv=kf,
        n_iter=config["machine_learning"]["n_iter"],
        n_jobs=-1,
        verbose=1,
    )
    # Fit model and search best params
    search.fit(X, Y)
    best = search.best_estimator_


    y_pred = []
    y_true = []
    for train_index, test_index in kf.split(X):
        x_train, x_test = X[train_index], X[test_index]
        y_train, y_test = Y[train_index], Y[test_index]
        best.fit(x_train, y_train)
        pred = best.predict(x_test)
        y_pred.extend(pred)
        y_true.extend(y_test)

    #y_pred = np.clip(y_pred, 0,100)
    metrics = calc_metrics(y_true, y_pred)
    dummy(X, Y)
    #print("ELASTICNET NEW METRICS:")
    #print(new_adjusted_metrics(y_true, y_pred))
    print("ELASTICNET OLD METRICS:")
    print(metrics)
    plot_results(y_true, y_pred, f'ElasticNet on {config["machine_learning"]["features_to_train_on"]}')

    params_str = ",".join(f"{key}={value}" for key, value in search.best_params_.items())
    line = f"{params_str},best_score={search.best_score_:.4f}"
    ml_logger.info(line)
    ml_logger.info("ElasticNet Results:")
    metric_line = ",".join(f"{key}={value:.4f}" for key, value in metrics.items())
    ml_logger.info(metric_line)

    os.chdir(config["data_generation"]["main_dir"])

    if config["machine_learning"]["save_model"]:
        dump(best, 'elasticnet_model.joblib')
        np.savetxt("X.csv", X, delimiter=",")
        np.savetxt("y.csv", Y, delimiter=",")

