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
from sklearn.ensemble import GradientBoostingRegressor
from sklearn.model_selection import LeaveOneOut
from visualization import visualize_model, analyze_feature_ranges,generate_feature_map, visualize_peptide_and_save_features
from sklearn.exceptions import ConvergenceWarning
from sklearn.utils._testing import ignore_warnings
from sklearn.linear_model import ElasticNet
from sklearn.pipeline import Pipeline

from joblib import dump, load

def calc_metrics(Y_test, y_pred):
    mse = mean_squared_error(Y_test, y_pred)
    print("MSE: ", mse)
    print("RMSE: ", np.sqrt(mse))
    print("MAE:", mean_absolute_error(Y_test, y_pred))
    print("R2: ", r2_score(Y_test, y_pred))


def plot_results(true_labels_for_testing, y_pred, model):
    plt.figure(figsize=(10, 6))
    plt.plot(true_labels_for_testing, 'o', label="Actual", markersize=6)
    plt.plot(y_pred, 'x', label="Predicted", markersize=6)


    for i, (true_val, pred_val) in enumerate(zip(true_labels_for_testing, y_pred)):
        plt.plot([i, i], [true_val, pred_val], 'k--', alpha=0.6)  # Vertical dashed line
    all_vals = np.concatenate([true_labels_for_testing, y_pred])
    ymin = all_vals.min() - 0.1
    ymax = all_vals.max() + 0.1
    plt.ylim(ymin, ymax)

    plt.title(f"{model} - Peptides")
    plt.xlabel("Peptide Index")
    plt.ylabel("Target Value")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

def peptide_csv_to_array(main_dir,names,feature):
    X = []
    for name in names:
        working_dir = os.path.join(main_dir, f"Peptide_{name}")
        if os.path.isdir(working_dir):
            os.chdir(working_dir)
            for file in os.listdir(working_dir):  # working in folder
                if file.endswith(f"{feature}.csv"):
                    data = pd.read_csv(file, header=None, index_col=None)
                    if feature == "BWDihedralNormalized":  # remove the last padded 0, then ensure that the boltzmann weighted ~0.9 -> 1
                        data = np.array(data)
                        data = data[:, :-1]
                        data[:, -1] = np.round(data[:, -1])
                        data = pd.DataFrame(data)
                    X.append(data.values.tolist())

    X = np.array(X)
    return X.reshape(len(X), -1)

def create_X(main_dir,names,features): #takes in csv file and reads into array
    X = []
    for feature in features:
        data = peptide_csv_to_array(main_dir,names,feature)
        X.append(data)
    X_new = np.hstack([arr for arr in X])
    return X_new


def sort_by_names_alphabetically(names,values) -> list:
    names_percents_dictionary = dict(zip(names, values))
    names_percents_dictionary = dict((k, names_percents_dictionary[k]) for k in natsorted(names_percents_dictionary))
    Y = []
    for value in names_percents_dictionary.values():
        Y.append(value)
    return Y

def create_Y(percents) -> np.array:
    Y =[]
    for string_percent in percents:
        Y.append(string_percent[0]/(string_percent[0]+string_percent[1]+string_percent[2]))
    return np.array(Y)



def run_RFC(X,Y):
    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.2, random_state=42)
    param_grid = {
        'n_estimators': [50,75],
        'max_depth': [None, 5, 10, 20],
        'min_samples_split': [2, 5, 10],
        'min_samples_leaf': [1, 2, 4],
        'max_features': ['sqrt'],
        'bootstrap': [True],
    }
    rfc = RandomForestClassifier( random_state=42)
    grid_search = GridSearchCV(
        estimator=rfc,
        param_grid=param_grid,
        cv=5,
        scoring='f1_macro',  # or 'f1', 'roc_auc' depending on your goal
        n_jobs=-1,  # use all available cores
        verbose=1
    )

    grid_search.fit(X_train, y_train)

    # Predict on test data
    y_pred = grid_search.predict(X_test)

    # Evaluate
    print("Test accuracy:", accuracy_score(y_test, y_pred))
    print("F1: ", f1_score(y_test, y_pred, average='macro'))
    print(classification_report(y_test, y_pred))
    plot_results(y_test, y_pred, 'rfc')



def true_errors(Y_test, y_pred):
    ranges = {
        "Excellent": 0,
        "Good": 0,
        "Fair": 0,
        "Poor": 0,
    }
    for Y, y in zip(Y_test, y_pred):
        if abs(Y-y) <= 0.1:
            ranges["Excellent"] += 1
        elif abs(Y-y) <= 0.20:
            ranges["Good"] += 1
        elif abs(Y-y) <= 0.3:
            ranges["Fair"] += 1
        else:
            ranges["Poor"] += 1
    return ranges

def compute_weighted_success(true_errors, weights):
    total = sum(true_errors.values())
    score = sum(true_errors[cat] * weights.get(cat, 0) for cat in true_errors)
    return score / total if total > 0 else 0


    # Custom weighted success scorer

def custom_success_metric(y_true, y_pred):
    errors = np.abs(y_pred - y_true)
    score = 0
    for err in errors:
        if err <= 0.1:
            score += 1.0
        elif err <= 0.2:
            score += 0.7
        elif err <= 0.3:
            score += 0.3

        # Poor contributes 0
    return score / len(y_true)


def run_RFR(X, Y, n_splits,test_size):
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_size, random_state=42)
    weighted_success_scorer = make_scorer(custom_success_metric, greater_is_better=True)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

    pipeline = Pipeline([
        ('model', RandomForestRegressor(random_state=42))
    ])

    param_grid = {
        'model__n_estimators': [50, 100, 200, 300, 400, 500, 750],
        'model__max_depth': [None, 5, 10, 20, 30, 40, 50, 75, 100],
        'model__min_samples_split': [2, 5, 10, 20, 50],
        'model__min_samples_leaf': [1, 2, 4, 10, 20],
        'model__max_features': ['sqrt', 'log2', 0.2, 0.4, 0.6, 0.8, 1.0],
        'model__bootstrap': [True],
        'model__max_leaf_nodes': [None, 10, 20, 50, 100],
        'model__criterion': ['squared_error', 'absolute_error', 'friedman_mse', 'poisson'],
    }
    from sklearn.model_selection import RandomizedSearchCV
    search = RandomizedSearchCV(
        estimator=pipeline,  # Your full pipeline with 'model' step
        param_distributions=param_grid,
        n_iter=200,  # Number of parameter settings sampled (adjust for speed)
        scoring='neg_mean_squared_error',
        cv=kf,  # 5-fold CV (adjust if desired)
        random_state=42,
        n_jobs=-1,
        verbose=1,
    )

    search.fit(X_train, Y_train)
    best = search.best_estimator_

    # Outputs
    print("Best params:", search.best_params_)
    #print("mean cv r2:", search.best_score_)
    loo = LeaveOneOut()
    y_pred = []
    y_true = []
    for train_index, test_index in loo.split(X_train):
        x_train, x_test = X_train[train_index], X_train[test_index]
        y_train, y_test = Y_train[train_index], Y_train[test_index]
        best.fit(x_train, y_train)
        y_pred.append(best.predict(x_test))
        y_true.append(y_test)
    calc_metrics(y_true, y_pred)
    plot_results(y_true, y_pred, 'random_forest ')
    #dump(best, 'random_forest.joblib')
    #np.savetxt("X.csv", X, delimiter=",")
   # np.savetxt("y.csv", Y, delimiter=",")


def calculate_cv_scores(kf,X_train,Y_train,best_estimator):
    scores = {
        "Excellent": 0,
        "Good": 0,
        "Fair": 0,
        "Poor": 0,
    }
    mean_CV_metric = 0
    print("Success rates on CV folds")
    for i, (train_index, test_index) in enumerate(kf.split(X_train)):
        x_train, x_test = X_train[train_index], X_train[test_index]
        y_train, y_test = Y_train[train_index], Y_train[test_index]
        best_estimator.fit(x_train, y_train)
        y_pred = best_estimator.predict(x_test)
        y_pred = np.clip(y_pred, 0, 1)
        temp_scores = true_errors(y_test, y_pred)
        scores = {key: scores[key] + temp_scores[key] for key in scores}
        plot_results(y_test, y_pred, i)
        temp_metric=custom_success_metric(y_test, y_pred)
        print(temp_scores,temp_metric )
        mean_CV_metric += temp_metric
    print("mean success rate on cv:", mean_CV_metric/kf.n_splits)

def dummy_RFR(X,Y,X_test,Y_test):

   # X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size= test_size,random_state=42)

    # Dummy model that always predicts the mean of y_train
    dummy = DummyRegressor(strategy='mean')
    dummy.fit(X, Y)

    # Predict and evaluate
    y_pred = dummy.predict(X_test)
    print("dummy regressor: ", custom_success_metric(Y_test, y_pred))
    print("dummy r2:", r2_score(Y_test, y_pred))

def plot_scores_distribution(scores):
    plt.bar(scores.keys(), scores.values())
    plt.xlabel('Category')
    plt.ylabel('Count')
    plt.title('Category Distribution')
    plt.show()


def plot_Y_distribution(Y):

    plt.hist(Y, bins=50, edgecolor='k')
    plt.title('Distribution of Y values')
    plt.xlabel('Y')
    plt.ylabel('Frequency')
    plt.show()

def run_SVR(X, Y, n_splits,test_size):
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_size, random_state=42)
    weighted_success_scorer = make_scorer(custom_success_metric, greater_is_better=True)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

    pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('model', SVR())
    ])

    param_grid = {
        'model__kernel': ['poly','rbf','linear'],
        'model__degree': [2, 3, 4, 5],  # Added degree=2 for comparison
        'model__C': [1, 5, 10, 15, 20, 30],  # Added smaller and larger C values for regularization strength
        'model__gamma': [0.0001, 0.0005, 0.001, 0.002, 0.005, 0.01],  # Broader gamma search
        'model__epsilon': [0.01, 0.05, 0.1, 0.15, 0.2],  # Added smaller epsilon for sensitivity
        'model__coef0': [0.0, 0.3, 0.5, 0.7, 1.0],  # Include zero for no independent term influence
    }
    from sklearn.model_selection import RandomizedSearchCV
    search = RandomizedSearchCV(
        estimator=pipeline,  # Your pipeline with 'model' step as SVR
        param_distributions=param_grid,
        n_iter=500,  # Number of parameter settings to sample
        scoring='r2',
        cv=5,
        random_state=42,
        n_jobs=-1,
        verbose=2,
    )

    search.fit(X_train, Y_train)
    best = search.best_estimator_

    # Outputs
    print("Best params:", search.best_params_)
    cv_scores = cross_val_score(best, X_train, Y_train, cv=kf, scoring='r2')
    for num in cv_scores:
        print(f"test result: {num}")
    print("mean success rate on cv:", cv_scores.mean())

    # Example for feature 0
    y_pred = best.predict(X_test)
    calc_metrics(Y_test, y_pred)
    plot_results(Y_test, y_pred, 'svr ')



def run_NN(X, Y, test_size, n_splits):
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_size, random_state=42)

    weighted_success_scorer = make_scorer(custom_success_metric, greater_is_better=True)

    pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('pca', PCA()),
        ('nn', MLPRegressor(random_state=42, early_stopping=True, validation_fraction=0.1, n_iter_no_change=30))
    ])

    param_grid = {
        'pca__n_components': [0.85, 0.90, 0.95, 0.99],
        'nn__hidden_layer_sizes': [(256,128,64,32), (200,150,100,50)],
        'nn__activation': ['relu', 'tanh'],
        'nn__solver': ['adam'],
        'nn__alpha': [0.0001, 0.001, 0.01],
        'nn__learning_rate_init': [0.001, 0.0005],
        'nn__max_iter': [2000],
        'nn__early_stopping': [True],
        'nn__n_iter_no_change': [30, 50]
    }

    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

    grid_search = GridSearchCV(
        estimator=pipeline,
        param_grid=param_grid,
        scoring='r2',
        cv=kf,
        n_jobs=-1,
        verbose=1
    )

    print("Starting GridSearchCV for automatic Neural Network optimization with PCA...")
    grid_search.fit(X_train, Y_train)
    print("GridSearchCV for Neural Network optimization finished.")

    best_pipeline = grid_search.best_estimator_
    y_pred = best_pipeline.predict(X_test)
    y_pred = np.clip(y_pred, 0, 1)

    print("\nBest params for Neural Network (pipeline) found by GridSearchCV:", grid_search.best_params_)

    print("\nMetrics on Test Set (Optimized Neural Network with PCA):")
    calc_metrics(Y_test, y_pred)

    test_case_errors = true_errors(Y_test, y_pred)
    print("\nTrue Errors (Optimized Neural Network with PCA):")
    print(test_case_errors)

    success_rate_test = custom_success_metric(Y_test, y_pred)
    print(f"\nSuccess rate on test set (Optimized Neural Network with PCA): {success_rate_test:.4f}")

    plot_results(Y_test, y_pred, 'Optimized Neural Network with PCA')

    print("\nNote: Direct feature_importances_ attribute is not available for MLPRegressor.")
    print("Consider using permutation importance or SHAP for feature importance analysis with NNs.")

    print("\nCalculating CV scores for the best Neural Network pipeline with PCA...")
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    calculate_cv_scores(kf, X_train, Y_train, best_pipeline)



def run_GBR(X, Y, test_size, n_splits):
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_size, random_state=42)
    weighted_success_scorer = make_scorer(custom_success_metric, greater_is_better=True)
    from xgboost import XGBRegressor

    pipeline = Pipeline([
        ('model', XGBRegressor(objective='reg:squarederror', random_state=42, verbosity=0))
    ])

    param_grid = {
        'model__n_estimators': [50, 100, 200, 300, 500, 750, 1000],
        'model__learning_rate': [0.001, 0.01, 0.05, 0.1, 0.2, 0.3],
        'model__max_depth': [2, 3, 4, 5, 7, 10, 15],
        'model__min_child_weight': [1, 3, 5, 10],
        'model__subsample': [0.4, 0.6, 0.8, 1.0],
        'model__colsample_bytree': [0.4, 0.6, 0.8, 1.0],
        'model__gamma': [0, 0.1, 0.3, 0.5, 1.0],  # For pruning
        'model__reg_alpha': [0.0, 0.001, 0.01, 0.1],  # L1 regularization
        'model__reg_lambda': [0.0, 0.01, 0.1, 1.0],  # L2 regularization
    }
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

    grid_search = RandomizedSearchCV(
        estimator=pipeline,
        param_distributions=param_grid,
        scoring='neg_mean_squared_error',
        cv=kf,
        n_iter=200,
        n_jobs=-1,
        verbose=1
    )

    grid_search.fit(X_train, Y_train)
    best = grid_search.best_estimator_

    loo = LeaveOneOut()
    y_pred = []
    y_true = []
    for train_index, test_index in loo.split(X_train):
        x_train, x_test = X_train[train_index], X_train[test_index]
        y_train, y_test = Y_train[train_index], Y_train[test_index]
        best.fit(x_train, y_train)
        y_pred.append(best.predict(x_test))
        y_true.append(y_test)
    print("r2 score:", r2_score(y_true, y_pred))
    test_case = true_errors(Y_test, y_pred)
    print("\n")
    print(test_case)

    plot_results(y_true, y_pred, 'gradient boosting')
    #visualize_model(best, X, Y)
    #calculate_cv_scores(kf, X_train, Y_train, best)






@ignore_warnings(category=ConvergenceWarning)
def run_elasticnet(X, Y, n_splits,test_size):
    # Split data
    # Custom scorer, replace with your own function if needed
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)
    # Pipeline: scaling + ElasticNet
    from sklearn.preprocessing import StandardScaler, PolynomialFeatures
    pipeline = Pipeline([
        ('model', ElasticNet(max_iter=10000))
    ])
    loo = LeaveOneOut()
    param_grid = {
        # Polynomial expansion (non-linear interactions)
        # ElasticNet regularization strength
        'model__alpha': [
            1e-8, 1e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3,
            0.01, 0.025, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5,
        ],

        # ElasticNet L1/L2 mixing - more granular values
        'model__l1_ratio': [
            0.0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45,
        ],

        # Convergence tolerance - expanded with smaller values
        'model__tol': [
            1e-8, 1e-7, 1e-6, 5e-6, 1e-5, 5e-5, 1e-4, 5e-4, 1e-3, 5e-3, 1e-2
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

    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)
    # GridSearchCV setup
    from sklearn.model_selection import RandomizedSearchCV

    search = RandomizedSearchCV(
        estimator=pipeline,  # Your pipeline with 'model' step as ElasticNet
        param_distributions=param_grid,
        scoring='neg_mean_squared_error',
        cv=kf,
        n_iter=200,
        n_jobs=-1,
        verbose=1,
    )
    # Fit model and search best params
    search.fit(X_train, Y_train)
    best = search.best_estimator_
    print("Best params:", search.best_params_)
    print("best score:", search.best_score_)
    # Outputs

    y_pred = []
    y_true = []
    for train_index, test_index in loo.split(X_train):
        x_train, x_test = X_train[train_index], X_train[test_index]
        y_train, y_test = Y_train[train_index], Y_train[test_index]
        best.fit(x_train, y_train)
        y_pred.append(best.predict(x_test))
        y_true.append(y_test)
    calc_metrics(y_true, y_pred)
    plot_results(y_true, y_pred, 'elastic net ')
   # dump(best, 'elasticnet_model.joblib')
    #np.savetxt("X.csv", X, delimiter=",")
    #np.savetxt("y.csv", Y, delimiter=",")


def create_Y_ROG(main_dir,names):
    Y= []
    for name in names:
        working_dir = os.path.join(main_dir, f"Peptide_{name}")
        if os.path.isdir(working_dir):
            os.chdir(working_dir)
            for file in os.listdir(working_dir):  # working in folder
                if file.endswith("RadiusOfGyration.csv"):
                    data = pd.read_csv(file, header=None, index_col=None)
                    Y.append(data.values.tolist()[0][0])
    return np.array(Y)