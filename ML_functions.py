import os
from natsort import natsorted
import matplotlib.pyplot as  plt
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn.model_selection import GridSearchCV ,cross_val_score, KFold
from sklearn.metrics import accuracy_score, classification_report, mean_squared_error, \
     f1_score,mean_absolute_error, r2_score
from sklearn.svm import SVR
from sklearn.model_selection import train_test_split
from sklearn.metrics import make_scorer
from sklearn.dummy import DummyRegressor
from sklearn.decomposition import PCA
from sklearn.ensemble import GradientBoostingRegressor




def calc_metrics(Y_test, y_pred):
    mse = mean_squared_error(Y_test, y_pred)
    print("MSE: ", mse)
    print("MAE,", mean_absolute_error(Y_test, y_pred))
    print("R2: ", r2_score(Y_test, y_pred))


def plot_results(true_labels_for_testing, y_pred, model):
    plt.figure(figsize=(10, 6))
    plt.plot(true_labels_for_testing, 'o', label="Actual", markersize=6)
    plt.plot(y_pred, 'x', label="Predicted", markersize=6)


    for i, (true_val, pred_val) in enumerate(zip(true_labels_for_testing, y_pred)):
        plt.plot([i, i], [true_val, pred_val], 'k--', alpha=0.6)  # Vertical dashed line

    plt.title(f"{model} - Peptides")
    plt.xlabel("Peptide Index")
    plt.ylabel("Target Value")
    plt.ylim(0, 1)
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

def plot_importances(best_estimator,X,top_n):
    importances = best_estimator.feature_importances_
    indices = np.argsort(importances)[::-1][:top_n]

    # Plot
    Xf = pd.DataFrame(X, columns=[f"Feature {i}" for i in range(X.shape[1])])
    plt.figure(figsize=(10, 6))
    plt.title("Feature Importances - Random Forest Regressor")
    plt.bar(range(top_n), importances[indices], align="center")
    plt.xticks(range(top_n), Xf.columns[indices], rotation=45)
    plt.ylabel("Importance")
    plt.tight_layout()
    plt.show()

def run_RFR(X,Y,test_size,n_splits):


    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_size,random_state=42)
    weighted_success_scorer = make_scorer(custom_success_metric,greater_is_better=True)

    param_grid = {
        'n_estimators': [25,75,50],
        'max_depth': [3,6,7,8,9],
        'max_features': ['sqrt'],
        'min_samples_split': [2,3,4,5],
        'min_samples_leaf': [1,2,3,4],
        'bootstrap': [True],

    }
    rf = RandomForestRegressor(random_state=42)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

    grid_search = GridSearchCV(
        estimator=rf,
        param_grid=param_grid,
        scoring=weighted_success_scorer,
        cv=kf,
        n_jobs=-1,
        verbose=1
    )


    grid_search.fit(X_train, Y_train)
    best = grid_search.best_estimator_

    y_pred = best.predict(X_test)
    print("Best params:", grid_search.best_params_)
    calc_metrics(Y_test, y_pred)
    test_case = true_errors(Y_test,y_pred)
    print("\n")
    print(test_case)
    print("success rate on test:",custom_success_metric(Y_test, y_pred))

    plot_results(Y_test, y_pred, 'random forest')
    plot_importances(best_estimator=best,X=X,top_n=10)

    calculate_cv_scores(kf,X_train,Y_train,best)

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
        temp_scores = true_errors(y_test, y_pred)
        scores = {key: scores[key] + temp_scores[key] for key in scores}
        plot_results(y_test, y_pred, i)
        temp_metric=custom_success_metric(y_test, y_pred)
        print(temp_scores,temp_metric )
        mean_CV_metric += temp_metric
    print("mean success rate on cv:", mean_CV_metric/kf.n_splits)

def dummy_RFR(X,Y,test_size):

    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size= test_size,random_state=42)

    # Dummy model that always predicts the mean of y_train
    dummy = DummyRegressor(strategy='mean')
    dummy.fit(X_train, y_train)

    # Predict and evaluate
    y_pred = dummy.predict(X_test)
    print("dummy regressor: ", custom_success_metric(y_test, y_pred))

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

def run_SVR(X,Y,test_size,n_splits):


    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=test_size,random_state=42)
    weighted_success_scorer = make_scorer(custom_success_metric,greater_is_better=True)
    kf = KFold(n_splits=n_splits, shuffle=True, random_state=42)

    param_grid = {
        'kernel': ['rbf','poly'],
        'C': [100,5,10,1],
        'epsilon': [ .1,.05],
    }

    svr = SVR()
    grid_search = GridSearchCV(svr, param_grid, scoring='r2', cv=kf, n_jobs=-1)
    grid_search.fit(X_train, Y_train)
    best = grid_search.best_estimator_

    y_pred = best.predict(X_test)
    print("Best params:", grid_search.best_params_)
    calc_metrics(Y_test, y_pred)
    test_case = true_errors(Y_test, y_pred)
    print("\n")
    print(test_case)
    print("success rate on test:", custom_success_metric(Y_test, y_pred))

    plot_results(Y_test, y_pred, 'svr')
    #plot_importances(best_estimator=best, X=X, top_n=10)

    calculate_cv_scores(kf, X_train, Y_train, best)






