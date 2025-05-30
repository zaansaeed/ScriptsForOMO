import os
import math
from natsort import natsorted
import matplotlib.pyplot as  plt
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, GradientBoostingRegressor, HistGradientBoostingRegressor
from sklearn.model_selection import GridSearchCV ,cross_val_score, KFold
from sklearn.metrics import accuracy_score, classification_report, mean_squared_error, r2_score, make_scorer, \
     f1_score,mean_absolute_error, r2_score, root_mean_squared_error
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.svm import SVR, SVC
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.metrics import make_scorer



def calc_metrics(Y_test, y_pred):
    mse = mean_squared_error(Y_test, y_pred)
    print("MSE: ", mse)
    print("MAE,", mean_absolute_error(Y_test, y_pred))
    print("Results on testing data:,", y_pred)
    print("True percents", Y_test)


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


def pad_square_dataframe_to_array(df,target_size,fill_value):
    size = df.shape[0]
    original = df.to_numpy()
    padded = np.full((target_size,target_size),fill_value)
    padded[:size,:size] = original
    return pd.DataFrame(padded).values.tolist()

def create_X(main_dir,feature): #takes in csv file and reads into array
    X= []
    for item in natsorted(os.listdir(main_dir)):
        working_dir = os.path.join(main_dir, item)
        if os.path.isdir(working_dir):
            os.chdir(working_dir)
            for file in os.listdir(working_dir):
                if file.endswith(f"{feature}.csv"):
                    data = pd.read_csv(file,header=None,index_col=None)
                    if feature == "BWDihedralNormalized": #remove the last padded 0, then ensure that the boltzmann weighted ~0.9 -> 1
                        data = np.array(data)
                        data = data[:,:-1]

                        data[:,-1] = np.round(data[:,-1])
                        data = pd.DataFrame(data)

                    X.append(data.values.tolist())

    X= np.array(X)
    return X.reshape(len(X),-1)

def create_Y(main_dir):
    Y =[]
    for item in os.listdir(main_dir):
        if item == "percent6-12-18.txt":
            with open(main_dir+"/percent6-12-18.txt",'r') as f:
                for line in f:
                    row = [float(num) for num in line.split()]
                    Y.append(row[0] / (row[0] + row[1] + row[2]))

    return Y




def create_YC(outputs,cutoff):
    Y =[]
    for item in outputs:
        if item >=cutoff:
            Y.append(1)
        else:
            Y.append(0)
    return Y




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
    accuracy = accuracy_score(y_test, y_pred)
    print("Accuracy: ", accuracy)
    print("F1: ", f1_score(y_test, y_pred, average='macro'))
    print(classification_report(y_test, y_pred))
    plot_results(y_test, y_pred, 'rfc')
    print(grid_search.best_params_)
    print(cross_val_score(grid_search.best_estimator_, X, Y, cv=5, scoring='f1_macro').mean())

def run_SVM(X,Y):
    scaler = StandardScaler()
    X = scaler.fit_transform(X)
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)
    param_grid = {
        'C': [1],
        'kernel': ['rbf'],
        'gamma': ['scale', 'auto'],
        'class_weight': ['balanced']
    }

    # Initialize SVM
    svm = SVC(probability=True, random_state=42)

    # Set up GridSearchCV
    grid_search = GridSearchCV(
        estimator=svm,
        param_grid=param_grid,
        scoring='f1',
        cv=5,
        n_jobs=-1,
        verbose=1
    )
    svm.fit(X_train, Y_train)
    y_pred = svm.predict(X_test)



    accuracy = accuracy_score(Y_test, y_pred)
    print("Accuracy: ", accuracy)
    print("F1: ", f1_score(Y_test, y_pred))
    print(classification_report(Y_test, y_pred))
    #print(grid_search.best_params_)
    print(cross_val_score(svm, X, Y, cv=5, scoring='f1_macro').mean())

    plot_results(Y_test, y_pred, grid_search)

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

def compute_weighted_success(category_counts, weights):
    total = sum(category_counts.values())
    score = sum(category_counts[cat] * weights.get(cat, 0) for cat in category_counts)
    return score / total if total > 0 else 0


    # Custom weighted success scorer
def custom_success_metric(y_true, y_pred):
    errors = np.abs((y_pred - y_true) / y_true) * 100
    score = 0
    for err in errors:
        if err <= 10:
            score += 1.0
        elif err <= 20:
            score += 0.7
        elif err <= 30:
            score += 0.3
        # Poor contributes 0
    return score / len(y_true)

def run_RFR(X,Y):

    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.3,random_state=42)
    weighted_success_scorer = make_scorer(custom_success_metric,greater_is_better=True)

    param_grid = {
        'n_estimators': [50,75],
        'max_depth': [8,9,10],
        'max_features': ['sqrt'],
        'min_samples_split': [2,3,4],
        'min_samples_leaf': [1, 2,3,4],
        'bootstrap': [True],
    }
    rf = RandomForestRegressor(random_state=42)
    kf = KFold(n_splits=5, shuffle=True, random_state=42)

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
    scores = cross_val_score(best, X_train, Y_train, cv=kf, scoring='r2')
    print(scores)
    print(f"Mean cross-validation score (Average R2 across 5 cv): {scores.mean()}")
    print("R2: ", r2_score(Y_test, y_pred))
    print("Best params:", grid_search.best_params_)
    print("Mean Squared Error: ", mean_squared_error(Y_test, y_pred))
    print("Root Mean Squared Error:", np.sqrt(mean_squared_error(Y_test, y_pred)))
    print("Mean Absolute Error,", mean_absolute_error(Y_test, y_pred))
    test_case = true_errors(Y_test,y_pred)
    print("\n")
    print(test_case)
    weights = {
        "Excellent": 1.0,
        "Good": 0.7,
        "Fair": 0.3,
        "Poor": 0.0
    }
    print("success rate:",compute_weighted_success(test_case,weights))

    plot_results(Y_test, y_pred, 'random forest')

    top_n= 10
    importances = grid_search.best_estimator_.feature_importances_
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
    scores = {
        "Excellent": 0,
        "Good": 0,
        "Fair": 0,
        "Poor": 0,
    }
    for i, (train_index, test_index) in enumerate(kf.split(X_train)):
        x_train, x_test = X_train[train_index], X_train[test_index]
        y_train, y_test = Y_train[train_index], Y_train[test_index]
        best.fit(x_train, y_train)
        y_pred = best.predict(x_test)
        temp_scores = true_errors(y_test,y_pred)
        scores = {key: scores[key]+temp_scores[key] for key in scores}
        plot_results(y_test, y_pred, i)
    print(scores)

    plt.bar(scores.keys(), scores.values())
    plt.xlabel('Category')
    plt.ylabel('Count')
    plt.title('Category Distribution')
    plt.show()

def run_SVR(X,Y):

    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2,random_state=41)

    pipeline = Pipeline([
        #('scaler', StandardScaler()),
        #('pca', PCA(n_components=20)),  # Let GridSearch decide n_components
        ('svr', SVR())
    ])

    param_grid = {
        'svr__kernel': ['rbf'],
        'svr__C': [100],
        'svr__epsilon': [ .1],
    }



    svr = SVR(C=1,epsilon=.1,kernel='rbf')
    grid_search = GridSearchCV(pipeline, param_grid, scoring='r2', cv=5, n_jobs=-1)
    grid_search.fit(X_train, Y_train)
    best_model = grid_search.best_estimator_

    y_pred = best_model.predict(X_test)
    y_pred = np.clip(y_pred, 0, 1)


    scores = cross_val_score(best_model, X_train, Y_train, cv=5, scoring='r2')
    print(f"Mean cross-validation score (Average R2 across 5 cv): {scores.mean()}")
    print("R2: ", r2_score(Y_test, y_pred))
    #print("Best params:", grid_search.best_params_)
    print("Mean Squared Error: ", mean_squared_error(Y_test, y_pred))
    print("Root Mean Squared Error:", np.sqrt(mean_squared_error(Y_test, y_pred)))
    print("Mean Absolute Error,", mean_absolute_error(Y_test, y_pred))


    plot_results(Y_test, y_pred, svr)


