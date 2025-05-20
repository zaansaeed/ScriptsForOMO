import os
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
import matplotlib.pyplot as  plt
import pandas as pd
import numpy as np
from sklearn.neural_network import MLPRegressor
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score
from sklearn.metrics import accuracy_score, classification_report, mean_squared_error, r2_score, make_scorer, mean_absolute_error
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.svm import SVR

def plot_results(true_labels_for_testing,y_pred,model):
    plt.plot(true_labels_for_testing, label="Actual", marker='o')
    plt.plot(y_pred, label="Predicted", marker='x')
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
    for item in os.listdir(main_dir):
        working_dir = os.path.join(main_dir, item)
        if os.path.isdir(working_dir):
            os.chdir(working_dir)
            for file in os.listdir(working_dir):
                if file.endswith(f"-{feature}.csv"):
                    data = pd.read_csv(file,header=None,index_col=None)
                    #X.append(pad_square_dataframe_to_array(data,5,0))
                    print(working_dir, data.shape)
                    X.append(data.values.tolist())

    return np.array(X).reshape(len(X),-1)

def create_outputs(main_dir):
    Y =[]
    for item in os.listdir(main_dir):
        if item == "percent6-12-18.txt":
            with open(main_dir+"/percent6-12-18.txt",'r') as f:
                for line in f:
                    row = [float(num) for num in line.split()]
                    Y.append(row)

    return Y

def six_over_target_percents(original_Y):
    new_Y = []
    for row in original_Y:
        new_Y.append(row[0]/(row[0]+row[1]+row[2]))
    return np.array(new_Y)

def create_Y(outputs,cutoff):
    Y =[]
    for item in outputs:
        if item > cutoff:
            Y.append(1)
        else:
            Y.append(0)
    return Y

def run_RFC(X,Y):

    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.4, random_state=42)

    param_grid = {
        'n_estimators': [1, 100, 150, 200, 400],  # Number of trees in the forest
        'max_depth': [10, 30, 40, None],  # Maximum depth of each tree
        'min_samples_split': [2, 5, 10],  # Minimum number of samples required to split an internal node
        'min_samples_leaf': [1, 2, 3],  # Minimum number of samples required to be at a leaf node
        'max_features': ['sqrt', 'log2'],  # Number of features to consider for each split
        'bootstrap': [True, False],  # Whether bootstrap sampling is used
    }
    model = RandomForestClassifier(random_state=42)
    grid_search = GridSearchCV(estimator=model, param_grid=param_grid, n_jobs=-1, cv=4, verbose=1, scoring='accuracy')  # f1
    grid_search.fit(X_train, Y_train)

    y_pred = grid_search.predict(X_test)
    accuracy = accuracy_score(Y_test, y_pred)
    print("Accuracy: ", accuracy)
    print(classification_report(Y_test, y_pred))

def run_RFR(X,Y):
    testing = X[20:40]  # ~5 testing data (PROTEINS R7C5-R8C1) AS TESTING DATA
    true_labels_for_testing = Y[20:40]
    X = np.delete(X, list(range(20, 40)), 0)
    Y = np.delete(Y, list(range(20, 40)), 0)
    #X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.1, random_state=42)

    param_grid = {
        'n_estimators': [200, 400, 600],
        'max_depth': [None, 40, 80],
        'min_samples_split': [2, 3],
        'min_samples_leaf': [1],  # Ensure deep trees are allowed
        'bootstrap': [True],
        'max_features': ['sqrt']
    }
    rf = RandomForestRegressor(random_state=42)
    grid_search = GridSearchCV(
        estimator=rf,
        param_grid=param_grid,
        scoring='neg_mean_absolute_error',
        cv=5,
        n_jobs=-1,
        verbose=1
    )


    rf.fit(X,Y)
    y_pred = rf.predict(testing)
    mse = mean_squared_error(true_labels_for_testing, y_pred)
    print("MSE: ", mse)
    print("MAE,", mean_absolute_error(true_labels_for_testing, y_pred))
    print("Results on testing data:,", y_pred)
    print("True percents", true_labels_for_testing)
    #print("best params:", grid_search.best_params_)
    scores = cross_val_score(grid_search, X, Y, cv=5, scoring='neg_mean_squared_error')
    print(f"Mean cross-validation score: {-scores.mean()}")
    plot_results(true_labels_for_testing, y_pred,"random forest")


def run_SVR(X,Y):
    scaler = StandardScaler()
    X = scaler.fit_transform(X)
    testing = X[20:40]  # ~5 testing data (PROTEINS R7C5-R8C1) AS TESTING DATA
    true_labels_for_testing = Y[20:40]
    X = np.delete(X, list(range(20, 40)), 0)
    Y = np.delete(Y, list(range(20, 40)), 0)
    param_grid = {
        'C': [1, 10, 50, 100],
        'epsilon': [0.01, 0.05, 0.1, 0.2],
        'gamma': [0.01, 0.1, 'scale', 'auto'],
        'kernel': ['rbf','poly','sigmoid'],
    }

    svr = SVR()
    grid_search = GridSearchCV(svr, param_grid, scoring='neg_mean_squared_error', cv=5, n_jobs=-1)
    svr.fit(X, Y)
    #best_model = grid_search.best_estimator_

    y_pred = svr.predict(testing)
    print("predicted: ", y_pred)
    print("Actual: ", true_labels_for_testing)
    print("mse: ", mean_squared_error(true_labels_for_testing, y_pred))
    scores = cross_val_score(svr, X, Y, cv=5, scoring='neg_mean_squared_error')
    print(f"Mean cross-validation score: {-scores.mean()}")
    #print("Best params:", grid_search.best_params_)

    plot_results(true_labels_for_testing, y_pred, "SVR (tuned)")


def run_LR(X,Y):

    testing = X[10:25]  # ~5 testing data (PROTEINS R7C5-R8C1) AS TESTING DATA
    true_labels_for_testing = Y[10:25]
    X = np.delete(X, list(range(10, 25)), 0)
    Y = np.delete(Y, list(range(10, 25)), 0)
    lr = make_pipeline(
        StandardScaler(),
        PolynomialFeatures(degree=3),
        PCA(n_components=15),
        ElasticNet(alpha=0.001, l1_ratio=0.8, max_iter=10000)
    )
    lr.fit(X, Y)
    y_pred = lr.predict(testing)
    print("predicted: ", np.clip(y_pred,0,1))
    print("Actual: ", true_labels_for_testing)
    print("MAE:", mean_absolute_error(true_labels_for_testing, y_pred))
    print("mse: ", mean_squared_error(true_labels_for_testing, y_pred))
    plot_results(true_labels_for_testing, np.clip(y_pred,0,1),"lr")
    scores = cross_val_score(lr, X, Y, cv=5, scoring='neg_mean_squared_error')
    print(f"Mean cross-validation score: {-scores.mean()}")

