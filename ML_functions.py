import os
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet
import matplotlib.pyplot as  plt
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, GradientBoostingRegressor, HistGradientBoostingRegressor
from sklearn.model_selection import train_test_split, GridSearchCV, cross_val_score, KFold
from sklearn.metrics import accuracy_score, classification_report, mean_squared_error, r2_score, make_scorer, mean_absolute_error
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.svm import SVR
from xgboost import XGBRegressor
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, ConstantKernel as C
from sklearn.gaussian_process.kernels import RBF, ConstantKernel, WhiteKernel
from sklearn.gaussian_process.kernels import RationalQuadratic

def calc_metrics(Y_test, y_pred):
    mse = mean_squared_error(Y_test, y_pred)
    print("MSE: ", mse)
    print("MAE,", mean_absolute_error(Y_test, y_pred))
    print("Results on testing data:,", y_pred)
    print("True percents", Y_test)

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
                    if feature == "BWDihedralNormalized": #remove the last padded 0, then ensure that the boltzmann weighted ~0.9 -> 1
                        data = np.array(data)
                        data = data[:,:-1]

                        data[:,-1] = np.round(data[:,-1])
                        print(data.shape,working_dir)
                        data = pd.DataFrame(data)

                    X.append(data.values.tolist())

    X= np.array(X)
    return X.reshape(len(X),-1)

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

    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

    param_grid = {
        'n_estimators': [ 200,600],  # Fewer trees to keep it lightweight
        'max_depth': [20, None],  # Shallow trees to prevent overfitting
        'min_samples_split': [2, 10],  # Try larger splits to regularize
        'min_samples_leaf': [2, 4],  # Larger leaves reduce model complexity
        'max_features': ['sqrt'],  # Limit number of features at each split
        'bootstrap': [True]  # Usually better for small data
    }
    rf = RandomForestRegressor(random_state=42)
    grid_search = GridSearchCV(
        estimator=rf,
        param_grid=param_grid,
        scoring='neg_mean_absolute_error',
        cv=8,
        n_jobs=-1,
        verbose=1
    )


    rf.fit(X_train,Y_train)
    y_pred = rf.predict(X_test)
    mse = mean_squared_error(Y_test, y_pred)
    print("MSE: ", mse)
    print("MAE,", mean_absolute_error(Y_test, y_pred))
    print("Results on testing data:,", y_pred)
    print("True percents", Y_test)
    #print("best params:", grid_search.best_params_)
    scores = cross_val_score(rf, X_train, Y_train, cv=5, scoring='neg_mean_squared_error')
    print(f"Mean cross-validation score: {-scores.mean()}")
    plot_results(Y_test, y_pred,"random forest")


def run_SVR(X,Y):
    scaler = StandardScaler()
    X = scaler.fit_transform(X)
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.2, random_state=42)

    param_grid = {
        'C': [50, 100,200],
        'epsilon': [0.01, 0.05, 0.1,0.2],
        'gamma': [0.01, 0.1, 'scale', 'auto'],
        'kernel': ['rbf'],
    }

    svr = SVR()
    grid_search = GridSearchCV(svr, param_grid, scoring='neg_mean_squared_error', cv=5, n_jobs=-1)
    svr.fit(X_train, Y_train)
    #best_model = svr.best_estimator_

    y_pred = svr.predict(X_test)
    print("predicted: ", y_pred)
    print("Actual: ", Y_test)
    print("mse: ", mean_squared_error(Y_test, y_pred))
    scores = cross_val_score(svr, X_train, Y_train, cv=5, scoring='neg_mean_squared_error')
    print(f"Mean cross-validation score: {-scores.mean()}")
    #print("Best params:", grid_search.best_params_)

    plot_results(Y_test, y_pred, "SVR (tuned)")

