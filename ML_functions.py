import os
import matplotlib.pyplot as  plt
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, GradientBoostingRegressor, HistGradientBoostingRegressor
from sklearn.model_selection import train_test_split, GridSearchCV, RandomizedSearchCV ,cross_val_score, KFold
from sklearn.metrics import accuracy_score, classification_report, mean_squared_error, r2_score, make_scorer, \
     f1_score,mean_absolute_error, r2_score, root_mean_squared_error
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.svm import SVR, SVC
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline


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
        if item >=cutoff:
            Y.append(1)
        else:
            Y.append(0)
    return Y




def run_RFC(X,Y):

    X_train, X_test, y_train, y_test = train_test_split(X, Y, test_size=0.3, random_state=42)

    pipeline = Pipeline([

        ('pca', PCA(n_components=10)),
        ('rf', RandomForestClassifier(random_state=42))
    ])

    pipeline.fit(X_train, y_train)

    # Predict on test data
    y_pred = pipeline.predict(X_test)



    # Evaluate
    print("Test accuracy:", accuracy_score(y_test, y_pred))
    accuracy = accuracy_score(y_test, y_pred)
    print("Accuracy: ", accuracy)
    print("F1: ", f1_score(y_test, y_pred, average='macro'))
    print(classification_report(y_test, y_pred))
    plot_results(y_test, y_pred, pipeline)

    from sklearn.dummy import DummyClassifier

    dummy = DummyClassifier(strategy="most_frequent")
    dummy.fit(X_train, y_train)
    print("Baseline accuracy:", dummy.score(X_test, y_test))
    print("Baseline f1:",f1_score(y_test, dummy.predict(X_test), average='macro'))

def run_SVM(X,Y):
    scaler = StandardScaler()
    X = scaler.fit_transform(X)
    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.3, random_state=42)
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
    grid_search.fit(X_train, Y_train)
    y_pred = grid_search.predict(X_test)



    accuracy = accuracy_score(Y_test, y_pred)
    print("Accuracy: ", accuracy)
    print("F1: ", f1_score(Y_test, y_pred))
    print(classification_report(Y_test, y_pred))
    print(grid_search.best_params_)
    plot_results(Y_test, y_pred, grid_search)

    from sklearn.dummy import DummyClassifier

    dummy = DummyClassifier(strategy="most_frequent")
    dummy.fit(X_train, Y_train)
    print("Baseline accuracy:", dummy.score(X_test, Y_test))
    print("Baseline f1:", f1_score(Y_test, dummy.predict(X_test), average='macro'))


def run_RFR(X,Y):

    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.4, random_state=42)

    param_grid = {
        'n_estimators': [150,1000],  # Fewer trees to keep it lightweight
        'max_depth': [10,40, None],  # Shallow trees to prevent overfitting
        'min_samples_split': [2, 10],  # Try larger splits to regularize
        'min_samples_leaf': [2, 4,5],  # Larger leaves reduce model complexity
        'max_features': ['sqrt'],  # Limit number of features at each split
        'bootstrap': [True]  # Usually better for small data
    }
    rf = RandomForestRegressor(random_state=42)
    grid_search = GridSearchCV(
        estimator=rf,
        param_grid=param_grid,
        scoring='neg_mean_squared_error',
        cv=6,
        n_jobs=-1,
        verbose=1
    )


    rf.fit(X_train,Y_train)
    y_pred = rf.predict(X_test)
    mse = mean_squared_error(Y_test, y_pred)

    print("Results on testing data:,", y_pred)
    print("True percents", Y_test)
   # print("best params:", grid_search.best_params_)
    scores = cross_val_score(rf, X_train, Y_train, cv=6, scoring='neg_mean_squared_error')
    print(f"Mean cross-validation score: {-scores.mean()}")
    print("MSE: ", mse)
    print("MAE,", mean_absolute_error(Y_test, y_pred))
    plot_results(Y_test, y_pred,"random forest")

    baseline = np.full_like(Y_test, np.mean(Y_train))
    baseline_mae = mean_absolute_error(Y_test, baseline)
    print("Baseline MAE (predicting mean):", baseline_mae)

def run_SVR(X,Y):


    X_train, X_test, Y_train, Y_test = train_test_split(X, Y, test_size=0.3, random_state=42)

    ''' pipeline = Pipeline([
        ('scaler', StandardScaler()),
        #('pca', PCA(n_components=15)),  # Let GridSearch decide n_components
        ('svr', SVR())
    ])

    param_grid = {
        'svr__kernel': ['rbf'],
        'svr__C': [1000],
        'svr__gamma': [ .001],
    }'''



    svr = SVR()
    #grid_search = GridSearchCV(pipeline, param_grid, scoring='r2', cv=5, n_jobs=-1)
    svr.fit(X_train, Y_train)
    #best_model = grid_search.best_estimator_

    y_pred = svr.predict(X_test)
    y_pred = np.clip(y_pred, 0, 1)


    scores = cross_val_score(svr, X_train, Y_train, cv=5, scoring='r2')
    print(f"Mean cross-validation score: {-scores.mean()}")
    #print("Best params:", grid_search.best_params_)
    print("mse: ", mean_squared_error(Y_test, y_pred))
    print("MAE,", mean_absolute_error(Y_test, y_pred))

    print("R2: ", r2_score(Y_test, y_pred))
    print("RMSE:", np.sqrt(mean_squared_error(Y_test, y_pred)))

    plot_results(Y_test, y_pred, "SVR (tuned)")

