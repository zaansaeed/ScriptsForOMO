import os
import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, GradientBoostingRegressor
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import accuracy_score, classification_report, mean_squared_error


def pad_square_dataframe_to_array(df,target_size,fill_value):
    size = df.shape[0]
    original = df.to_numpy()
    padded = np.full((target_size,target_size),fill_value)
    padded[:size,:size] = original
    return pd.DataFrame(padded).values.tolist()

def create_X(main_dir): #takes in csv file and reads into array
    X= []
    for item in os.listdir(main_dir):
        working_dir = os.path.join(main_dir, item)
        if os.path.isdir(working_dir):
            os.chdir(working_dir)
            for file in os.listdir(working_dir):
                if file.endswith("-BWdistances.csv"):
                    data = pd.read_csv(file,header=None,index_col=None)
                    X.append(pad_square_dataframe_to_array(data,5,0))
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

#def run_RFR(X,Y):

