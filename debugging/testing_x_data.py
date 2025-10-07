import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split, RandomizedSearchCV, KFold, cross_val_score
from sklearn.metrics import r2_score, mean_squared_error
import os
# -------------------------------------------------
# 1. Load data
# -------------------------------------------------
os.chdir("/Users/zaan/PycharmProjects/ScriptsForOMO/debugging")
X = pd.read_csv("X.csv", header=None)
y = pd.read_csv("y.csv", header=None).squeeze()  # convert to 1D array

# -------------------------------------------------
# 2. Extract every 6th column of each 16-column block
# -------------------------------------------------
selected_cols = [i + 5 for i in range(0, X.shape[1], 16) if i + 5 < X.shape[1]]
X_sel = X.iloc[:, selected_cols]

#print each row of x_sel

# Give columns readable names
X_sel.columns = [f"Block{(i//16)+1}_Col6" for i in range(0, X.shape[1], 16) if i + 5 < X.shape[1]]

print(f"[Info] Selected {X_sel.shape[1]} features from {X.shape[1]} total columns.")
print(X_sel.head())

# -------------------------------------------------
# 3. Split into training and testing
# -------------------------------------------------
X_train, X_test, y_train, y_test = train_test_split(X_sel, y, test_size=0.2, random_state=42)

# -------------------------------------------------
# 4. Hyperparameter tuning with RandomizedSearchCV
# -------------------------------------------------
param_dist = {
    "n_estimators": np.arange(100, 801, 50),
    "max_depth": [None] + list(np.arange(3, 21)),
    "min_samples_split": [2, 5, 10],
    "min_samples_leaf": [1, 2, 4],
    "max_features": ["auto", "sqrt", "log2", 0.5, 0.75, 1.0]
}

rfr = RandomForestRegressor(random_state=42, n_jobs=-1)
cv = KFold(n_splits=min(5, len(X_sel)), shuffle=True, random_state=42)

rand_search = RandomizedSearchCV(
    estimator=rfr,
    param_distributions=param_dist,
    n_iter=50,
    scoring="neg_mean_squared_error",
    cv=cv,
    random_state=42,
    n_jobs=-1,
    verbose=1
)
rand_search.fit(X_train, y_train)

best_rfr = rand_search.best_estimator_

print("\n=== Best Hyperparameters ===")
print(rand_search.best_params_)
print(f"Best CV RMSE: {np.sqrt(-rand_search.best_score_):.3f}")

# -------------------------------------------------
# 5. Evaluate on the test set
# -------------------------------------------------
y_pred = best_rfr.predict(X_test)

r2 = r2_score(y_test, y_pred)
rmse = np.sqrt(mean_squared_error(y_test, y_pred))

print("\n=== Test Set Performance ===")
print(f"R²: {r2:.3f}")
print(f"RMSE: {rmse:.3f}")

# -------------------------------------------------
# 6. Feature importances
# -------------------------------------------------
importances = pd.Series(best_rfr.feature_importances_, index=X_sel.columns).sort_values(ascending=False)

plt.figure(figsize=(10, 6))
sns.barplot(x=importances.values, y=importances.index, palette="crest")
plt.title("Feature Importance (Random Forest)")
plt.xlabel("Importance")
plt.ylabel("Feature")
plt.tight_layout()
plt.show()

# -------------------------------------------------
# 7. Cross-validation performance summary
# -------------------------------------------------
cv_r2 = cross_val_score(best_rfr, X_sel, y, cv=cv, scoring="r2", n_jobs=-1)
cv_rmse = np.sqrt(-cross_val_score(best_rfr, X_sel, y, cv=cv, scoring="neg_mean_squared_error", n_jobs=-1))

print("\n=== Cross-Validation Results ===")
print(f"R² (mean ± std):   {cv_r2.mean():.3f} ± {cv_r2.std():.3f}")
print(f"RMSE (mean ± std): {cv_rmse.mean():.3f} ± {cv_rmse.std():.3f}")


# 1. Pairwise partial dependence (2-way interactions)
from sklearn.inspection import partial_dependence, PartialDependenceDisplay

PartialDependenceDisplay.from_estimator(best_rfr, X_sel, ['Block6_Col6'])
PartialDependenceDisplay.from_estimator(best_rfr, X_sel, [('Block2_Col6','Block6_Col6')])
plt.show()

# 2. SHAP values (per-sample feature impact)
import shap
explainer = shap.TreeExplainer(best_rfr)
shap_values = explainer(X_sel)
shap.summary_plot(shap_values, X_sel, plot_type="bar")

