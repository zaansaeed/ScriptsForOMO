import os
import ast
import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors

from sklearn.model_selection import train_test_split, KFold, cross_val_score, RandomizedSearchCV
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_squared_error

import matplotlib.pyplot as plt

# -------------------------
# Utils
# -------------------------
def load_data(file_path: str) -> pd.DataFrame:
    return pd.read_csv(file_path)

def here() -> str:
    # Directory of this script
    return os.path.dirname(os.path.abspath(__file__))

def save_bar_chart(series: pd.Series, title: str, ylabel: str, outfile: str):
    plt.figure()
    series.plot(kind="bar")
    plt.title(title)
    plt.ylabel(ylabel)
    plt.xlabel("")
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    plt.close()

# -------------------------
# 0) Paths & chdir
# -------------------------
os.chdir(here())
MONOMER_CSV = os.path.join("/Users/zaansaeed/PycharmProjects/pythonProject/ScriptsForOMO/permeability/monomer_list.csv")
PEPTIDE_CSV = os.path.join("/Users/zaansaeed/PycharmProjects/pythonProject/ScriptsForOMO/permeability/processed_peptides.csv")
# If you want to save plots:
# PLOTS_DIR = os.path.join("/Users/zaan/PycharmProjects/ScriptsForOMO/permeability", "plots")
# os.makedirs(PLOTS_DIR, exist_ok=True)

# -------------------------
# 1) Load monomer list and build {symbol: (smiles, logP)}
# -------------------------
monomer_df = load_data(MONOMER_CSV)

monomers_list = {}
for symbol, smile in zip(monomer_df["Symbol"], monomer_df["replaced_SMILES"]):
    mol = Chem.MolFromSmiles(smile)
    if mol is None:
        continue
    logP = float(Descriptors.MolLogP(mol))
    monomers_list[symbol] = (smile, logP)

if not monomers_list:
    raise ValueError("monomers_list is empty. Check monomer_list.csv content.")

# -------------------------
# 2) Load peptides & compute per-position logP
# -------------------------
peptides = load_data(PEPTIDE_CSV)
if peptides["Sequence"].dtype == object:
    peptides["Sequence"] = peptides["Sequence"].apply(ast.literal_eval)

LOGP_DICTIONARY = {}
for seq, permeability, pid in zip(peptides["Sequence"], peptides["Permeability"], peptides["ID"]):
    if not isinstance(seq, (list, tuple)) or len(seq) != 6:
        continue

    logP_values = []
    bad = False
    for monomer in seq:
        if monomer not in monomers_list:
            bad = True
            break
        smile = monomers_list[monomer][0]
        mol = Chem.MolFromSmiles(smile)
        if mol is None:
            bad = True
            break
        logP_values.append(float(Descriptors.MolLogP(mol)))
    if bad:
        continue

    LOGP_DICTIONARY[tuple(seq)] = [logP_values, float(permeability), pid]

if not LOGP_DICTIONARY:
    raise ValueError("LOGP_DICTIONARY is empty. Check sequences/monomer symbols match your monomer list.")

# -------------------------
# 3) Build tidy feature table
# -------------------------
rows = []
for seq, (logP_values, permeability, pid) in LOGP_DICTIONARY.items():
    row = {"Permeability": permeability, "Sequence": "-".join(seq), "ID": pid}
    for i, v in enumerate(logP_values, start=1):
        row[f"Pos_{i}_logP"] = v
    rows.append(row)

df = pd.DataFrame(rows)

num_cols = ["Permeability"] + [f"Pos_{i}_logP" for i in range(1, 7)]
df[num_cols] = df[num_cols].apply(pd.to_numeric, errors="coerce")
df = df.dropna(subset=num_cols).reset_index(drop=True)

# ✅ Exclude invalid permeability rows (e.g., sentinel -10)
df = df[df["Permeability"] != -10].reset_index(drop=True)

# Drop constant columns (rare but safe)
const_cols = [c for c in num_cols if df[c].nunique(dropna=True) <= 1]
if const_cols:
    print("[Info] Dropping constant columns (no variance):", const_cols)
    df = df.drop(columns=const_cols)
    num_cols = [c for c in num_cols if c not in const_cols]

print(f"[Info] Final dataset shape: {df.shape}")
print(df.head())

# -------------------------
# 4) Correlations (optional print)
# -------------------------
feature_cols = [c for c in num_cols if c != "Permeability"]
pearson = df[["Permeability"] + feature_cols].corr(method="pearson")["Permeability"].drop("Permeability", errors="ignore")
spearman = df[["Permeability"] + feature_cols].corr(method="spearman")["Permeability"].drop("Permeability", errors="ignore")

print("\n=== Pearson r with Permeability ===")
print(pearson)
print("\n=== Spearman ρ with Permeability ===")
print(spearman)

# -------------------------
# 5) Split data
# -------------------------
X = df[feature_cols].values
y = df["Permeability"].values

# If dataset is small, use a slightly larger test split
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=max(0.2, min(0.3, 1.0 / max(2, len(df)))), random_state=42
)

# -------------------------
# 6) Hyperparameter tuning (RandomizedSearchCV)
#    Use neg_mean_squared_error so we can report RMSE directly
# -------------------------
from sklearn.pipeline import Pipeline

# Build a simple pipeline that only contains the RandomForestRegressor
rfr_pipeline = Pipeline([
    ("model", RandomForestRegressor(random_state=42, n_jobs=-1))
])

param_dist = {
    "model__n_estimators": np.linspace(200, 1000, num=9, dtype=int),
    "model__max_depth": [None] + list(np.arange(4, 26, 2)),
    "model__min_samples_split": [2, 4, 6, 8, 10, 20, 40],
    "model__min_samples_leaf": [1, 2, 3, 4, 5, 8, 10],
    "model__max_features": ["sqrt", "log2", 0.3, 0.5, 0.7],
    "model__bootstrap": [True, False],
    "model__max_samples": [None, 0.6, 0.8, 1.0],
    "model__criterion": ["squared_error", "friedman_mse"]
}

cv = KFold(n_splits=min(5, len(df)), shuffle=True, random_state=42)

rand_search = RandomizedSearchCV(
    estimator=rfr_pipeline,
    param_distributions=param_dist,
    n_iter=20,
    scoring="neg_mean_squared_error",
    cv=cv,
    random_state=42,
    n_jobs=-1,
    verbose=1
)


rand_search.fit(X_train, y_train)
best_rfrr = rand_search.best_estimator_
best_rfr = best_rfrr.named_steps["model"]
best_cv_rmse = np.sqrt(-rand_search.best_score_)
print("\n=== Hyperparameter Tuning Results ===")
print("Best params:", rand_search.best_params_)
print(f"Best CV RMSE (neg MSE scoring): {best_cv_rmse:.3f}")

# -------------------------
# 7) Evaluate tuned model + feature importances + chart
# -------------------------
y_pred = best_rfr.predict(X_test)
r2 = r2_score(y_test, y_pred)
rmse = np.sqrt(mean_squared_error(y_test, y_pred))

print("\n=== Tuned Random Forest (Test Set) ===")
print(f"R² (test):  {r2:.3f}")
print(f"RMSE (test): {rmse:.3f}")
importances = pd.Series(best_rfr.feature_importances_, index=feature_cols).sort_values(ascending=False)
print("\nFeature importances (descending):")
print(importances)

# -------------------------
# 8) Cross-Validation summary on full data (optional)
# -------------------------
cv_r2 = cross_val_score(best_rfr, X, y, cv=cv, scoring="r2", n_jobs=-1)
cv_rmse = np.sqrt(-cross_val_score(best_rfr, X, y, cv=cv, scoring="neg_mean_squared_error", n_jobs=-1))

print("\n=== Cross-Validation with Tuned RF ===")
print(f"R² (mean ± std):   {cv_r2.mean():.3f} ± {cv_r2.std():.3f}")
print(f"RMSE (mean ± std): {cv_rmse.mean():.3f} ± {cv_rmse.std():.3f}")

# -------------------------
# 9) Predicted vs True chart (1:1 line)
# -------------------------
plt.figure(figsize=(7, 7))
plt.scatter(y_test, y_pred, alpha=0.7, edgecolors='k', s=70)

# 1:1 perfect prediction line
lims = [min(y_test.min(), y_pred.min()), max(y_test.max(), y_pred.max())]
plt.plot(lims, lims, 'r--', lw=2, label='Ideal: y = x')

plt.xlabel("True Permeability", fontsize=13)
plt.ylabel("Predicted Permeability", fontsize=13)
plt.title(f"Predicted vs. True Permeability\nR² = {r2:.3f}, RMSE = {rmse:.3f}", fontsize=14)
plt.legend()
plt.grid(True, linestyle='--', alpha=0.5)
plt.tight_layout()
plt.show()

import joblib

# Create output folder if it doesn’t exist
output_dir = "saved_model"
os.makedirs(output_dir, exist_ok=True)

# Ensure data is in DataFrame form before saving
pd.DataFrame(X).to_csv(os.path.join(output_dir, "X_copy.csv"), index=False)
pd.DataFrame(y).to_csv(os.path.join(output_dir, "Y_copy.csv"), index=False, header=["Target"])

# Save trained Random Forest model
model_path = os.path.join(output_dir, "random_forest_model.joblib")
joblib.dump(best_rfrr, model_path)

print(f"\n=== Saved Outputs ===")
print(f"X_copy.csv  → {os.path.join(output_dir, 'X_copy.csv')}")
print(f"Y_copy.csv  → {os.path.join(output_dir, 'Y_copy.csv')}")
print(f"Model saved → {model_path}")