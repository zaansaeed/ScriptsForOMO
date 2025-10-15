import shap
import joblib
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# === 1. Load saved model and data ===

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
folder = os.path.join(os.getcwd())

# Adjust these filenames if different in your setup
model = joblib.load(folder+"/saved_model/random_forest_model.joblib")  # or model.joblib
model = model.named_steps['model']
X = pd.read_csv(folder + "/saved_model/X_copy.csv")                 # features
y = pd.read_csv(folder + "/saved_model/Y_copy.csv")                 # target

# === 2. Create SHAP explainer ===
explainer = shap.TreeExplainer(model)
shap_values = explainer.shap_values(X)

# Handle classifier (multiple classes) vs regressor
if isinstance(shap_values, list):
    shap_values = shap_values[1]  # use the positive class (e.g., "high permeability")
# === 3. Summary plot (global importance) ===
plt.title("Global Feature Importance (SHAP Summary Plot)")
shap.summary_plot(shap_values, X, show=False)
plt.tight_layout()
plt.savefig(folder + "/shap_summary_plot.png", dpi=300)
plt.close()

# === 4. Dependence plots (pairwise feature interactions) ===
# Pick your most important features automatically
important_features = X.columns[
    abs(shap_values).mean(0).argsort()[::-1][:5]
]

for f in important_features:
    shap.dependence_plot(f, shap_values, X, show=False)
    plt.tight_layout()
    plt.savefig(folder + f"/shap_dependence_{f}.png", dpi=300)
    plt.close()

# === 5. Optional: explain a single prediction ===
sample_idx = 0  # change to inspect another sample
shap.force_plot(
    explainer.expected_value,
    shap_values[sample_idx, :],
    X.iloc[sample_idx, :],
    matplotlib=True
)
plt.tight_layout()
plt.savefig(folder + "/shap_force_plot.png", dpi=300)
plt.close()

print("✅ SHAP analysis complete! Plots saved in:", folder)

