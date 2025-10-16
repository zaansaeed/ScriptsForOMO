# inverse_design.py
import os
from pathlib import Path
import joblib
import numpy as np
import pandas as pd
from scipy.optimize import differential_evolution

# === 0) Paths (same pattern as your shap_analysis.py) ===
os.chdir(os.path.dirname(os.path.abspath(__file__)))
folder = os.path.join(os.getcwd())

MODEL_PATH = folder + "/saved_model/random_forest_model.joblib"
X_PATH     = folder + "/saved_model/X_copy.csv"
Y_PATH     = folder + "/saved_model/Y_copy.csv"  # not strictly required, but available
MONOMER_LIST_CSV = folder + "/monomer_logP.csv"

# === 1) Load model and data ===
pipe = joblib.load(MODEL_PATH)
# If your saved object is a Pipeline with a 'model' step, extract it; otherwise use the object directly
model = getattr(pipe, "named_steps", {}).get("model", pipe)

X = pd.read_csv(X_PATH)
# y = pd.read_csv(Y_PATH)  # optional

FEATURES = list(X.columns)  # e.g., 6 columns for your 6 side chains





def target_region(
    target,
    constraints=None,           # {"Feature_2": 1.25} or {"Feature_3": (0.0, 2.0)}
    eps=0.05,                   # tolerance on |y_pred - target|
    n_samples=20000,            # how many candidates to try
    jitter_frac=0.5,           # add small noise to broaden coverage
    include_shap=False          # set True if you want SHAP overlay (needs TreeExplainer)
):
    #create random nmber generator
    rng = np.random.default_rng(0)
    #create base dataset, using sampling with replacement if n_samples > len(X)
    base = X.sample(n=min(len(X), n_samples), replace=len(X) < n_samples, random_state=0).copy()

    # small Gaussian jitter (per-feature) to explore nearby space
    std = X.std().replace(0, 1e-9).values #std of each feature, avoid zero std
    jitter = rng.normal(0, 1, size=(len(base), X.shape[1])) * (jitter_frac * std) # jitter taken from gaussian distruvtion and scaled by std and jitter_frac for realistic
    cand = base.values + jitter #add this jitter to base values, creating candidate set
    cand = pd.DataFrame(cand, columns=FEATURES)

    # enforce constraints
    if constraints:
        for raw_key, v in constraints.items():
            col = raw_key
            if col is None or col not in cand.columns:
                print(f"[target_region] Skipping constraint for '{raw_key}' — no matching feature in FEATURES={FEATURES}.")
                continue
            if isinstance(v, (tuple, list)):
                lo, hi = float(v[0]), float(v[1])
                cand[col] = cand[col].clip(lo, hi) #replace all candidates col values outside of (lo, hi) with lo or hi
            else:
                cand[col] = float(v) #replace all candidates col values with this value

    # keep candidates within robust bounds (1st–99th pct) to avoid crazy outliers
    ql, qh = X.quantile(0.01), X.quantile(0.99)
    for c in FEATURES:
        cand[c] = cand[c].clip(float(ql[c]), float(qh[c])) #clip each feature to be within 1st and 99th percentile of original feature values

    # predict and filter to those near the target
    try:
        yhat = pipe.predict(cand[FEATURES])
    except Exception:
        yhat = model.predict(cand[FEATURES])
    yhat = np.asarray(yhat).ravel()

    mask = np.abs(yhat - target) <= eps # only keep candidates within eps of target (boolean mask)
    hit = cand.loc[mask].copy() #only keep candidates within eps of target, selects only rows where mask is True
    hit["y_pred"] = yhat[mask]

    if hit.empty:
        return None, {"msg": f"No candidates within ±{eps} of target. Try increasing eps, n_samples, or jitter_frac."}

    monomers = pd.read_csv(MONOMER_LIST_CSV)
    # summarize feasible ranges
    def preview_list(lst, k=8):
        lst = list(lst)
        return lst if len(lst) <= k else lst[:k] + [f"...(+{len(lst)-k} more)"]

    summary = pd.DataFrame({
        "min": hit[FEATURES].min(),
        "q25": hit[FEATURES].quantile(0.25),
        "median": hit[FEATURES].median(),
        "q75": hit[FEATURES].quantile(0.75),
        "max": hit[FEATURES].max(),
        "possible_monomers" : [
        preview_list(
            monomers.loc[
                (monomers["logP"] >= hit[c].min()) &
                (monomers["logP"] <= hit[c].max()),
                "Symbol"
            ].dropna().drop_duplicates().tolist(),
            k=3
        )
        for c in FEATURES
    ]
    })

    out = {
        "n_candidates": int(len(cand)),
        "n_hits": int(len(hit)),
        "hit_rate": float(len(hit)/len(cand)),
        "target": float(target),
        "eps": float(eps)
    }

    # optional: SHAP on the hit set to explain which features “made” the target
    shap_summary = None
    if include_shap:
        try:
            import shap
            tree_model = getattr(pipe, "named_steps", {}).get("model", model)
            explainer = shap.TreeExplainer(tree_model)
            # If your model expects preprocessed inputs, ensure cand is what the model trained on
            phi = explainer.shap_values(hit[FEATURES].values)
            shap_summary = {
                "mean_abs_shap": dict(zip(FEATURES, np.mean(np.abs(phi), axis=0).tolist()))
            }
        except Exception as e:
            shap_summary = {"error": f"SHAP failed: {e}"}

    examples = hit.sample(min(10, len(hit)), random_state=1)
    examples = examples[FEATURES + ["y_pred"]]




    return {"feasible_ranges": summary, "examples": examples}, {**out, "shap": shap_summary}



# === 7) Example usage ===
if __name__ == "__main__":
    TARGET = -5
    CONSTRAINTS =  None #e.g., {"2": (0, 1), "5": 1.5}
    y = pd.read_csv(Y_PATH).squeeze()

    res, meta = target_region(TARGET, constraints=CONSTRAINTS, eps=0.06, n_samples=30000, include_shap=True)
    if res is None:
        print(meta["msg"])
    else:
        print("=== Feasible ranges (hits only) ===")
        print(res["feasible_ranges"].to_string())
        print("\n=== Example designs (first 10 sampled hits) ===")
        print(res["examples"].head(10).to_string(index=False))
        print("\nMETA:", meta)

    
    