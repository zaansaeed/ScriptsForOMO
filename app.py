import streamlit as st
import yaml, subprocess, sys
from pathlib import Path
import matplotlib.pyplot as plt


# --------------------------------------------------------------------------- #
#  PATHS
# --------------------------------------------------------------------------- #
BASE_DIR    = Path(__file__).resolve().parent
cfg_path    = BASE_DIR / "config.yaml"
main_script = BASE_DIR / "main.py"           # still a Path object

# --------------------------------------------------------------------------- #
#  LOAD CONFIG (fallback to empty dicts instead of crashing)
# --------------------------------------------------------------------------- #
try:
    with cfg_path.open() as f:
        cfg = yaml.safe_load(f) or {}
except FileNotFoundError:
    cfg = {}

cfg.setdefault("rerun", {})
cfg.setdefault("data_generation", {})
cfg.setdefault("machine_learning", {})

# --------------------------------------------------------------------------- #
#  UI
# --------------------------------------------------------------------------- #
st.title("Pipeline Config")

## --- Rerun flags ---------------------------------------------------------- ##
with st.expander("Recalculation flags", expanded=True):
    for k, v in cfg["rerun"].items():
        cfg["rerun"][k] = st.checkbox(k, v)

## --- Paths ---------------------------------------------------------------- ##
dg = cfg["data_generation"]
with st.expander("Paths", expanded=True):
    dg["schrodinger_path"] = st.text_input(
        "Schrödinger install path", dg.get("schrodinger_path", "")
    )
    dg["main_dir"] = st.text_input("Main data directory", dg.get("main_dir", ""))
    dg["smiles_input_file"] = st.text_input("SMILES input file name", dg.get("smiles_input_file", ""))
    dg["names_input_file"] = st.text_input("Names input file name", dg.get("names_input_file", ""))
    dg["target_input_file"] = st.text_input("Target input file name", dg.get("target_input_file", ""))


## --- confsearch setup --------------------------------------------- ##
cf = cfg["data_generation"]["conf_search_settings"]
with st.expander("ConfSearch settings", expanded=False):

    # Ensure the dict exists even if the YAML was missing it
    cfg.setdefault("data_generation", {}).setdefault("conf_search_settings", {})
    cf = cfg["data_generation"]["conf_search_settings"]

    # 0️⃣  Wait time for ConfSearch (seconds) --------------------------------
    cf["wait_time_for_conf_search"] = st.number_input(
        "Wait time for ConfSearch (seconds)",
        min_value=0,
        max_value=3600,
        step=1,
        value=int(cf.get("wait_time_for_conf_search", 120)),
    )
    # 1️⃣  Method ------------------------------------------------------------

    cf["CONFSEARCH_METHOD"] = st.selectbox(
        "Method",
        options=["MCMM", "LBFGS", "TNC"],          # add others if you support them
        index=["MCMM", "LBFGS", "TNC"].index(cf.get("CONFSEARCH_METHOD", "MCMM")),
    )

    # 2️⃣  Force field -------------------------------------------------------
    cf["FORCE_FIELD"] = st.selectbox(
        "Force field",
        options=["OPLS_2005", "OPLS4", "AMBER99SB"],
        index=["OPLS_2005", "OPLS4", "AMBER99SB"].index(
            cf.get("FORCE_FIELD", "OPLS_2005")
        ),
    )

    # 3️⃣  Solvent -----------------------------------------------------------
    cf["SOLVENT"] = st.selectbox(
        "Solvent",
        options=["None", "Water", "DMSO", "Ethanol"],
        index=["None", "Water", "DMSO", "Ethanol"].index(cf.get("SOLVENT", "None")),
    )

    # 4️⃣  Dielectric constant (float) --------------------------------------
    cf["DIELECTRIC_CONSTANT"] = st.number_input(
        "Dielectric constant",
        min_value=1.0,
        max_value=80.0,
        step=0.1,
        value=float(cf.get("DIELECTRIC_CONSTANT", 1.0)),
    )

    # 5️⃣  Charges from ------------------------------------------------------
    cf["CHARGES_FROM"] = st.selectbox(
        "Charges from",
        options=["Force field", "LigPrep", "QM"],
        index=["Force field", "LigPrep", "QM"].index(cf.get("CHARGES_FROM", "Force field")),
    )

    # 6️⃣  Cut-off (Å) ------------------------------------------------------- ###update to handle float/values
    cf["CUTOFF"] = st.text_input(
        "Cut-off (Å) " ,value=cf.get("CUTOFF",""),
    )

    # 7️⃣  Minimisation method ----------------------------------------------
    cf["MINI_METHOD"] = st.selectbox(
        "Minimisation method",
        options=["PRCG", "SD", "LBFGS"],
        index=["PRCG", "SD", "LBFGS"].index(cf.get("MINI_METHOD", "PRCG")),
    )

    # 8️⃣  Maximum iterations (int) -----------------------------------------
    cf["MAXIMUM_ITERATIONS"] = st.number_input(
        "Maximum iterations",
        min_value=100,
        max_value=100_000,
        step=100,
        value=int(cf.get("MAXIMUM_ITERATIONS", 2500)),
    )

    # 9️⃣  Converge on -------------------------------------------------------
    cf["CONVERGE_ON"] = st.selectbox(
        "Converge on",
        options=["Gradient", "Energy"],
        index=["Gradient", "Energy"].index(cf.get("CONVERGE_ON", "Gradient")),
    )

    # 10️⃣  Convergence threshold (kcal/mol/Å) ------------------------------
    cf["CONVERGENCE_THRESHOLD"] = st.number_input(
        "Convergence threshold",
        min_value=0.001,
        max_value=10.0,
        step=0.001,
        format="%.3f",
        value=float(cf.get("CONVERGENCE_THRESHOLD", 0.05)),
    )

    # 11️⃣  Outconfs per search (int) ---------------------------------------
    cf["OUTCONFS_PER_SEARCH"] = st.number_input(
        "Out-confs per search",
        min_value=10,
        max_value=50_000,
        step=10,
        value=int(cf.get("OUTCONFS_PER_SEARCH", 10_000)),
    )

    # 12️⃣  ConfSearch steps (int) ------------------------------------------
    cf["CONFSEARCH_STEPS"] = st.number_input(
        "ConfSearch steps",
        min_value=100,
        max_value=1_000_000,
        step=100,
        value=int(cf.get("CONFSEARCH_STEPS", 1000)),
    )

    # 13️⃣  Steps per rotatable (int) ---------------------------------------
    cf["CONFSEARCH_STEPS_PER_ROTATABLE"] = st.number_input(
        "Steps per rotatable bond",
        min_value=10,
        max_value=10_000,
        step=10,
        value=int(cf.get("CONFSEARCH_STEPS_PER_ROTATABLE", 100)),
    )

    # 14️⃣  Energy window (kcal/mol) ----------------------------------------
    cf["ENERGY_WINDOW"] = st.number_input(
        "Energy window (kcal/mol)",
        min_value=0.0,
        max_value=200.0,
        step=0.1,
        value=float(cf.get("ENERGY_WINDOW", 104.6)),
    )

    # 15️⃣  Torsion sampling level ------------------------------------------
    cf["CONFSEARCH_TORSION_SAMPLING"] = st.selectbox(
        "Torsion sampling",
        options=["Coarse", "Intermediate", "Fine"],
        index=["Coarse", "Intermediate", "Fine"].index(
            cf.get("CONFSEARCH_TORSION_SAMPLING", "Intermediate")
        ),
    )

## --- ML ------------------------------------------------------------------- ##
ml = cfg["machine_learning"]
with st.expander("Machine-learning options", expanded=True):
    models = ml["models"]
    scorers = ["neg_mean_squared_error", "r2", "neg_mean_absolute_error","neg_root_mean_squared_error","custom_scorer"]
    ml["model_name"] = st.selectbox(label ="Model", options=models, index=models.index(ml.get("model_name", "ElasticNet")))

    ml["n_iter"] = st.number_input(
        "n_iter", min_value=1, max_value=10_000, value=ml.get("n_iter", 100)
    )
    ml["features_to_train_on"] = st.multiselect(
        "Features",
        ["side_chain_descriptors", "BWdistances", "molecular_descriptors","BWDihedralNormalized"],
        default=ml.get("features_to_train_on", []),
    )
    ml["save_model"] = st.checkbox(
        "Save model", value=ml.get("save_model", False) )

# --------------------------------------------------------------------------- #
#  SAVE + RUN
# --------------------------------------------------------------------------- #
if st.button("Save & Run"):

    # 1) write the YAML
    with cfg_path.open("w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)

    # 2) start the pipeline with the SAME interpreter that runs Streamlit
    proc = subprocess.Popen(
        [sys.executable, str(main_script)],  # sys.executable is the venv’s python
        cwd=BASE_DIR,                       # so main.py sees the repo root
    )

    st.success(f"Pipeline started! (PID {proc.pid})")
    st.write(f"Config saved to **{cfg_path.relative_to(BASE_DIR)}**")
