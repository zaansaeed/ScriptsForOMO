import streamlit as st
import yaml, subprocess, sys
from pathlib import Path

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
with st.expander("Rerun flags", expanded=False):
    for k, v in cfg["rerun"].items():
        cfg["rerun"][k] = st.checkbox(k, v)

## --- Paths ---------------------------------------------------------------- ##
dg = cfg["data_generation"]
with st.expander("Paths", expanded=True):
    dg["schrodinger_path"] = st.text_input(
        "Schrödinger install path", dg.get("schrodinger_path", "")
    )
    dg["main_dir"] = st.text_input("Main data directory", dg.get("main_dir", ""))

## --- ML ------------------------------------------------------------------- ##
ml = cfg["machine_learning"]
with st.expander("Machine-learning options", expanded=True):
    ml["model_name"] = st.text_input("Model name", ml.get("model_name", ""))
    ml["n_iter"] = st.number_input(
        "n_iter", min_value=1, max_value=10_000, value=ml.get("n_iter", 100)
    )
    ml["features_to_train_on"] = st.multiselect(
        "Features",
        ["side_chain_descriptors", "BWdistances", "molecular_descriptors"],
        default=ml.get("features_to_train_on", []),
    )

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
