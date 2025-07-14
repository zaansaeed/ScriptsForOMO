import streamlit as st, yaml, subprocess

cfg_path = 'config.yaml'
main_script = 'main.py'

cfg = yaml.safe_load(open(cfg_path))          # load
st.title('Pipeline Config')

# Rerun flags
for k in cfg['rerun']:
    cfg['rerun'][k] = st.checkbox(k, cfg['rerun'][k])

# Paths
dg = cfg['data_generation']
dg['schrodinger_path'] = st.text_input('schrodinger_path', dg['schrodinger_path'])
dg['main_dir']          = st.text_input('main_dir', dg['main_dir'])

# ML
ml = cfg['machine_learning']
ml['model_name'] = st.text_input('model', ml['model_name'])
ml['n_iter']     = st.number_input('n_iter', 1, 10_000, ml['n_iter'])
features = st.multiselect('features', ['side_chain_descriptors','BWdistances','molecular_descriptors'],
                          default=ml['features_to_train_on'])
ml['features_to_train_on'] = features

if st.button('Save & Run'):
    yaml.safe_dump(cfg, open(cfg_path, 'w'), sort_keys=False)
    subprocess.Popen(['python', main_script])
    st.success('Pipeline started!')