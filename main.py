import yaml
import functions as funcs
import logging
import os
from natsort import natsorted
import ML_functions as ML


def setup_logging_data(directory_of_log):
    logger = logging.getLogger("data_logger")
    logger.setLevel(logging.INFO)
    logger.propagate = False  # Prevents double logging if root logger is also configured

    log_path = os.path.join(directory_of_log, "data_processing.log")
    handler = logging.FileHandler(log_path, mode='w')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)

    # Avoid adding the handler twice if this function is called again
    if not logger.handlers:
        logger.addHandler(handler)

    return logger

def setup_logging_ml(directory_of_log):
    logger = logging.getLogger("ml")
    logger.setLevel(logging.INFO)
    logger.propagate = False

    log_path = os.path.join(directory_of_log, "ml.log")
    handler = logging.FileHandler(log_path, mode='a')
    formatter = logging.Formatter("%(asctime)s,%(levelname)s,%(message)s")
    handler.setFormatter(formatter)

    if not logger.handlers:
        logger.addHandler(handler)

    return logger

def main():
    print("Starting the pipeline...")
    with open("config.yaml", "r") as f:
        config = yaml.safe_load(f)
        funcs.init_config(config)
        ML.init_config(config)
    
    main_dir = os.path.abspath(config["data_generation"]["main_dir"])

    data_logger = setup_logging_data(main_dir)
    ml_logger = setup_logging_ml(main_dir)
    data_logger.info("Starting data generation...")


    schrodinger_path = config["data_generation"]["schrodinger_path"]
    smiles_input_file = config["data_generation"]["smiles_input_file"]
    names_input_file = config["data_generation"]["names_input_file"]
    targets_lines = config["data_generation"]["target_input_file"]

    os.chdir(main_dir)

    smiles_lines = funcs.load_lines(smiles_input_file)
    names_lines = funcs.load_lines(names_input_file)
    targets_lines = funcs.load_lines(targets_lines)

    main_dictionary = {}
    for i,name in enumerate(names_lines):
        main_dictionary[name] = (smiles_lines[i],targets_lines[i])
    main_dictionary = dict(natsorted(main_dictionary.items()))

    for name, (smile, target) in main_dictionary.items():
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smile)
        if "*" not in smile and funcs.find_n_terminus(mol) is not None and len(funcs.add_amides(mol)) >= config["data_generation"]["number_of_residues"]-1: ###temporary line to skip some peptides, skip ones that dont have n_terminus
            if not os.path.exists(main_dir+f"/Peptide_{name}"):
                os.mkdir(main_dir+f"/Peptide_{name}")
                data_logger.info(f"Created {main_dir+f'/Peptide_{name}'}")
            print(name)
            working_dir = main_dir+f"/Peptide_{name}"
            logging.info(f"Processing {name} in {working_dir}")
            funcs.create_target_file(name,target,working_dir,config["machine_learning"]["target_name"])
            funcs.smile_to_mae(smile, name,working_dir)
            funcs.run_confSearch(working_dir,config["data_generation"]["conf_search_settings"]["wait_time_for_conf_search"])
            funcs.mae_to_pdb(working_dir)
            funcs.pdb_to_xyz(working_dir)
            funcs.xyz_to_individual_xyz(name,working_dir)
            funcs.extract_energies_to_csv(name,working_dir)
            data_logger.info(f"Finished processing {name}")
            data_logger.info("------------------------------------------------------------------------")

            ########## feature genereation
            funcs.boltzmann_weight_distances(name,working_dir)
            funcs.boltzmann_weight_dihedrals(name,working_dir)
            funcs.create_new_descriptor("side_chain_descriptors",name,working_dir)
            ##########
        else:
            file_path = f"{main_dir}/skipped_peptides.txt"
            with open(file_path, "a+") as outfile:
                outfile.seek(0)  # Move to start to read contents
                contents = [line.strip() for line in outfile.readlines()]
                if name not in contents:
                    outfile.write(f"{name}\n")  # Write without leading newline

    main_dictionary = ML.filter_names(main_dictionary)
    names = [name for name in main_dictionary.keys()]

    features = config["machine_learning"]["features_to_train_on"]
    target_value = config["machine_learning"]["target_name"]

    X, Y = ML.create_model_data(names,features,main_dir,target_value)
    X = X.tolist()
    Y = Y.tolist()

    for i in reversed(range(len(X))):
        if Y[i] == -10.0:
            del X[i]
            del Y[i]
    import numpy as np
    X = np.array(X)
    Y = np.array(Y)
    ml_logger.info(f"X shape: {X.shape}, Y shape: {Y.shape}")
    ml_logger.info(f"Features: {features}")
    ml_logger.info(f"Target value: {target_value}")

    model_map = {
        "RandomForestRegressor": ML.run_RFR,
        "ElasticNet": ML.run_ElasticNet,
        "SVR": ML.run_SVR
    }
    
    if config["machine_learning"]["model_name"] in model_map:
        ml_logger.info(f"Starting {config["machine_learning"]["model_name"]} training...")
        model_function = model_map[config["machine_learning"]["model_name"]]
        model_function(X, Y)
if __name__ == "__main__":
    main()





