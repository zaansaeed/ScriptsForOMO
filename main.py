import yaml
import functions as fs
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
    ##initialize the config file
    with open("config.yaml", "r") as f:
        config = yaml.safe_load(f)
        fs.init_config(config)
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

    smiles_lines = fs.load_lines(smiles_input_file)
    names_lines = fs.load_lines(names_input_file)
    targets_lines = fs.load_lines(targets_lines)

    main_dictionary = {}
    for i,name in enumerate(names_lines):
        main_dictionary[name] = (smiles_lines[i],targets_lines[i])
    main_dictionary = dict(natsorted(main_dictionary.items()))

    for name, (smile, target) in main_dictionary.items():
        if name != "pNP-43a" and name != "BICyP22":
            if not os.path.exists(main_dir+f"/Peptide_{name}"):
                os.mkdir(main_dir+f"/Peptide_{name}")
                data_logger.info(f"Created {main_dir+f'/Peptide_{name}'}")
            
            working_dir = main_dir+f"/Peptide_{name}"
            logging.info(f"Processing {name} in {working_dir}")
            fs.create_target_file(name,target,working_dir,config["machine_learning"]["target_name"])
            fs.smile_to_mae(smile, name,working_dir)
            fs.run_confSearch(working_dir,config["data_generation"]["wait_time_for_conf_search"])
            fs.mae_to_pdb(working_dir)
            fs.pdb_to_xyz(working_dir)
            fs.xyz_to_individual_xyz(name,working_dir)
            fs.extract_energies_to_csv(name,working_dir)
            ######
            fs.boltzmann_weight_distances(name,working_dir)
            fs.boltzmann_weight_dihedrals(name,working_dir)
            fs.create_new_descriptor("side_chain_descriptors",name,working_dir)
            ######
            data_logger.info(f"Finished processing {name}")
            data_logger.info("------------------------------------------------------------------------")

    main_dictionary = ML.filter_names(main_dictionary)
    names = [name for name in main_dictionary.keys()]

    features = config["machine_learning"]["features_to_train_on"]
    target_value = config["machine_learning"]["target_name"]

    X, Y = ML.create_model_data(names,features,main_dir,target_value)
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





