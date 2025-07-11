import yaml
import functions as fs
import logging
from main2 import *
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
    handler = logging.FileHandler(log_path, mode='w')
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)

    if not logger.handlers:
        logger.addHandler(handler)

    return logger

def main():

    with open("config.yaml", "r") as f:
        config = yaml.safe_load(f)
        fs.init_config(config)
    main_dir = os.path.abspath(config["data_generation"]["main_dir"])

    data_logger = setup_logging_data(main_dir)
    ml_logger = setup_logging_ml(main_dir)
    data_logger.info("Starting data generation...")


    schrodinger_path = config["data_generation"]["schrodinger_path"]
    smiles_input_file = config["data_generation"]["smiles_input_file"]
    names_input_file = config["data_generation"]["names_input_file"]
    percents_input_file = config["data_generation"]["target_input_file"]

    os.chdir(main_dir)

    smiles_lines = fs.load_lines(smiles_input_file)
    names_lines = fs.load_lines(names_input_file)
    percents_lines = fs.load_lines(percents_input_file)

    main_dictionary = {}
    for i,name in enumerate(names_lines):
        main_dictionary[name] = (smiles_lines[i],percents_lines[i])

    main_dictionary = dict(natsorted(main_dictionary.items()))

    for name, (smile, percent) in main_dictionary.items():
        if name != "pNP-43a" and name != "BICyP22":
            if not os.path.exists(main_dir+f"/Peptide_{name}"):
                os.mkdir(main_dir+f"/Peptide_{name}")
                data_logger.info(f"Created {main_dir+f'/Peptide_{name}'}")
            working_dir = main_dir+f"/Peptide_{name}"
            logging.info(f"Processing {name} in {working_dir}")
            fs.create_target_file(name,percent,working_dir)
            fs.smile_to_mae(smile, name,working_dir)
            fs.run_confSearch(working_dir,120)
            fs.mae_to_pdb(working_dir)
            fs.pdb_to_xyz(working_dir)
            fs.xyz_to_individual_xyz(name,working_dir)
            fs.extract_energies_to_csv(name,working_dir)
            ######
            fs.boltzmann_weight_distances(name,working_dir)
            fs.extract_boltzmann_weighted_dihedrals_normalized(name,working_dir)
            #fs.create_new_descriptor("side_chain_descriptors",name,working_dir)
            data_logger.info(f"Finished processing {name}")
            data_logger.info("------------------------------------")

    main_dictionary = ML.filter_names(main_dictionary)

    names = [name for name in main_dictionary.keys()]
    #BWdistances, side_chain_descriptors, BWDihedralNormalized, molecular_descriptors
    features = ["side_chain_descriptors"]
    target_value = "target"

    X, Y = ML.create_model_data(names,features,main_dir,target_value)
    print(X.shape)
    print(Y.shape)


    ML.run_elasticnet(X,Y,5,.2)


if __name__ == "__main__":
    main()





