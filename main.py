import os
import functions_for_smile_to_xyz as fs
import ML_functions as ML

schrodinger_path ="/opt/schrodinger/suites2024-3/"
# Define the working directory (where results will be stored)
#/Users/zaan/zasaeed@g.hmc.edu - Google Drive/My Drive/OMO Lab - Peptide Cyclization - Zaan Saeed/Data/NewPeptideLibrary
main_dir = os.path.abspath("/Users/zaansaeed/Peptides")
smiles_input_file = "all_peptides.smi"
names_input_file = "all_names.txt"


os.chdir(main_dir)
with open(smiles_input_file, "r") as f:
    smiles_lines = f.readlines()
    smiles_lines = [line.strip() for line in smiles_lines]
with open(names_input_file, "r") as f:
    names_lines = f.readlines()
    names_lines =[name.strip() for name in names_lines]


update_matrices = 0
for i, name in enumerate(names_lines): #processing : smiles -> xyzs
    working_dir =main_dir+f"/Peptide_{name}"
    print(name)

    if update_matrices == 1:
        if working_dir in os.listdir(main_dir):
            fs.boltzmann_weight_energies(name, working_dir, update_matrices)
    else:
        fs.smile_to_mae(smiles_lines[i], name)
        fs.run_confsearch(name,working_dir)
        fs.mae_to_pdb(name,working_dir)
        fs.pdb_to_xyz(name,working_dir)
        fs.xyz_to_individual_xyz(name,working_dir)
        fs.extract_energies_to_csv(name,working_dir)
        fs.boltzmann_weight_energies(name,working_dir,update_matrices)


#outputs = ML.six_over_target_percents(ML.create_outputs(main_dir))

#X = ML.create_X(main_dir) #ready for input
#Y = ML.create_Y(outputs,.75) #ready for input

#testing = X[0]


#ML.run_RFR(X,outputs)
#ML.run_SVR(X,outputs)
#ML.run_LR(X,outputs)
