import os
import functions_for_smile_to_xyz as fs
import ML_functions as ML

schrodinger_path ="/opt/schrodinger/suites2024-3/"
# Define the working directory (where results will be stored)
#/Users/zaan/zasaeed@g.hmc.edu - Google Drive/My Drive/OMO Lab - Peptide Cyclization - Zaan Saeed/Data/NewPeptideLibrary
main_dir = os.path.abspath("/Users/zaan/zasaeed@g.hmc.edu - Google Drive/Shared drives/OMO Lab/Projects/OMO Lab - Zaan Saeed/Data/Peptides")
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
fs.extract_boltzmann_weighted_dihedrals_normalized()








#feature: BWdistances or BWdihedrals or BWDihedralNormalized

X = ML.create_X(main_dir, "BWDihedralNormalized") #ready for input
import pandas as pd
df = pd.DataFrame(X)
#df.to_csv("/Users/zaan/zasaeed@g.hmc.edu - Google Drive/Shared drives/OMO Lab/Projects/OMO Lab - Zaan Saeed/Data/Peptides/X1.csv", header=False, index=False)

Y = ML.six_over_target_percents(ML.create_outputs(main_dir))
#Y = create_Y(Y,.75 )


#print("Count of 0s:", Y.count(0))
#print("Count of 1s:", Y.count(1))


#ML.run_RFC(X,Y)
#ML.run_SVM(X,Y)

ML.run_RFR(X,Y)
#ML.run_SVR(X,Y)
