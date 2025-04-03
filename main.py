import os
import functions_for_smile_to_xyz as fs


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



for i, name in enumerate(names_lines):
    working_dir =main_dir+f"/Peptide_{name}"
    print(working_dir)
    print(name)

    fs.smile_to_mae(smiles_lines[i], name)
    fs.run_confsearch(name,working_dir)
    fs.mae_to_pdb(name,working_dir)
    fs.pdb_to_xyz(name,working_dir)
    fs.xyz_to_individual_xyz(name,working_dir)
    fs.extract_energies_to_csv(name,working_dir)
    fs.boltzmann_weight_energies(name,working_dir)
