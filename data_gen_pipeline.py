import functions as fs
from main import *


def main():
    schrodinger_path ="/opt/schrodinger/suites2024-3/"
    main_dir = os.path.abspath("/Users/zaansaeed/Desktop/NewPeptides")
    smiles_input_file = "all_peptides.smi"
    names_input_file = "all_names.txt"
    # Define the working directory (where results will be stored)
    #/Users/zaan/zasaeed@g.hmc.edu - Google Drive/My Drive/OMO Lab - Peptide Cyclization - Zaan Saeed/Data/NewPeptideLibrary

    os.chdir(main_dir)
    with open(smiles_input_file, "r") as f:
        smiles_lines = f.readlines()
        smiles_lines = [line.strip() for line in smiles_lines]
    with open(names_input_file, "r") as f:
        names_lines = f.readlines()
        names_lines = [name.strip() for name in names_lines]


    update_matrices = False


    for i, name in enumerate(names_lines): #processing : smiles -> xyzs
        if not os.path.exists(main_dir+f"/Peptide_{name}"):
            os.mkdir(main_dir+f"/Peptide_{name}")
        working_dir = main_dir+f"/Peptide_{name}"
        print(name)
        fs.smile_to_mae(smiles_lines[i], name,working_dir)
        fs.run_confSearch(name,working_dir)
        fs.mae_to_pdb(name,working_dir)
        fs.pdb_to_xyz(name,working_dir)
        fs.xyz_to_individual_xyz(name,working_dir)
        fs.extract_energies_to_csv(name,working_dir)
        fs.boltzmann_weight_energies(name,working_dir,update_matrices)


    fs.extract_boltzmann_weighted_dihedrals_normalized(main_dir)


if __name__ == "__main__":
    main()





