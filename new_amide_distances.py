from array import array
from tkinter.font import names

from rdkit import Chem
from ML_functions import *
from natsort import natsorted
from collections import defaultdict

from functions import boltzmann
from functions import addAmides

def distance_between_two_atoms(mol,ID1,ID2):
    conf = mol.GetConformer()
    atom_1_pos = conf.GetAtomPosition(ID1)
    atom_2_pos = conf.GetAtomPosition(ID2)
    return np.linalg.norm(atom_1_pos - atom_2_pos)

def getAmideDistances(amideGroups,mol):
    distance_matrix = []
    for i,amid1 in enumerate(amideGroups):
        distances_for_one_amid = []
        for j,amid2 in enumerate(amideGroups):
            if i == j:
                distances_for_one_amid.append([0,0,0])
            else:
                amid1C=amid1.getC()
                amid2C=amid2.getC()
                amid1N=amid1.getN()
                amid2N=amid2.getN()
                amid1O=amid1.getO()
                amid2O=amid2.getO()
                distances_for_one_amid.append([distance_between_two_atoms(mol,amid1C,amid2C),distance_between_two_atoms(mol,amid1O,amid2O),distance_between_two_atoms(mol,amid1N,amid2N)])
        distance_matrix.append(distances_for_one_amid)
    return distance_matrix

def load_xyz_coords(mol, xyz_path):
    conf = Chem.Conformer(mol.GetNumAtoms())

    with open(xyz_path, 'r') as f:
        lines = f.readlines()
        lines = lines[2:]

        for i, line in enumerate(lines):
            tokens = line.strip().split()
            if len(tokens) < 4:
                continue
            x, y, z = map(float, tokens[1:4])
            conf.SetAtomPosition(i, (x, y, z))

    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol


def create_new_descriptor(descriptor_name,directory_of_peptides):
    main_dir = os.path.abspath(directory_of_peptides)
    os.chdir(main_dir)
    for folder in natsorted(os.listdir(main_dir)):
        if os.path.isdir(folder):
            os.chdir(folder) #currenlty working in Peptide_{name}
            if not os.path.exists(f"{descriptor_name}.csv"):
                working_dir = os.getcwd()
                name = folder.split("_")[1]
                smiles_string = open(f"{name}.smi").read().strip()
                print(name)
                peptide_descriptors = []
                mol = Chem.MolFromSmiles(smiles_string)
                mol = Chem.AddHs(mol)

                amideGroups = addAmides(mol)
                for conformation_xyz in natsorted(os.listdir(f"{name}_Conformations")):
                    if conformation_xyz.endswith('.xyz'):  # working within 1 conformer
                        mol.RemoveAllConformers()
                        mol = load_xyz_coords(mol, f"{working_dir}/{name}_Conformations/{conformation_xyz}")

                        #here, put the function of what you want to calcualte for each conformation
                        peptide_descriptors.append(getAmideDistances(amideGroups,mol))

                peptide_boltzmann = boltzmann(peptide_descriptors,working_dir,name)
                peptide_boltzmann = peptide_boltzmann.reshape(len(peptide_boltzmann),-1)
                df = pd.DataFrame(peptide_boltzmann)
                df.to_csv(f'{descriptor_name}.csv', index=False, header=False)
            os.chdir(main_dir)

def smart_parse_token(token: str):
    try:
        return float(token)
    except ValueError:
        return token

def read_file(file_name) -> list[str]:
    with open(file_name, 'r') as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines]
        if '\t' in lines[0]: #if its gonna be a list of lists
            lines = [[smart_parse_token(x) for x in line.strip().split('\t')] for line in lines]
        return lines



def remove_duplicates(smiles_lines_all,names_lines_all,percents_all) -> tuple[list,list,list]:
    smile_to_indices = defaultdict(list)
    for i, smile in enumerate(smiles_lines_all):
        smile_to_indices[smile].append(i)


    indices_to_keep = set()
    for indices in smile_to_indices.values():
        if len(indices) == 1:
            indices_to_keep.add(indices[0])
        else:
            best_index = max(indices, key=lambda x: percents_all[x][0]/100) #first index of the [6,12,18]
            indices_to_keep.add(best_index)

    indices_to_remove = [i for i in range(len(smiles_lines_all)) if i not in indices_to_keep]

    smiles_lines_all = [j for i,j in enumerate(smiles_lines_all) if i not in indices_to_remove]
    names_lines_all = [j for i,j in enumerate(names_lines_all) if i not in indices_to_remove]
    percents_all = [j for i,j in enumerate(percents_all) if i not in indices_to_remove]
    return smiles_lines_all, names_lines_all, percents_all


def main():
    main_dir = "/Users/zaan/zasaeed@g.hmc.edu - Google Drive/Shared drives/OMO Lab/Projects/OMO Lab - Zaan Saeed/Data/Peptides"
    os.chdir(main_dir)

    smiles_lines_all = read_file("all_peptides.smi") #list of smiles
    names_lines_all = read_file("all_names.txt") # list of names
    percents_all = read_file("percent6-12-18.txt")# list of lists of percents
    smiles_lines_all = sort_by_names_alphabetically(names_lines_all,smiles_lines_all) #sorts bynames
    percents_all = sort_by_names_alphabetically(names_lines_all,percents_all) #sorts by names
    names_lines_all = sort_by_names_alphabetically(names_lines_all,names_lines_all) #changes order of names_lines, must be last


    smiles_final, names_final, percents_final = remove_duplicates(smiles_lines_all,names_lines_all,percents_all) #sorted, with no duplicates

    X = create_X(main_dir,names_final,"BWdistances")
    Y = create_Y(percents_final)

    ### removal of 1s and 0s

    X = X.tolist()
    Y =Y.tolist()

    for i in reversed(range(len(X))):
        if Y[i] == 0 or Y[i] == 1:
            del X[i]
            del Y[i]
    X = np.array(X)
    Y = np.array(Y)
    plot_Y_distribution(Y)
    run_RFR(X, Y, 0.3, 5)



if __name__ == "__main__":
    main()









