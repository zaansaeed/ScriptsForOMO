from array import array
from tkinter.font import names

import numpy as np
from rdkit import Chem
from ML_functions import *
from natsort import natsorted
from collections import defaultdict
from rdkit.Chem import Descriptors3D, Descriptors, Crippen, Lipinski
import inspect
from functions import get_amide_distances
from functions import boltzmann
from functions import add_amides

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
            peptide_name  = folder.split('_')[1]
            if not os.path.exists(f"{peptide_name}_{descriptor_name}.csv"): #change to/from not
                working_dir = os.getcwd()
                name = folder.split("_")[1]
                smiles_string = open(f"{name}.smi").read().strip()
                print(name)
                peptide_descriptors = []
                mol = Chem.MolFromSmiles(smiles_string)
                mol = Chem.AddHs(mol)

                for conformation_xyz in natsorted(os.listdir(f"{name}_Conformations")):
                    if conformation_xyz.endswith('.xyz'):  # working within 1 conformer
                        mol.RemoveAllConformers()
                        mol = load_xyz_coords(mol, f"{working_dir}/{name}_Conformations/{conformation_xyz}")
                        #here, put the function of what you want to calcualte for each conformation
                        peptide_descriptors.append(compute_global_descriptors(mol))

                peptide_boltzmann = boltzmann(peptide_descriptors,working_dir,name)
                print(peptide_boltzmann.shape)
                peptide_boltzmann = peptide_boltzmann.reshape(len(peptide_boltzmann),-1)
                print(peptide_boltzmann)
                df = pd.DataFrame(peptide_boltzmann)
                df.to_csv(f'{peptide_name}_{descriptor_name}.csv', index=False, header=False)
            os.chdir(main_dir)

def side_chain_descriptors(amidegroups):
    descriptors_for_peptide = []

    for amide_group in amidegroups:
        mol = amide_group.getResidue1()[0]
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic() and not atom.IsInRing():
                atom.SetIsAromatic(False)
        Chem.SanitizeMol(mol)
        results = {
            "TPSA": Descriptors.TPSA(mol),  # Topological Polar Surface Area
            "Radius" :Descriptors3D.RadiusOfGyration(mol)

        }
        descriptors_for_peptide.append(list(results.values()))
        if amide_group.getResidue2() is not None:## CHECKING THE LAST AMIDE GROUP, IT IS THE ONLY ONE WITH 2 RESIDUES
            mol = amide_group.getResidue2()[0]
            for atom in mol.GetAtoms():
                if atom.GetIsAromatic() and not atom.IsInRing():
                    atom.SetIsAromatic(False)
            Chem.SanitizeMol(mol)
            results = {
                "TPSA": Descriptors.TPSA(mol),  # Topological Polar Surface Area
                "Radius": Descriptors3D.RadiusOfGyration(mol)

            }
            descriptors_for_peptide.append(list(results.values()))
    print(len(descriptors_for_peptide))
    return descriptors_for_peptide




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

#computes your specified descriptor for each conformation
def compute_global_descriptors(mol):
    descriptor_funcs = {
        "Radius": Descriptors3D.RadiusOfGyration(mol)
    }
    return list(descriptor_funcs.values())

def main():
    main_dir = "/Users/zaansaeed/Peptides"
    os.chdir(main_dir)

    create_new_descriptor('RadiusOfGyration', main_dir)



    smiles_lines_all = read_file("all_peptides.smi") #list of smiles
    names_lines_all = read_file("all_names.txt") # list of names
    percents_all = read_file("percent6-12-18.txt")# list of lists of percents
    smiles_lines_all = sort_by_names_alphabetically(names_lines_all,smiles_lines_all) #sorts by names
    percents_all = sort_by_names_alphabetically(names_lines_all,percents_all) #sorts by names
    names_lines_all = sort_by_names_alphabetically(names_lines_all,names_lines_all) #changes order of names_lines, must be last


    smiles_final, names_final, percents_final = remove_duplicates(smiles_lines_all,names_lines_all,percents_all) #sorted, with no duplicates

    X = create_X(main_dir,names_final,["BWDihedralNormalized","BWdistances"])

    Y = create_Y_ROG(main_dir,names_final)
    print(X.shape,Y.shape)
    ### removal of 1s and 0s



    print(X.shape,Y.shape)
    plot_Y_distribution(Y)

    #run_elasticnet(X,Y,5,0.2)
    #run_RFR(X,Y,5,0.2)


    run_SVR(X,Y,5,0.2)
    #run_NN(X,Y,0.3,5)



if __name__ == "__main__":
    main()



"""
    n_bins = 10
    bin_edges = np.linspace(Y.min(), Y.max(), n_bins + 1)
    y_bins = np.digitize(Y, bin_edges) - 1  # -1 to make bins 0-indexed
    # Set maximum samples per bin
    max_per_bin = 6
    keep_indices = []

    # Process each bin
    for bin_num in range(n_bins):
        bin_mask = (y_bins == bin_num)
        bin_indices = np.where(bin_mask)[0]

        if len(bin_indices) > max_per_bin:
            # Randomly sample from this bin
            np.random.seed(42)  # For reproducibility
            selected_indices = np.random.choice(bin_indices, size=max_per_bin, replace=False)
        else:
            # Keep all samples in this bin
            selected_indices = bin_indices

        keep_indices.extend(selected_indices)

    # Convert to numpy array and sort
    keep_indices = np.array(keep_indices)
    keep_indices = np.sort(keep_indices)

    X_test_indices = np.array([8, 13, 44, 79, 75, 98, 90, 94, 65])
    # Ensure test set is separate from training set
    X_test = X[X_test_indices]
    Y_test = Y[X_test_indices]
    # Create balanced dataset
    X = X[keep_indices]
    Y = Y[keep_indices]

    X = np.array(X)
    X_test = np.array(X_test)
    Y = np.array(Y)
    Y_test = np.array(Y_test)
"""





