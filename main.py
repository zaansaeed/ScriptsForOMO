import os

from sklearn.preprocessing import StandardScaler

import functions_for_smile_to_xyz as fs
import ML_functions as ML
from natsort import natsorted
import numpy as np

schrodinger_path ="/opt/schrodinger/suites2024-3/"
main_dir = os.path.abspath("/Users/zaansaeed/Peptides")
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

Y = ML.create_Y(main_dir)
from collections import defaultdict

name_to_indices = defaultdict(list)
for i, name in enumerate(smiles_lines):
    name_to_indices[name].append(i)

indices_to_keep = set()
for indices in name_to_indices.values():
    if len(indices) == 1:
        indices_to_keep.add(indices[0])
    else:
        best_index = max(indices, key=lambda x: Y[x])
        indices_to_keep.add(best_index)

indices_to_remove = [i for i in range(len(Y)) if i not in indices_to_keep]



update_matrices = 0


for i, name in enumerate(names_lines): #processing : smiles -> xyzs
    working_dir = main_dir+f"/Peptide_{name}"

    if update_matrices == 1:
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

X_dihedrals = ML.create_X(main_dir, "BWDihedralNormalized") #ready for input
X_distances = ML.create_X(main_dir, "BWdistances")
scalar = StandardScaler()


X_dihedrals = X_dihedrals.reshape(len(X_dihedrals), -1)
X_distances = X_distances.reshape(len(X_distances), -1)


X_combined = np.hstack((X_distances,X_dihedrals))

X_distances = [X for i, X in enumerate(X_distances) if i not in indices_to_remove]
Y = [Y for i, Y in enumerate(Y) if i not in indices_to_remove]

for i in reversed(range(len(X_distances))):
    if Y[i]==0 or Y[i]==1:
        del X_distances[i]
        del Y[i]
X_distances =np.array(X_distances)




#Y = ML.create_YC(Y,.7 )



#ML.run_RFC(X_dihedrals,Y)
#ML.run_SVM(X,Y)

ML.run_RFR(X_distances,Y)

import matplotlib.pyplot as plt
for x in range(X_distances.shape[1]):
    feature_values = X_distances[:, x]  # shape (N,)

    plt.scatter(feature_values, Y, alpha=0.6)
    plt.xlabel(f'Feature at index {x}')
    plt.ylabel('Target variable')
    plt.title(f'Scatter plot of Feature {x} vs Target')
    plt.show()
