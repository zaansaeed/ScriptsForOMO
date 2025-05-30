
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, rdFreeSASA, ChemicalFeatures
from scipy import stats
import pandas as pd
import numpy as np
from ML_functions import *
from natsort import natsorted
from functions_for_smile_to_xyz import boltzmann
from functions_for_smile_to_xyz import AmideGroup, addAmides

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



main_dir = os.path.abspath("/Users/zaan/zasaeed@g.hmc.edu - Google Drive/Shared drives/OMO Lab/Projects/OMO Lab - Zaan Saeed/Data/Peptides")
os.chdir(main_dir)

for folder in natsorted(os.listdir(main_dir)):
    if os.path.isdir(folder):
        os.chdir(folder) #currenlty working in Peptide_{name}
        if os.path.exists("NewBWDistances.csv"):
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
                    print(conformation_xyz)
                    mol.RemoveAllConformers()
                    mol = load_xyz_coords(mol, f"{working_dir}/{name}_Conformations/{conformation_xyz}")
                    peptide_descriptors.append(getAmideDistances(amideGroups,mol))
            peptide_boltzmann = boltzmann(peptide_descriptors,working_dir,name)
            peptide_boltzmann = peptide_boltzmann.reshape(len(peptide_boltzmann),-1)
            df = pd.DataFrame(peptide_boltzmann)
            df.to_csv('NewBWDistances.csv', index=False, header=False)
        os.chdir(main_dir)






X = create_X(main_dir, "BWDihedralNormalized")

percents = create_Y(main_dir)

os.chdir(main_dir)
with open("all_peptides.smi", "r") as f:
    smiles_lines = f.readlines()
    smiles_lines = [line.strip() for line in smiles_lines]

with open("all_names.txt", "r") as f:
    names_lines = f.readlines()
    names_lines = [name.strip() for name in names_lines]

names_percents_dictionary = dict(zip(names_lines,percents))
names_percents_dictionary = dict((k,names_percents_dictionary[k]) for k in natsorted(names_percents_dictionary))

Y = []
for value in names_percents_dictionary.values():
    Y.append(value)

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


X = [j for i, j in enumerate(X) if i not in indices_to_remove]
Y = [j for i, j in enumerate(Y) if i not in indices_to_remove]

for i in reversed(range(len(X))):
    if  Y[i]==0 or Y[i]==1:
        del X[i]
        del Y[i]

X =np.array(X)
Y=np.array(Y)
#run_RFC(X,Y)
run_RFR(X,Y)
import matplotlib.pyplot as plt

plt.hist(Y, bins=50, edgecolor='k')
plt.title('Distribution of Y values')
plt.xlabel('Y')
plt.ylabel('Frequency')
plt.show()
plt.clf()
