import numpy as np
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, rdFreeSASA, ChemicalFeatures
from scipy import stats
import pandas as pd
from ML_functions import *
from natsort import natsorted

from ScriptsForOMO.ML_functions import create_Y
from functions import boltzmann

from mordred import Calculator, descriptors

import numbers
def compute_global_descriptors(mol):
    calc = Calculator(descriptors,ignore_3D=False )
    calc.descriptors = [d for d in calc.descriptors if d.require_3D == True]
    desc = calc(mol)
    desc_dict = desc.asdict()
    clean_desc = {k: v for k,v in desc_dict.items() if isinstance(v,numbers.Number)}
    return list(clean_desc.values())

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



main_dir = os.path.abspath("/Users/zaansaeed/Peptides")
os.chdir(main_dir)
all_records = []
if not os.path.exists(main_dir+'/boltzmann_descriptors.csv'):
    for folder in natsorted(os.listdir(main_dir)):
        if os.path.isdir(folder) :
            os.chdir(folder)
            working_dir = os.getcwd()
            name = folder.split("_")[1]
            smiles_string = open(f"{name}.smi").read().strip()
            print(name)
            peptide_descriptors = []
            mol = Chem.MolFromSmiles(smiles_string)
            mol = Chem.AddHs(mol)
            for conformation_xyz in natsorted(os.listdir(f"{name}_Conformations")):
                if conformation_xyz.endswith('.xyz'):  # working within 1 conformer
                    print(conformation_xyz)
                    mol.RemoveAllConformers()
                    mol = load_xyz_coords(mol, f"{working_dir}/{name}_Conformations/{conformation_xyz}")
                    peptide_descriptors.append(compute_global_descriptors(mol))
            peptide_boltzmann= boltzmann(peptide_descriptors,working_dir,name)
            print(peptide_boltzmann)
            all_records.append(peptide_boltzmann)
            os.chdir(main_dir)
    all_records = np.array(all_records)
    df = pd.DataFrame(all_records)
    df.to_csv(main_dir+'/boltzmann_descriptors.csv', index=False, header=False)

X = []
data = pd.read_csv(main_dir+'/boltzmann_descriptors.csv',header=None,index_col=None)
for d in data.values.tolist():
    X.append(d)
X = np.array(X)

percents = create_Y(main_dir)

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
   if Y[i]==1 or Y[i]==0:
        del X[i]
        del Y[i]


X = np.array(X)


run_RFR(X,Y,.2)