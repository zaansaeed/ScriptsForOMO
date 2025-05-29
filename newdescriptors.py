from array import array
from venv import create

import numpy as np
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, rdFreeSASA, ChemicalFeatures
from natsort import natsorted
from ML_functions import *
from functions_for_smile_to_xyz import boltzmann

def compute_global_descriptors(mol):
    """
    Computes whole-molecule descriptors (2D + 3D) for a single Mol.
    """
    desc = {}
    desc['MolLogP'] = Descriptors.MolLogP(mol)
    desc['MolMR'] = Descriptors.MolMR(mol)
    desc['TPSA'] = rdMolDescriptors.CalcTPSA(mol)
    desc['NumHAcceptors'] = Descriptors.NumHAcceptors(mol)
    desc['NumHDonors'] = Descriptors.NumHDonors(mol)

    # --- STERIC (2D) ---
    desc['FractionCSP3'] = Descriptors.FractionCSP3(mol)
    desc['NumRotatableBonds'] = Descriptors.NumRotatableBonds(mol)
    desc['LabuteASA'] = Descriptors.LabuteASA(mol)
    desc['BertzCT'] = Descriptors.BertzCT(mol)
    desc['Kappa1'] = Descriptors.Kappa1(mol)
    desc['Kappa2'] = Descriptors.Kappa2(mol)
    desc['Kappa3'] = Descriptors.Kappa3(mol)
    desc['BalabanJ'] = Descriptors.BalabanJ(mol)

    # --- STERIC (3D) ---
    desc['RadiusOfGyration'] = rdMolDescriptors.CalcRadiusOfGyration(mol)
    desc['InertialShapeFactor'] = rdMolDescriptors.CalcInertialShapeFactor(mol)
    desc['PMI1'] = rdMolDescriptors.CalcPMI1(mol)
    desc['PMI2'] = rdMolDescriptors.CalcPMI2(mol)
    desc['PMI3'] = rdMolDescriptors.CalcPMI3(mol)
    desc['NPR1'] = rdMolDescriptors.CalcNPR1(mol)
    desc['NPR2'] = rdMolDescriptors.CalcNPR2(mol)
    desc['Asphericity'] = rdMolDescriptors.CalcAsphericity(mol)
    desc['Eccentricity'] = rdMolDescriptors.CalcEccentricity(mol)
    desc['SpherocityIndex'] = rdMolDescriptors.CalcSpherocityIndex(mol)

    return list(desc.values())

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
            for conformation_xyz in natsorted(os.listdir(f"{name}_Conformations")):
                if conformation_xyz.endswith('.xyz'):  # working within 1 conformer
                    mol = Chem.MolFromSmiles(smiles_string)
                    mol = Chem.AddHs(mol)
                    mol = load_xyz_coords(mol, f"{working_dir}/{name}_Conformations/{conformation_xyz}")
                    peptide_descriptors.append(compute_global_descriptors(mol))
            peptide_boltzmann= boltzmann(peptide_descriptors,working_dir,name)
            print(peptide_boltzmann)
            all_records.append(peptide_boltzmann)
            os.chdir(main_dir)
    all_records = np.array(all_records)
    df = pd.DataFrame(all_records)
    df.to_csv(main_dir+'/boltzmann_descriptors.csv', index=False, header=False)

X_22 = []
data = pd.read_csv(main_dir+'/boltzmann_descriptors.csv',header=None,index_col=None)
for d in data.values.tolist():
    X_22.append(d)

X_22 = np.array(X_22)
os.chdir(main_dir)

Y = X_22[:,13]

with open("all_peptides.smi", "r") as f:
    smiles_lines = f.readlines()
    smiles_lines = [line.strip() for line in smiles_lines]

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

indices_to_remove = [i for i in range(len(X_22)) if i not in indices_to_keep]
Y = [j for i, j in enumerate(Y) if i not in indices_to_remove]

X_22 =np.delete(X_22,13,1)

X_dihedrals = create_X(main_dir, "BWDihedralNormalized") #ready for input
X_distances = create_X(main_dir, "BWdistances")

X_dihedrals = X_dihedrals.reshape(len(X_dihedrals), -1)
X_distances = X_distances.reshape(len(X_distances), -1)

X = np.hstack((X_dihedrals, X_distances,X_22))

X = [j for i, j in enumerate(X) if i not in indices_to_remove]
X=np.array(X)


print(len(X),len(X[0]),len(Y))


import matplotlib.pyplot as plt

run_RFR(X,Y)


for x in range(X.shape[1]):
    feature_values = X[:, x]  # shape (N,)

    plt.scatter(feature_values, Y, alpha=0.6)
    plt.xlabel(f'Feature at index {x}')
    plt.ylabel('Target variable')
    plt.title(f'Scatter plot of Feature {x} vs Target')
    plt.show()



