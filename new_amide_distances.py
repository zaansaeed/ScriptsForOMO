
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
            if i == j and amid1.getH() is not None:
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



main_dir = os.path.abspath("/Users/zaansaeed/Peptides")
os.chdir(main_dir)
all_records = []
if not os.path.exists(main_dir+'/amide_distances.csv'):
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

            amideGroups = addAmides(mol)
            for conformation_xyz in natsorted(os.listdir(f"{name}_Conformations")):
                if conformation_xyz.endswith('.xyz'):  # working within 1 conformer
                    print(conformation_xyz)
                    mol.RemoveAllConformers()
                    mol = load_xyz_coords(mol, f"{working_dir}/{name}_Conformations/{conformation_xyz}")
                    peptide_descriptors.append(getAmideDistances(amideGroups,mol))
            peptide_boltzmann = boltzmann(peptide_descriptors,working_dir,name)
            print(peptide_boltzmann)
            print(peptide_boltzmann.shape)
            all_records.append(peptide_boltzmann)
            os.chdir(main_dir)
    all_records = np.array(all_records)
    all_record = np.reshape(len(all_records),-1)
    df = pd.DataFrame(all_records)
    df.to_csv(main_dir+'/boltzmann_descriptors.csv', index=False, header=False)

X = []
data = pd.read_csv(main_dir+'/boltzmann_descriptors.csv',header=None,index_col=None)
for d in data.values.tolist():
    X.append(d)
Y = create_Y(main_dir)


run_SVR(np.array(X),Y)