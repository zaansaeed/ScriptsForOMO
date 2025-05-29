
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, rdFreeSASA, ChemicalFeatures
from scipy import stats
import pandas as pd
from ML_functions import *
from natsort import natsorted
from functions_for_smile_to_xyz import boltzmann

from mordred import Calculator, descriptors

import numbers
def compute_global_descriptors(mol):
    calc = Calculator(descriptors,ignore_3D=False )
    calc.descriptors = [d for d in calc.descriptors if d.require_3D == True]
    desc = calc(mol)
    desc_dict = desc.asdict()
    clean_desc = {k: v for k,v in desc_dict.items() if isinstance(v,numbers.Number)}
    print(len(clean_desc))
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



main_dir="/Users/zaan/zasaeed@g.hmc.edu - Google Drive/Shared drives/OMO Lab/Projects/OMO Lab - Zaan Saeed/Data/Peptides"
os.chdir(main_dir)
all_records = []
if os.path.exists(main_dir+'/boltzmann_descriptors.csv'):
    for folder in os.listdir(main_dir):
        if os.path.isdir(folder) :
            os.chdir(folder)
            working_dir = os.getcwd()
            name = folder.split("_")[1]
            smiles_string = open(f"{name}.smi").read().strip()
            print(name)
            peptide_descriptors = []
            for conformation_xyz in natsorted(os.listdir(f"{name}_Conformations")):
                if conformation_xyz.endswith('.xyz'):  # working within 1 conformer
                    print(conformation_xyz)
                    mol = Chem.MolFromSmiles(smiles_string)
                    mol = Chem.AddHs(mol)
                    mol = load_xyz_coords(mol, f"{working_dir}/{name}_Conformations/{conformation_xyz}")
                    compute_global_descriptors(mol)
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
print(X[0])
Y = six_over_target_percents(create_outputs(main_dir))
print(Y)

run_SVR(np.array(X),Y)