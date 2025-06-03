from array import array
from venv import create

import numpy as np
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, rdFreeSASA, ChemicalFeatures
from natsort import natsorted
from ML_functions import *
from functions import boltzmann

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




