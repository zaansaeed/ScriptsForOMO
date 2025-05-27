from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor, GradientBoostingRegressor, HistGradientBoostingRegressor
from sklearn.model_selection import train_test_split, GridSearchCV, RandomizedSearchCV ,cross_val_score, KFold
from sklearn.metrics import accuracy_score, classification_report, mean_squared_error, r2_score, make_scorer, \
     f1_score,mean_absolute_error, r2_score, root_mean_squared_error
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.svm import SVR, SVC
from sklearn.svm import SVC
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
import numpy as np
from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors, rdFreeSASA, ChemicalFeatures
from scipy import stats
import pandas as pd
from ML_functions import *
from functions_for_smile_to_xyz import boltzmann

def compute_global_descriptors(mol):
    """
    Computes whole-molecule descriptors (2D + 3D) for a single Mol.
    """
    descriptors=  {
        # 2D topological/descriptive
        'MolWt':           Descriptors.MolWt(mol),
        'LogP':            Descriptors.MolLogP(mol),
        'TPSA':            rdMolDescriptors.CalcTPSA(mol),
        'NumHDonors':      rdMolDescriptors.CalcNumHBD(mol),
        'NumHAcceptors':   rdMolDescriptors.CalcNumHBA(mol),
        'NumRotatable':    rdMolDescriptors.CalcNumRotatableBonds(mol),
        'NumRings':        rdMolDescriptors.CalcNumRings(mol),
        'LabuteASA':       rdMolDescriptors.CalcLabuteASA(mol),
        # 3D conformation‐dependent
        'RadiusOfGyration': rdMolDescriptors.CalcRadiusOfGyration(mol, confId=0),
        'Eccentricity':     rdMolDescriptors.CalcEccentricity(mol, confId=0),
        'InertiaTensor1':   rdMolDescriptors.CalcInertialShapeFactor(mol, confId=0),
    }
    return list(descriptors.values())

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
if not os.path.exists(main_dir+'/boltzmann_descriptors.csv'):
    for folder in os.listdir(main_dir):
        if os.path.isdir(folder) :
            os.chdir(folder)
            working_dir = os.getcwd()
            name = folder.split("_")[1]
            smiles_string = open(f"{name}.smi").read().strip()
            print(name)
            peptide_descriptors = []
            for conformation_xyz in os.listdir(f"{name}_Conformations"):
                if conformation_xyz.endswith('.xyz'):  # working within 1 conformer
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