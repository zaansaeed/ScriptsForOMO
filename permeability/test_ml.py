import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
import functions as funcs
import ML_functions as ML

import numpy as np

data = pd.read_csv('helm_to_smiles_with_permeability.csv')

array_of_features = []

Y= []  # Reshape Y to be a 2D array with one column
for id, smiles, permeability in zip(data["ID"], data["SMILES"], data["Permeability"]):
     mol = Chem.MolFromSmiles(smiles)
     mol = Chem.AddHs(mol)
     Chem.SanitizeMol(mol)
     AllChem.EmbedMolecule(mol)
     AllChem.MMFFOptimizeMolecule(mol)
     try:
          amide_groups = funcs.add_amides(mol)
          array_of_features.append(funcs.side_chain_descriptors(amide_groups, mol))
          Y.append(permeability)  # Append permeability to Y
          print(f"Processed {id} with permeability {permeability}")
     except Exception as e:
          print(f"Error processing {id}: {e}")

X = np.array(array_of_features)
Y = np.array(Y).reshape(-1, 1)  # Reshape Y to be a 2D array with one column
X = X.reshape(X.shape[0], -1)  # Flatten the features if necessary
print(X.shape, Y.shape)
ML.run_ElasticNet(X, Y)


