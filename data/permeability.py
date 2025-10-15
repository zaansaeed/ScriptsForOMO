import pandas as pd
from rdkit import Chem

# Load and preprocess data
data = pd.read_csv('permeability_data.csv')
data = data[["ID", "HELM", "Permeability"]]
data["HELM"] = data["HELM"].str.extract(r"\{([^}]*)\}")

# Prepare output lists
ids = []
smiles_list = []
permeabilities = []

for id, helm_seq, permeability in zip(data["ID"], data["HELM"], data["Permeability"]):
    helm_string = f"PEPTIDE1{{{helm_seq}}}$$$$"
    
    try:
        mol = Chem.MolFromHELM(helm_string)  # replace with your own HELM parser if needed
    except:
        mol = None

    if mol is not None:
        try:
            Chem.SanitizeMol(mol)
            smiles = Chem.MolToSmiles(mol)
            ids.append(id)
            smiles_list.append(smiles)
            permeabilities.append(permeability)
        except:
            continue  # skip invalid molecules

# Create and save final DataFrame
output_df = pd.DataFrame({
    "ID": ids,
    "SMILES": smiles_list,
    "Permeability": permeabilities
})

output_df.to_csv("helm_to_smiles_with_permeability.csv", index=False)
