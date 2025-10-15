import pandas as pd
from lxml import etree
from rdkit import Chem
from rdkit.Chem import rdChemReactions
from rdkit.Chem import molzip
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import Draw
import copy
import ast

data = pd.read_csv("monomer_list.csv")

mol_dict = {}
for symbol, smiles,typ in zip(data["Symbol"], data["CXSMILES"],data["Monomer_Type"]):
    mol = Chem.MolFromSmiles(smiles)
    mol.SetProp("Type",str(typ))
    if mol is not None:
        mol_dict[symbol] = mol

print(len(mol_dict))
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')



def combine_fragments(m1, m2):
    m1 = Chem.Mol(m1)
    m2 = Chem.Mol(m2)
    for atm in m1.GetAtoms():
        if atm.HasProp('atomLabel') and atm.GetProp('atomLabel') == '_R2':
            atm.SetAtomMapNum(1)
    for atm in m2.GetAtoms():
        if atm.HasProp('atomLabel') and atm.GetProp('atomLabel') == '_R1':
            atm.SetAtomMapNum(1)
    return molzip(m1, m2)

def make_peptide(sequence):
    sequence = ast.literal_eval(sequence)
    monomerlist = []
    for sym in sequence:
        monomerlist.append(mol_dict[sym])
    for idx, monomer in enumerate(monomerlist):
        if Chem.MolToSmiles(monomer).count("*") == 1:
            continue
        if idx==0:
            res = monomer
        else:
            res = combine_fragments(res, monomer)
    return res

def cap_terminal(m):
    m = Chem.Mol(m)
    n_term = Chem.MolFromSmiles('[H][*:1]')
    c_term = Chem.MolFromSmiles('CO[*:2]')
    for atm in m.GetAtoms():
        if atm.HasProp('atomLabel') and atm.GetProp('atomLabel') == '_R1':
            atm.SetAtomMapNum(1)
        if atm.HasProp('atomLabel') and atm.GetProp('atomLabel') == '_R2':
            atm.SetAtomMapNum(2)
    res = molzip(m, n_term)
    res = molzip(res, c_term)
    return res


data = pd.read_csv("permeability_data.csv")

data = data[["ID","Sequence","Permeability"]]

proteins = []

for ID, sequence, permeability in zip(data["ID"], data["Sequence"],data["Permeability"]):
    mol = make_peptide(sequence)
    mol = cap_terminal(mol)
    Chem.SanitizeMol(mol)
    smiles = Chem.MolToSmiles(mol)
    proteins.append({
        "ID": ID,
        "Sequence": sequence,
        "Permeability": permeability,
        "SMILES": smiles
    })
df = pd.DataFrame(proteins)
df.to_csv("processed_peptides.csv", index=False)



