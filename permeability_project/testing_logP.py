import rdkit
from rdkit import Chem
import rdkit.Chem.Descriptors as Descriptors
import pandas as pd
import numpy as np
import os, sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import functions as funcs

from visualization import visualize_molecule_with_highlights

pdb_file_of_template = os.getcwd()+"/debugging/template-1.pdb"
peptide = funcs.add_double_bonds_to_pdb(pdb_file_of_template)
amide_groups = funcs.add_amides(peptide)
pb = Chem.MolToMolBlock(peptide)



group = amide_groups[-1]
residue_ids = group.getResidue2()


residue = funcs.make_submol_from_atom_ids(peptide, residue_ids)
for atom in residue.GetAtoms():
    if atom.GetIsAromatic() and not atom.IsInRing():
        atom.SetIsAromatic(False)
Chem.SanitizeMol(residue)


results = Descriptors.MolLogP(residue)
print(f"[Info] logP of residue: {results}")
