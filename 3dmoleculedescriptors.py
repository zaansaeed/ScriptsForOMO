from rdkit.Chem import AllChem, Draw, Descriptors3D, Descriptors
from rdkit.Chem import rdMolDescriptors
from rdkit import Chem
import sys, py3Dmol
from rdkit.Chem.Draw import IPythonConsole



m = Chem.MolFromSmiles('CC(C)C[C@H](C(N[C@H](Cc1ccccc1)C(N[C@H](CC(C)C)C(NC(C)(C)C(N[C@H](CC(C)C)C(N[C@@H](CC(O[n]1nnc2cccnc12)=O)C(NC(C)c1ccccc1)=O)=O)=O)=O)=O)=O)N')
mol=Chem.AddHs(m)
AllChem.EmbedMolecule(mol)

vals = Descriptors3D.CalcMolDescriptors3D(mol)
