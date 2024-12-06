from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors
from rdkit.Chem.Draw import SimilarityMaps

# Step 1: Create the molecule
smiles = 'CC(C)C[C@H](C(N[C@H](Cc1ccccc1)C(N[C@H](CC(C)C)C(NC(C)(C)C(N[C@H](CC(C)C)C(N[C@@H](CC(O[n]1nnc2cccnc12)=O)C(NC(C)c1ccccc1)=O)=O)=O)=O)=O)=O)N'
mol = Chem.MolFromSmiles(smiles)

# Step 2: Compute Gasteiger charges as a descriptor
AllChem.ComputeGasteigerCharges(mol)
contribs = [mol.GetAtomWithIdx(i).GetDoubleProp('_GasteigerCharge') for i in range(mol.GetNumAtoms())]

# Step 3: Visualize the similarity map based on Gasteiger charges
d2d = Draw.MolDraw2DCairo(800, 800)  # Increase the canvas size (800x800 pixels)

# Step 4: Generate the similarity map
_ = SimilarityMaps.GetSimilarityMapFromWeights(
    mol,
    contribs,
    draw2d=d2d,
    colorMap='jet',
    contourLines=10
)

# Step 5: Save the result
d2d.FinishDrawing()
with open("descriptor_visualization_fixed.png", "wb") as f:
    f.write(d2d.GetDrawingText())
print("Image saved as 'descriptor_visualization_fixed.png'")

vals = Descriptors.CalcMolDescriptors(mol)
print(vals['TPSA'])