from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit import Chem
from rdkit.Chem import rdDistGeom
from rdkit.Chem import rdMolAlign
import rdkit
import os
from natsort import natsorted



def load_xyz_coords(mol, xyz_path) -> None:
    """
    Loads atomic coordinates for a molecular structure from an XYZ file and updates
    the input molecule object with the new coordinates. The function extracts
    coordinates from the XYZ file, constructs a conformer, and assigns it to the
    input RDKit molecule.

    :param mol: The molecule object that will have its coordinates updated.
    :type mol: Chem.Mol
    :param xyz_path: Path to the XYZ file containing molecular coordinates in the format
                     where the first two lines are skipped, and subsequent lines contain
                     atomic symbol followed by x, y, z coordinates.
    :type xyz_path: str
    :return: The input molecule with updated coordinates.
    :rtype: Chem.Mol
    """
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

    mol.AddConformer(conf, assignId=True)

def add_conformations_to_mol(mol,conformation_folder):
    for conformation in natsorted(os.listdir(conformation_folder)):
        load_xyz_coords(mol, os.path.join(conformation_folder, conformation))





mol = Chem.MolFromPDBFile('/Users/zaan/Library/CloudStorage/GoogleDrive-zasaeed@g.hmc.edu/Shared drives/OMO Lab/Projects/OMO Lab - ML - Zaan Saeed/Data/Peptides/Peptide_C7R1/C7R1-out-template.pdb',removeHs=False)
add_conformations_to_mol(mol,'/Users/zaan/Library/CloudStorage/GoogleDrive-zasaeed@g.hmc.edu/Shared drives/OMO Lab/Projects/OMO Lab - ML - Zaan Saeed/Data/Peptides/Peptide_C7R1/C7R1_Conformations')


dists = []

import concurrent.futures

def calculate_rmsd(i, j):
    rmsd = rdMolAlign.GetBestRMS(mol, mol, i, j)  # align the i-th conformation to the j-th conformation
    return (i, j, rmsd)


# Using ThreadPoolExecutor to parallelize RMSD calculation
"""with concurrent.futures.ThreadPoolExecutor() as executor:
    # Submit tasks for all pairs of conformers
    futures = [executor.submit(calculate_rmsd, i, j) for i in range(393) for j in range(i)]

    # Retrieve results as they finish
    for future in concurrent.futures.as_completed(futures):
        i, j, rmsd = future.result()
        dists.append(rmsd)
        print(f"RMSD between conformation {i + 1} and conformation {j + 1}: {rmsd}")
"""

conformer_10 = mol.GetConformer(10)
print(conformer_10.GetPositions())