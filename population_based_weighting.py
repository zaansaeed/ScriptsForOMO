
from rdkit import Chem
from rdkit.Chem import rdMolAlign
import numpy as np
import os
from natsort import natsorted
import concurrent.futures




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


n_conformers = 393  # Adjust based on your number of conformers
rmsd_matrix = np.zeros((n_conformers, n_conformers))  # This will store the RMSD values in a lower triangular matrix

def calculate_rmsd(i, j):
    # Calculate RMSD for conformers i and j
    rdMolAlign.AlignMol(mol,mol,i,j)
    rmsd = rdMolAlign.GetBestRMS(mol, mol, i, j,numThreads=10,)
    return rmsd

# Using ThreadPoolExecutor to parallelize RMSD calculation

for i in range(n_conformers):
    for j in range(i):
        rsmd = calculate_rmsd(i, j)
        rmsd_matrix[i, j] = rsmd
        print(f"conformation {i}-{j}: {rsmd}")
np.savetxt("rmsd_matrix.txt", rmsd_matrix, delimiter=",")
