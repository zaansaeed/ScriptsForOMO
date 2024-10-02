from rdkit import Chem
from ase import Atoms
from ase.io import write


def sdf_to_pdb(nameOfJob, sdf_file, pdb_folder):
    #molecules contained in sdf file
    supplier = Chem.SDMolSupplier(sdf_file)
    #interating through every molecule in sdf file
    for i, mol in enumerate(supplier):
        if mol is None: continue
        print(f"Molecule: {i+1} has {mol.GetNumAtoms()} atoms" )

        Chem.MolToPDBFile(mol,pdb_folder + "/" + nameOfJob +"_Molecule_" + str(i)+ ".pdb")

#path of SDF file data
sdf = ''
#path of folder where to put the XYZ files
pdb = ''

sdf_to_pdb("LinearPeptide_to_PDB",sdf,pdb)