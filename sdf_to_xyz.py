from rdkit import Chem
from ase import Atoms
from ase.io import write


def sdf_to_xyz(nameOfJob, sdf_file, xyz_folder):
    #molecules contained in sdf file
    supplier = Chem.SDMolSupplier(sdf_file)
    #interating through every molecule in sdf file
    for i, mol in enumerate(supplier):
        if mol is None: continue
        print(f"Molecule: {i+1} has {mol.GetNumAtoms()} atoms" )

        atoms_symbols = [atom.GetSymbol() for atom in mol.GetAtoms()]  # gets each element in molecule
        atom_positions = mol.GetConformer().GetPositions() # coordinates of atoms


        atoms = Atoms(symbols=atoms_symbols, positions=atom_positions) #creates atoms object with given coordinates
                                                                       #and elemental symbols, # of position arrays is equal
                                                                       #and corresponds to # of element symbols

        write(xyz_folder + "/" + nameOfJob +"_Molecule_" + str(i)+ ".xyz", atoms, format ="xyz") #writes .xyz file of molecule to desires folder

        print(f"Converted molecule {i+1} to XYZ format")

#path of SDF file data
sdf = '/Users/zaansaeed/Desktop/OMO Research/9:20_SDF_to_XYZ/SDF_Data_LinearPeptides/64peptideslinearstructures.sdf'
#path of folder where to put the XYZ files
xyz = '/Users/zaansaeed/Desktop/OMO Research/9:20_SDF_to_XYZ'

sdf_to_xyz("LinearPeptide_to_XYZ",sdf,xyz)