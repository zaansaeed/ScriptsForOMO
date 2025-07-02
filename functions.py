import csv
import os
import subprocess
import time
from tokenize import group

from rdkit.Chem import AllChem, FragmentOnBonds
from rdkit import Chem
import pandas as pd
from collections import deque
import math
import numpy as np
from natsort import natsorted



main_dir = os.path.abspath("/Users/zaansaeed/Desktop/NewPeptides")

def waiting_for_file(working_dir,name,wait_time) -> None:
    """
    Monitors a given directory for a file to exist and ensures the file size
    stabilizes before exiting the function. This function is useful when a file is
    being created or updated by another process, and you want to wait until it
    is fully written before proceeding.

    :param working_dir: The directory path to observe for the specified file.
    :type working_dir: str
    :param name: The name of the file to monitor within the directory.
    :type name: str
    :param wait_time: Time in seconds to wait between file size checks.
    :type wait_time: float
    :return: None
    """
    os.chdir(working_dir)
    prev_size = -1
    while not os.path.exists(name) or os.path.getsize(name) != prev_size:
        if os.path.exists(name):
            new_size = os.path.getsize(name)
            if new_size == prev_size:
                break
            prev_size = new_size
        print(f"waiting... for {name}")
        time.sleep(wait_time)


def split_xyz(working_dir,input_file,name,counter) -> int:
    """
    Splits a single XYZ file containing all conformation coordinates into
    separate XYZ files, each representing one conformation.

    :param working_dir: Directory where the output XYZ files will be stored
    :type working_dir: str
    :param input_file: Path to the input XYZ file containing all conformations
    :type input_file: str
    :param name: Base name for the generated XYZ output files
    :type name: str
    :return: None

    """
    temp_count = counter-1
    with open(input_file, "r") as f:
        lines = f.readlines()
    with open(input_file, "r") as f:
        num_coords = int(lines[0].split()[0])
        num_cords = num_coords + 2  # num of cords + 2 because of the length and name of peptide
        conformations = []  #will contain the n many arrays of coordinates of n conformations
        while True:
            lines = [f.readline() for i in range(num_cords)]
            if not any(lines):
                break
            conformations.append(lines)
    for i, conf in enumerate(conformations):
        output_file = os.path.join(working_dir, f"{name}_Conformation_{temp_count}.xyz")
        temp_count+=1
        with open(output_file, "w") as f:
            f.writelines(conf)
    return temp_count



def smile_to_mae(smile_string,name,working_dir) -> None:
    """
    Converts a SMILES string to Maestro (MAE) format and organizes the output
    into a specific folder structure. This function ensures that the required
    files (SMILES file and Maestro file) are present, generating them if
    necessary. It uses Schrodinger's LigPrep tool for the conversion process
    when the Maestro file is not already available.

    :param smile_string: The SMILES string representing the chemical structure
        that will be converted into a Maestro (.mae) file.
    :type smile_string: str
    :param name: The name used for creating a unique folder and identifying
        the associated files. This ensures proper organization and file
        management.
    :type name: str
    :return: None
    """

    def split_mae_structconvert(mae_file, structconvert_path):
        base_name = os.path.splitext(mae_file)[0]
        output_prefix = base_name + "_model.mae"

        subprocess.run([
            structconvert_path,
            "-split-nstructures", "1",
            mae_file,
            output_prefix
        ], check=True)

        print(f"Split files created with prefix {output_prefix}_*.mae")

    os.chdir(working_dir)
    if not os.path.exists(f"{name}.smi"): #if there is no smile file, create it -contains smile string
        temp_smi = os.path.join(working_dir, f"{name}.smi")
        with open(temp_smi, "w") as temp_file:
            temp_file.write(f"{smile_string}\n")
    if not os.path.exists(f"{name}.mae"): #if there is no .mae file, create it - ready for maestro input
        subprocess.run(["/opt/schrodinger/suites2024-3/ligprep", "-ismi", f"{name}.smi", "-omae", f"{name}.mae"])
        waiting_for_file(working_dir,f"{name}.mae",20) # wait until file is created

        split_mae_structconvert(f"{name}.mae","/opt/schrodinger/suites2024-3/utilities/structconvert") #split the .mae file into multiple files, one for each model
        #if tat exists


    print('smile ->mae')

    os.chdir(main_dir)



import glob
def get_split_files():
    files = glob.glob("*_model*.mae")
    filtered = [f for f in files if "out" not in f]
    return sorted([os.path.splitext(f)[0] for f in filtered])

def run_confSearch(name,working_dir) -> None:
    """
    Generates and executes a conformational search configuration for a specified molecular input file.

    This function sets up the required input file for a conformational search using
    Schrödinger's Macromodel and executes the conformational search. It monitors the
    creation and size of the output file to ensure completion and reports when done.

    :param name: The base name of the input .mae file and other related conformational
        search files (e.g., "molecule" for "molecule.mae").
    :param working_dir: The directory where the input and output files for conformational
        search are located or will be created.
    :return: None
    """
    os.chdir(working_dir)

    split_files = get_split_files()
    for file in split_files:
        if not os.path.exists(file):  # if the conf search file does not exist, we will create it.
            conf_search_file = os.path.join(working_dir, f"{file}")
            with open(conf_search_file, "w") as f:
                f.write(f"INPUT_STRUCTURE_FILE {file}.mae\n")
                f.write("JOB_TYPE CONFSEARCH\n")
                f.write("CONFSEARCH_METHOD MCMM\n")
                f.write("FORCE_FIELD OPLS_2005\n")
                f.write("SOLVENT None\n")
                f.write("DIELECTRIC_CONSTANT 1.0\n")
                f.write("CHARGES_FROM Force field\n")
                f.write("CUTOFF None\n")
                f.write("MINI_METHOD PRCG\n")
                f.write("MAXIMUM_ITERATION 2500\n")
                f.write("CONVERGE_ON Gradient\n")
                f.write("CONVERGENCE_THRESHOLD 0.05\n")
                f.write("OUTCONFS_PER_SEARCH 10000\n")
                f.write("CONFSEARCH_STEPS 1000\n")
                f.write("CONFSEARCH_STEPS_PER_ROTATABLE 100\n")
                f.write("ENERGY_WINDOW 104.6\n")
                f.write("CONFSEARCH_TORSION_SAMPLING Intermediate\n")
        if not os.path.exists(f"{file}-out.mae"):
            subprocess.run([
                "/opt/schrodinger/suites2024-3/macromodel",
                f"{file}"])
            waiting_for_file(working_dir, f"{file}-out.mae", 120)

    print("DONE RUNNING CONFSEARCH")

    os.chdir(main_dir)


def mae_to_pdb(name,working_dir) -> None:
    """
    Converts a Maestro (.mae) file to a Protein Data Bank (.pdb) file using the
    Schrödinger structconvert utility. The function navigates to the specified
    working directory, performs the file conversion, and ensures the output file
    is created before returning to the original directory.

    :param name: The name of the file without extension. This will be appended
        with '-out.mae' and '-out.pdb' for input and output file names,
        respectively.
    :type name: str
    :param working_dir: The working directory where the input file is located.
        It is also used for temporary operations during the file conversion.
    :type working_dir: str
    :return: None
    """
    #convert to pdb

    os.chdir(working_dir)

    split_files = get_split_files()
    for file in split_files:
        if not os.path.exists(f"{file}-out.pdb"):
            subprocess.run(["/opt/schrodinger/suites2024-3/utilities/structconvert", f"{file}-out.mae", f"{file}-out.pdb"])
            waiting_for_file(working_dir, f"{file}-out.pdb", 30)

    print("DONE CONVERTING MAE TO PDB")

    os.chdir(main_dir)


def pdb_to_xyz(name,working_dir) -> None:
    """
    Converts a PDB file to an XYZ file using Open Babel. This function changes the current
    working directory to the specified `working_dir` and triggers the conversion process.
    If the conversion output file is not found or its size changes during the process,
    it revalidates periodically until the file has stabilized or exists.

    :param name: Name of the file (without extension) to be converted.
    :type name: str
    :param working_dir: Directory where the PDB input file exists and where the XYZ
        output file will be created.
    :type working_dir: str
    :returns: None
    """
    #convert to xyz
    os.chdir(working_dir)
    split_files = get_split_files()
    for f in split_files:
        if not os.path.exists(f"{f}-out.xyz"):
            os.system(f"obabel -ipdb {f}-out.pdb -O {f}-out.xyz")

            waiting_for_file(working_dir,f"{f}-out.xyz",15)


    print("DONE CONVERTING PDB TO XYZ")

    os.chdir(main_dir)


def xyz_to_individual_xyz(og_name,working_dir) -> None:
    """
    Converts a combined XYZ file to individual XYZ files and saves them in a newly
    created directory named after the provided name followed by "_Conformations".
    If the target directory already exists, the function will skip the directory
    creation step. The process switches directories to perform the operation and
    returns to the initial working directory at the end.

    :param name: The base name used to create the target directory and individual XYZ
        files.
    :type name: str
    :param working_dir: The current working directory where the operation begins
        and where the combined XYZ file is expected to be located.
    :type working_dir: str
    :return: None
    """
    os.chdir(working_dir)
    split_files = get_split_files()
    conf_counter=1
    if not os.path.exists(f"{og_name}_Conformations"):  # create folder with all conforamtions
        os.mkdir(f"{og_name}_Conformations")
        temp_working_dir = working_dir + f"/{og_name}_Conformations"
        for f in split_files:
            os.chdir(temp_working_dir)
            num_confs = split_xyz(temp_working_dir,working_dir+f"/{f}-out.xyz",og_name,conf_counter)
            conf_counter+=num_confs
        os.chdir(working_dir)
        print("DONE CONVERTING XYZ TO INDIVIDUAL XYZ")

    os.chdir(main_dir)


def extract_energies_to_csv(og_name,working_dir) -> None:
    """
    Extracts energy values from a log file and saves them into a CSV file. Processes
    each line in a specified log file to identify and extract conformation energy
    values. These values are then arranged into a pandas DataFrame object, and
    stored in a CSV file with a specified naming convention.

    This function assumes that the log file includes lines containing the text
    "Conformation " and that the energy values are located in the fourth column
    after splitting the line.

    If a CSV file with the expected name already exists in the working directory,
    the function will not process the log file and terminates its operations without
    duplicating the file.

    :param name: The base name of the log file (without the file extension). The same
                 name is used as the base name for the output CSV file.
    :type name: str
    :param working_dir: Absolute or relative path to the working directory where the
                        input log file is located and the output CSV file will be saved.
    :type working_dir: str
    :return: This function does not return any value.
    :rtype: None
    """
    os.chdir(working_dir)
    names = get_split_files()
    all_energies = []
    conf_counter = 1

    for name in names:
        log_file = f"{name}.log"
        print(log_file)

        # Extract energies directly from .log file
        with open(log_file, "r") as f:
            for line in f:
                if "Conformation " in line:
                    try:
                        energy = float(line.split()[3])  # adjust if different column
                        all_energies.append((f"{og_name}_Conformation_{conf_counter}", energy))
                        conf_counter += 1
                    except ValueError:
                        continue

    # Sort by energy
    all_energies.sort(key=lambda x: x[1])

    # Save to total_energies.csv
    df = pd.DataFrame(all_energies, columns=["Name", "Energy"])
    df.to_csv(f"{og_name}_total_energies.csv", index=False)


    os.chdir(main_dir)





#+====================================================================================

def bfs_traversal(mol, starting_id) -> list[int]:
    """
    Performs a breadth-first search (BFS) traversal on a molecular graph starting from
    the specified atom ID (N-terminus) and returns a list containing the IDs of nitrogen atoms visited
    in the order they are encountered during the traversal.

    The function uses a queue-based BFS algorithm to explore the molecular structure,
    and only detects and collects nitrogen atoms ('N') encountered during the traversal.

    :param mol: A molecular object representing the molecular graph. Typically requires
        the object to support accessing atom information through methods like
        `GetAtomWithIdx` and `GetNeighbors`.
    :type mol: RDKit Mol

    :param starting_id: The index of the starting atom in the traversal. This is the
        point from which the BFS exploration begins.
    :type starting_id: int

    :return: A list of indices representing the nitrogen atoms encountered in BFS order
        during the traversal.
    :rtype: list[int]
    """
    visited = set()
    queue = deque([starting_id])
    bfs_order = []
    while queue:
        atom_id = queue.popleft()
        if atom_id not in visited:
            visited.add(atom_id)
            bfs_order.append(atom_id)
        atom = mol.GetAtomWithIdx(atom_id)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in visited:
                queue.append(neighbor.GetIdx())

    bfs_nitrogens = []
    for id in bfs_order:
        if mol.GetAtomWithIdx(id).GetSymbol() == 'N':
            bfs_nitrogens.append(id)
    return bfs_nitrogens

class AmideGroup:
    """
    Represents an amide group within a peptide, including its atomic components and relationships.

    This class identifies and stores the relevant atom IDs (Carbon, Oxygen, Nitrogen, and Hydrogen)
    of an amide group from a peptide structure. The primary purpose of this class is to facilitate
    access to these atom IDs, allowing further analysis or manipulation of the amide group within
    the peptide.
    """

    def __init__(self,IDS, group_num, peptide,amide_groups):
        # Recheck if this is the best way to find N and other IDs, based on how amide groups are found
        self.group_num = group_num
        self.N = IDS[0]

        nitrogen = peptide.GetAtomWithIdx(self.N)
        neighbors = nitrogen.GetNeighbors()
        hydrogen_id = None

        for neighbor in neighbors:
            if neighbor.GetSymbol() == 'H':
                hydrogen_id = neighbor.GetIdx()
            if neighbor.GetSymbol() == 'C':
                for neighbor2 in neighbor.GetNeighbors():
                    if neighbor2.GetSymbol() == 'O':
                        self.C=neighbor.GetIdx()
                        self.O=neighbor2.GetIdx()

        self.Residue1= None
        self.Residue2=None #only the last amide group will have 2 residues
        self.H = hydrogen_id
        self.atom_IDs = (self.C,self.O,self.N)
        #if group num is 1, then cut bond on Nitrogen, thats residue,
        # cut bond from nitrogen to oxygen and nitrogen to carbon, (second to last id )
        # ('NCC(=O)N'))# normal case, 2 carbons
        #     'NCCC(=O)N'))# abnormal case, 3 carbons
        atom1, atom2, atom3, atom4 = 0,0,0,0
        bond1,bond2 = None,None
        if group_num == 1: #if its the first amide, aka N-termini case
            atom1 = self.N #last nitrogen, which is where we want to cut
            atom2 = self.C #carbon to which this nitrogen is connected to
            for bond in peptide.GetBonds():
                if {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()} == {atom1, atom2}:
                    bond1 = bond.GetIdx()
            mol_frag = FragmentOnBonds(peptide,[bond1],addDummies=False)

            frags_ids = Chem.GetMolFrags(mol_frag,asMols=False)

            frags = []
            for atom_ids in frags_ids:
                # PathToSubmol keeps atom indices + conformers
                frag_mol = Chem.PathToSubmol(mol_frag, atom_ids, useQuery=False)
                if self.C in atom_ids:
                    self.Residue1 = (frag_mol,atom_ids)
                    break
        else: # if this is another amide, not the first one
            prev_amide = amide_groups[group_num-2]
            atom1 = prev_amide.getN() #the residue AFTER these atoms
            atom2 = prev_amide.getC()
            atom3 = self.N #the residue BEFORE these atoms (remember, the amide group youre currently on is  the 0th index of atom_ids, or self.N)
            atom4 = self.C #
            for bond in peptide.GetBonds():
                if {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()} == {atom1, atom2}:
                    bond1 = bond.GetIdx()
                if {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()} == {atom3, atom4}:
                    bond2 = bond.GetIdx()
            if bond1 == bond2: #trivial case, should never happen
                self.Residue1 = None
                return
            mol_frag = FragmentOnBonds(peptide,[bond1,bond2],addDummies=False)
            frags_ids = Chem.GetMolFrags(mol_frag,asMols=False)
            frags = []
            for atom_ids in frags_ids:
                # PathToSubmol keeps atom indices + conformers

                frag_mol = Chem.PathToSubmol(mol_frag, atom_ids, useQuery=False)
                if self.C in atom_ids:
                    self.Residue1 = (frag_mol, atom_ids)
                    break

        if self.group_num == 5:  # check if its the last one after residue1 is estavblished
            atom1 = self.N
            atom2 = self.C
            atom3 = IDS[-1]
            atom4 = IDS[-3]
            for bond in peptide.GetBonds():
                if {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()} == {atom1, atom2}:
                    bond1 = bond.GetIdx()
                if {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()} == {atom3, atom4}:
                    bond2 = bond.GetIdx()
            if bond1 == bond2:  # trivial case, should never happen
                self.Residue2 = None
                return
            mol_frag = FragmentOnBonds(peptide, [bond1, bond2], addDummies=False)
            frags_ids = Chem.GetMolFrags(mol_frag, asMols=False)
            frags = []
            for atom_ids in frags_ids:
                # PathToSubmol keeps atom indices + conformers
                frag_mol = Chem.PathToSubmol(mol_frag, atom_ids, useQuery=False)
                if self.N in atom_ids:
                    self.Residue2 = (frag_mol, atom_ids)
                    break






    def getResidue1(self):
        return self.Residue1
    def getResidue2(self):
        return self.Residue2
    def getIDs(self):
        return self.atom_IDs
    def getC(self):
        return self.C
    def getO(self):
        return self.O
    def getN(self):
        return self.N
    def getH(self):
        return self.H


def find_n_terminus(input_peptide) -> int:
    """
    Identifies and returns the index of the N-terminus in a given peptide structure. The N-terminus
    residue can exist in either a normal or an abnormal form. The function searches the peptide
    structure for these substructural matches and determines the N-terminus accordingly.

    :param input_peptide: The molecular representation of the peptide structure to be analyzed.
                         This is an RDKit Mol object that must represent a peptide structure.
    :return: The zero-based index of the atom constituting the N-terminus of the peptide.
    :rtype: int
    :raises Exception: If no N-terminus structure is found in the input peptide.
    """
    n_terminus_residue_normal = input_peptide.GetSubstructMatches(Chem.MolFromSmarts('[NH2]C[C](=O)'))
    n_terminus_residue_abnormal = input_peptide.GetSubstructMatches(Chem.MolFromSmarts('[NH2]CC[C](=O)'))
    try:
        if n_terminus_residue_normal:
            return n_terminus_residue_normal[0][0]
        else:
            return n_terminus_residue_abnormal[0][0]
    except IndexError:
        raise Exception("No N-terminus found in the input peptide.")

def add_amides(input_peptide) -> list[AmideGroup]:
    """
    Identifies amide groups in the given peptide molecule and returns them as a list of AmideGroup
    objects. This function searches for specific substructure patterns indicative of amide groups.
    The function also ensures that only unique amide groups, identified by their atom IDs, are
    added to the resulting list. It processes the molecule in a breadth-first traversal order
    starting from the N-terminus.

    :param input_peptide: The input molecule, represented as an RDKit molecule object, in which
                          amide groups are to be identified.
    :type input_peptide: Chem.Mol
    :return: A list of AmideGroup objects representing the identified amide groups.
    :rtype: list[AmideGroup]
    """
    amide_groups, used_ids = [],[]
    matches1 = input_peptide.GetSubstructMatches(Chem.MolFromSmiles('NCC(=O)N'))# normal case, 2 carbons
    matches2 = input_peptide.GetSubstructMatches(Chem.MolFromSmiles('NCCC(=O)N'))# abnormal case, 3 carbons
    matches = matches1 + matches2

    n_terminus = find_n_terminus(input_peptide)
    bfs_order = bfs_traversal(input_peptide, n_terminus)

    i = 1
    for nID in bfs_order:
        for match in matches:
            if nID in match and n_terminus not in match:
                amide = AmideGroup(match, i, input_peptide,amide_groups)
                amide_ids = amide.getIDs()
                if amide_ids not in used_ids:

                    amide_groups.append(amide)
                    used_ids.append(amide_ids)
                    i += 1
    return amide_groups


def get_amide_distances(amide_groups,peptide) -> list[list[float]]:
    """
    Computes the distances between amide groups in a protein structure based on their
    atomic coordinates. This function calculates a pairwise distance matrix where each
    entry represents the distance between the Hydrogen atom of one amide group's Nitrogen
    and the Oxygen atom of another amide group's Carbon. For amide groups without a Nitrogen
    Hydrogen atom (e.g., proline), the distance is set to -1.0.

    :param amide_groups: A list of amide group instances, where each instance should have
        methods `getH()` to return the Hydrogen atom index (or None if not present) and
        `getO()` to return the Oxygen atom index.
    :type amide_groups: list
    :param atom_coordinates: A dictionary mapping atomic indices to their 3D coordinates as
        numpy arrays. The keys correspond to atomic indices, and the values are numpy arrays
        of shape (3,) representing x, y, z coordinates.
    :return: A 2D distance matrix where the entry at [i][j] represents the distance between
        the Hydrogen atom of the ith amide group and the Oxygen atom of the jth amide group.
        Entries for self-distances (i.e., when i == j) are set to 0.0 if a Hydrogen atom exists
        for the amide group. For any missing Hydrogen, the distance is -1.0 (e.g., prolines).
    :rtype: list[list[float]]
    """

    distance_matrix = [[0.0 for _ in range(len(amide_groups))] for _ in range(len(amide_groups))]
    for i,amid1 in enumerate(amide_groups):
        for j,amid2 in enumerate(amide_groups):
            if i == j and amid1.getH() is not None:
                distance_matrix[i][j] = 0.0
            else:
                amid1H = amid1.getH()
                if amid1H is None:
                    distance_matrix[i][j] = -1.0 #prolines have no nitrogen Hydrogen, so the distance is non existant
                else:
                    amid2O = amid2.getO()
                    amid1H_pos = peptide.GetConformer().GetAtomPosition(amid1H)
                    amid2O_pos = peptide.GetConformer().GetAtomPosition(amid2O)
                    distance_matrix[i][j] = np.linalg.norm(np.array(amid1H_pos)-np.array(amid2O_pos))
    return distance_matrix


##ENSURE THIS WORKS
def xyz_to_array(xyz_file) -> list[list[float]]:
    """
    Convert the content of an `.xyz` file into a list representation where the
    coordinates of atoms are stored along with the number of atoms in the file.
    The `.xyz` format is commonly used for molecular structure files containing
    atomic positions.

    This function reads the file, extracts the number of atoms, and iterates through
    the lines containing atomic coordinate data, converting them into a list of
    float values.

    :param xyz_file: Path to the input `.xyz` file containing atomic coordinate data.
    :type xyz_file: str
    :return: A list where the first element is an integer indicating the number
        of atoms, and subsequent elements are sublists, each containing the
        x, y, and z coordinates of one atom.
    :rtype: list[list[float]]
    """
    #oth index of coordinates is num atoms
    #ith index of coordinates is ith atom ID i.e. index 1 is the first atom ID 1

    with open(xyz_file, 'r', encoding= 'utf-8') as file:
        lines = file.readlines()
    num_atoms = int(lines[0].strip())
    coordinates = []
    coordinates.append(num_atoms)

    for line in lines[2:]:
        parts = line.split()
        x, y, z = map(float, parts[1:4])
        coordinates.append([x, y, z])
    return coordinates


def boltzmann_weight_energies(og_name,working_dir, update_matrices) -> None:
    """
    Calculate and store the Boltzmann-weighted energies based on molecular conformations.

    This function performs multiple operations including creating molecular
    representations, calculating amide group distances for conformations, computing
    Boltzmann-weighted matrices of distances, and saving results to a CSV file. The
    process involves reading input SMILES files, generating molecular conformations,
    and determining their Boltzmann-weighted contributions.

    :param name: The name identifier for the molecular system. This is used to
                 reference the input SMILES (.smi) file and naming the output
                 files for distances.
    :type name: str
    :param working_dir: The file path to the working directory where input files
                        are located and output files are saved. Subdirectories for
                        intermediate data may also be created under this path.
    :type working_dir: str
    :param update_matrices: A flag indicating whether to force the update of the
                            Boltzmann-weighted matrix and overwrite the existing
                            result file.
    :type update_matrices: bool
    :return: None
    """

    os.chdir(working_dir)
    names = get_split_files()

    if os.path.exists(f"{og_name}-BWdistances.csv"):
        with open(f"{og_name}.smi", "r") as f:
            smiles_string = f.readlines()[0]

        peptide = Chem.AddHs(Chem.MolFromSmiles(smiles_string))
        amideGroups = add_amides(peptide)

        temp_working_dir = working_dir + f"/{og_name}_Conformations"
        os.chdir(temp_working_dir) #working in conforamtions folder
        distances = []
        for conformation_xyz in natsorted(os.listdir(temp_working_dir)):
            peptide.RemoveAllConformers()
            peptide = load_xyz_coords(peptide,f"{conformation_xyz}")
            distances.append(get_amide_distances(amideGroups,peptide))
        os.chdir(working_dir)

        boltzmann_matrix = boltzmann_weighted_average(distances, working_dir,og_name)


        df = pd.DataFrame(boltzmann_matrix)
        df.to_csv(working_dir+f'/{og_name}-BWdistances.csv', index=False, header=False)



    os.chdir(main_dir)




################################################
################################################

def load_xyz_coords(mol, xyz_path) -> Chem.Mol:
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

    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol

def calculate_dihedrals(residue,mol):
    def get_dihedral_angle(p1, p2, p3, p4):
        """Calculate the dihedral angle between four points in 3D space."""
        # Convert to numpy arrays
        p1, p2, p3, p4 = map(np.array, (p1, p2, p3, p4))

        # Bond vectors
        b1 = p2 - p1
        b2 = p3 - p2
        b3 = p4 - p3

        # Normal vectors
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)

        # Normalize normals
        n1 /= np.linalg.norm(n1)
        n2 /= np.linalg.norm(n2)

        # Normalize b2
        b2_unit = b2 / np.linalg.norm(b2)

        # Compute the angle
        x = np.dot(n1, n2)
        y = np.dot(np.cross(n1, n2), b2_unit)

        angle = np.arctan2(y, x)
        return np.degrees(angle)



    if len(residue) == 5: #n-terminus case (normal) (5000,5000,psi) [NH2]C[C](=O)[N]
        temp_dihedrals = []
        temp_dihedrals.append(5000)
        temp_dihedrals.append(5000)
        p1= mol.GetConformer().GetAtomPosition(residue[0])
        p2 = mol.GetConformer().GetAtomPosition(residue[1])
        p3= mol.GetConformer().GetAtomPosition(residue[2])
        p4 = mol.GetConformer().GetAtomPosition(residue[4])
        temp_dihedrals.append(get_dihedral_angle(p1,p2,p3,p4))
        return temp_dihedrals

    if len(residue) == 6: #n-terminus case (abnormal) (5000,theta,psi) [NH2]CC[C](=O)[N]
        temp_dihedrals = []
        temp_dihedrals.append(5000)
        p1= mol.GetConformer().GetAtomPosition(residue[0])
        p2= mol.GetConformer().GetAtomPosition(residue[1])
        p3=mol.GetConformer().GetAtomPosition(residue[2])
        p4 = mol.GetConformer().GetAtomPosition(residue[3])
        temp_dihedrals.append(get_dihedral_angle(p1,p2,p3,p4))
        p1= mol.GetConformer().GetAtomPosition(residue[1])
        p2= mol.GetConformer().GetAtomPosition(residue[2])
        p3=mol.GetConformer().GetAtomPosition(residue[3])
        p4 = mol.GetConformer().GetAtomPosition(residue[5])

        temp_dihedrals.append(get_dihedral_angle(p1,p2,p3,p4))
        return temp_dihedrals
    if len(residue) == 7: #normal   (phi,5000,psi) C(=O)NCC(=O)N
        temp_dihedrals = []
        p1= mol.GetConformer().GetAtomPosition(residue[0])
        p2 = mol.GetConformer().GetAtomPosition(residue[2])
        p3 = mol.GetConformer().GetAtomPosition(residue[3])
        p4 = mol.GetConformer().GetAtomPosition(residue[4])
        temp_dihedrals.append(get_dihedral_angle(p1,p2,p3,p4))
        temp_dihedrals.append(5000)
        p1= mol.GetConformer().GetAtomPosition(residue[2])
        p2 = mol.GetConformer().GetAtomPosition(residue[3])
        p3 = mol.GetConformer().GetAtomPosition(residue[4])
        p4 = mol.GetConformer().GetAtomPosition(residue[6])
        temp_dihedrals.append(get_dihedral_angle(p1,p2,p3,p4))
        return temp_dihedrals

    if len(residue)==8: #abnormal case (phi,theta,psi) C(=O)NCCC(=O)N
        temp_dihedrals = []
        p1= mol.GetConformer().GetAtomPosition(residue[0])
        p2 = mol.GetConformer().GetAtomPosition(residue[2])
        p3 = mol.GetConformer().GetAtomPosition(residue[3])
        p4 = mol.GetConformer().GetAtomPosition(residue[4])
        temp_dihedrals.append(get_dihedral_angle(p1,p2,p3,p4))
        p1= mol.GetConformer().GetAtomPosition(residue[2])
        p2 = mol.GetConformer().GetAtomPosition(residue[3])
        p3 = mol.GetConformer().GetAtomPosition(residue[4])
        p4 = mol.GetConformer().GetAtomPosition(residue[5])
        temp_dihedrals.append(get_dihedral_angle(p1,p2,p3,p4))
        p1= mol.GetConformer().GetAtomPosition(residue[3])
        p2 = mol.GetConformer().GetAtomPosition(residue[4])
        p3 = mol.GetConformer().GetAtomPosition(residue[5])
        p4 = mol.GetConformer().GetAtomPosition(residue[7])
        temp_dihedrals.append(get_dihedral_angle(p1,p2,p3,p4))
        return temp_dihedrals
    return



def extract_boltzmann_weighted_dihedrals_normalized(main_dir):
    os.chdir(main_dir)
    for folder in natsorted(os.listdir(main_dir)):
        if os.path.isdir(folder):
            os.chdir(folder)
            working_dir = os.getcwd()
            name = folder.split("_")[1]
            print(name)
            if not os.path.exists(f"{name}-BWDihedralNormalized.csv"):
                smiles_string = open(f"{name}.smi").read().strip() #generate the smiles string, currently working in Peptide _XXXX folder
                peptide_normalized_dihedrals = []
                mol = Chem.MolFromSmiles(smiles_string)
                mol = Chem.AddHs(mol)
                #add embeding?
                n_terminus = find_n_terminus(mol)
                nitrogen_order = bfs_traversal(mol, n_terminus)
                n_terminus_residue_normal = mol.GetSubstructMatches(Chem.MolFromSmarts('[NH2]C[C](=O)[N]'))
                n_terminus_residue_abnormal = mol.GetSubstructMatches(Chem.MolFromSmarts('[NH2]CC[C](=O)[N]'))
                normal_residues = mol.GetSubstructMatches(Chem.MolFromSmiles('C(=O)NCC(=O)N'))
                abnormal_residues = mol.GetSubstructMatches(Chem.MolFromSmiles('C(=O)NCCC(=O)N'))
                all_residues = normal_residues + abnormal_residues + n_terminus_residue_normal + n_terminus_residue_abnormal
                all_residues = list(all_residues)
                ############## get rid of asparagine
                query1 = Chem.MolFromSmarts('C(=O)[N]C[C](=O)[NH2]')
                query2 = Chem.MolFromSmarts('C(=O)[N]CC[C](=O)[NH2]')
                query3 = Chem.MolFromSmarts('[NH2]C[C](=O)[NH2]')
                query4 = Chem.MolFromSmarts('[NH2]CC[C](=O)[NH2]')
                asparagines = mol.GetSubstructMatches(query1) + mol.GetSubstructMatches(
                    query2) + mol.GetSubstructMatches(query3) + mol.GetSubstructMatches(query4)
                for asparagine in asparagines:
                    for residue in all_residues:
                        if asparagine[-1] in residue:
                            all_residues.remove(residue)
                ##########################
                ordered_residues = []
                for id in nitrogen_order:
                    for residue in all_residues:
                        if id in residue:
                            ordered_residues.append(residue)
                            all_residues.remove(residue)

                for conformation_xyz in natsorted(os.listdir(f"{name}_Conformations")):
                    if conformation_xyz.endswith('.xyz'): #working within 1 conformer
                        mol.RemoveAllConformers()
                        mol = load_xyz_coords(mol, f"{working_dir}/{name}_Conformations/{conformation_xyz}")
                        conformation_dihedrals = [] #contains (phi,theta,psi) 6 times, 1 for each residue
                        for residue in ordered_residues:
                            conformation_dihedrals.append(calculate_dihedrals(residue,mol))
                        #convert each angle to sin/cos components, and add flag = 1 if row contains 5000
                        flag = (1,0)
                        conformation_normalized_dihedrals = []
                        for i in range(len(conformation_dihedrals)): #num of residues
                            residue_normalized_dihedrals_and_flag = []
                            flag = (1,0)
                            for angle in conformation_dihedrals[i]: #working on converting (phi,theta,psi) -> ((sin,cos)...(sin,cos), flag)
                                if angle > 1000:
                                    flag = (0,0)
                                    residue_normalized_dihedrals_and_flag.append((0,0))
                                else:
                                    residue_normalized_dihedrals_and_flag.append((math.sin(math.radians(angle)),math.cos(math.radians(angle))))
                            residue_normalized_dihedrals_and_flag.append(flag)

                            conformation_normalized_dihedrals.append(residue_normalized_dihedrals_and_flag)

                    peptide_normalized_dihedrals.append(conformation_normalized_dihedrals)



                #boltzmann weight the n many conformation [(sin,cos)...(sin,cos),flag]
                boltzmann_matrix = boltzmann_weighted_average(peptide_normalized_dihedrals, working_dir,name)


                boltzmann_matrix = boltzmann_matrix.reshape(len(boltzmann_matrix),-1)
                print(boltzmann_matrix)
                df = pd.DataFrame(boltzmann_matrix)
                df.to_csv(working_dir+f'/{name}-BWDihedralNormalized.csv', index=False, header=False)
                print("dihedral calculation done for " +name )
                os.chdir(main_dir)

            os.chdir(main_dir)



def boltzmann_weighted_average(values, working_dir, name):
    # Load total energies
    pattern = os.path.join(working_dir, f"*{name}*energies*.csv")
    matching_files = glob.glob(pattern)
    filepath = matching_files[0]
    energy_df = pd.read_csv(filepath)

    first_col = energy_df.columns[0]  # Get the name of the first column dynamically
    second_col = energy_df.columns[1]

    energy_dict = dict(zip(energy_df[first_col], energy_df[second_col]))

    energy_dict = natsorted(energy_dict.items(), key=lambda x: x[0])
    energy_vals = np.array([x[1] for x in energy_dict])
    # Boltzmann weighting
    R = 8.314e-3  # kJ/mol·K
    T = 298       # K
    weights = np.exp(-energy_vals / (R * T))
    weights /= np.sum(weights)

    # Weighted sum
    weighted_sum = np.zeros_like(values[0], dtype=float)
    for i, arr in enumerate(values):
        weighted_sum += weights[i] * np.array(arr)

    return weighted_sum













