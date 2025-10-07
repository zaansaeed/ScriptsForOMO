import os
import subprocess
import time
from rdkit.Chem import FragmentOnBonds
from rdkit import Chem
from rdkit.Chem import Descriptors3D, Descriptors
import pandas as pd
from collections import deque
import math
import numpy as np
import random
from natsort import natsorted
import glob
import logging
from collections import defaultdict, Counter



logger = logging.getLogger("data_logger")
config = None
schrodinger_path = None
struct_convert_path = None
macro_model_path = None
lig_prep_path = None
main_dir = None
number_residues = None

def init_config(cfg) -> dict:
    global config, schrodinger_path, macro_model_path, struct_convert_path, lig_prep_path, main_dir, number_residues
    config = cfg
    schrodinger_path = config["data_generation"]["schrodinger_path"]
    struct_convert_path = schrodinger_path + "/utilities/structconvert"
    macro_model_path = schrodinger_path + "/macromodel"
    lig_prep_path = schrodinger_path + "/ligprep"
    main_dir = config["data_generation"]["main_dir"]
    number_residues = config["data_generation"]["number_of_residues"]
    return config



def create_target_file(name,value,working_dir,target_file_name) -> None:
    os.chdir(working_dir)
    if not os.path.exists(f"{name}_{target_file_name}.txt"):
        with open(f"{name}_{target_file_name}.txt", "w") as f:
            f.write(f"{value}")
        logger.info(f"Created {name}_{target_file_name}.txt")


def load_lines(filepath) -> list[str]:
    with open(filepath, "r") as f:
        return [line.strip() for line in f]


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
        logger.info(f"Waiting... for {name}")
        time.sleep(wait_time)


def split_xyz(working_dir,input_file,name,counter) -> int:
    """
    Takes an xyz with multple molecules and splits each into its own xyz file in working_dir ( a folder called {name} conformations. 
    Returns the number of molecules in the xyz file for indexing purposes. Assumes that xyz file has the following format (typical format): 
    
    #  numAtoms
    # Name
    # Atom_1 X Y Z
    # etc..
    # Atom_numAtoms X Y Z

    :param working_dir: 
    :param input_file: 
    :param name: 
    :param counter: 
    :return: int
    """"""
    """
    temp_count = counter-1
    with open(input_file, "r") as f:
        lines = f.readlines()
    with open(input_file, "r") as f:
        num_coords = int(lines[0].split()[0])
        num_cords = num_coords + 2  # num of cords + 2 because of the length and name of peptide
        conformations = []  #will contain the n many arrays of coordinates of n conformations
        while True:
            lines = [f.readline() for _ in range(num_cords)]
            if not any(lines):
                break
            conformations.append(lines)
    for _, conf in enumerate(conformations):
        output_file = os.path.join(working_dir, f"{name}_Conformation_{temp_count}.xyz")
        temp_count+=1
        with open(output_file, "w") as f:
            f.writelines(conf)
    return temp_count


def split_mae_structconvert(mae_file, structconvert_path) -> None:
    """
    Splits a .mae file with multiple structures into multiple .mae files with 1 structure each. Appends "_model" to the end of the file name.
    :param mae_file:
    :param structconvert_path:
    :return:
    """

    base_name = os.path.splitext(mae_file)[0]
    output_prefix = base_name + "_model.mae"
    try:
        subprocess.run([
            structconvert_path,
            "-split-nstructures", "1",
            mae_file,
            output_prefix
        ], check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"Error running structconvert: {e}")
        raise RuntimeError(f"Failed to split MAE file {mae_file}. Please check the input file or the structconvert installation.")

    logger.info(f"Split files created with prefix {output_prefix}_*.mae")


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

    :param working_dir: The working directory where the output files will be
    :return: None
    :param smile_string:
    :param name:
    :param working_dir:
    :return:
    """

    os.chdir(working_dir)
    if not os.path.exists(f"{name}.smi"): #if there is no smile file, create it - contains only the smile string
        temp_smi = os.path.join(working_dir, f"{name}.smi")
        with open(temp_smi, "w") as temp_file:
            temp_file.write(f"{smile_string}\n")

    if not os.path.exists(f"{name}.mae"): #if there is no .mae file, create it - ready for maestro input
        try:
            subprocess.run([lig_prep_path, "-ismi", f"{name}.smi", "-omae", f"{name}.mae"],check=True)
            waiting_for_file(working_dir,f"{name}.mae",20) # wait until file is created
            split_mae_structconvert(f"{name}.mae",struct_convert_path) #split the .mae file into multiple files, one for each model, if a tautomer exists
            logger.info(f"Created {name}.mae")

        except subprocess.CalledProcessError as e:
            logger.error(f"Error running LigPrep: {e}")
            raise RuntimeError(f"Failed to convert SMILES to MAE for {name}. Please check the input SMILES string or the LigPrep installation.")
        


def get_split_files() -> list[str]:
    """
    Returns a list of all the names that end with "_model" in the current working directory.
    :return:
    """
    files = glob.glob("*_model*.mae")
    filtered = [f for f in files if "out" not in f]
    return sorted([os.path.splitext(f)[0] for f in filtered])


def run_confSearch(working_dir,wait_time) -> None:
    """
    Generates and executes a conformational search configuration for a specified molecular input file.

    This function sets up the required input file for a conformational search using
    Schrödinger's Macromodel and executes the conformational search. It monitors the
    creation and size of the output file to ensure completion and reports when done.

    :param working_dir: The directory where the input and output files for conformational
        search are located or will be created.
    :param wait_time: The time in seconds to wait for the conformational search to complete.
    :return: None
    """
    os.chdir(working_dir)

    names = get_split_files()  # get all the names of the files that end with "_model"
    for name in names:
        if not os.path.exists(name):  # if the conf search input file does not exist, we will create it.
            conf_search_file = os.path.join(working_dir, f"{name}")
            with open(conf_search_file, "w") as f:
                f.write(f"INPUT_STRUCTURE_FILE {name}.mae\n")
                f.write("JOB_TYPE CONFSEARCH\n")
                f.write(f"CONFSEARCH_METHOD {config["data_generation"]["conf_search_settings"]["CONFSEARCH_METHOD"]}\n")
                f.write(f"FORCE_FIELD {config["data_generation"]["conf_search_settings"]["FORCE_FIELD"]}\n")
                f.write(f"SOLVENT {config["data_generation"]["conf_search_settings"]["SOLVENT"]}\n")
                f.write(f"DIELECTRIC_CONSTANT {config["data_generation"]["conf_search_settings"]["DIELECTRIC_CONSTANT"]}\n")
                f.write(f"CHARGES_FROM {config["data_generation"]["conf_search_settings"]["CHARGES_FROM"]}\n")
                f.write(f"CUTOFF {config["data_generation"]["conf_search_settings"]["CUTOFF"]}\n")
                f.write(f"MINI_METHOD {config["data_generation"]["conf_search_settings"]["MINI_METHOD"]}\n")
                f.write(f"MAXIMUM_ITERATION {config["data_generation"]["conf_search_settings"]["MAXIMUM_ITERATIONS"]}\n")
                f.write(f"CONVERGE_ON {config["data_generation"]["conf_search_settings"]["CONVERGE_ON"]}\n")
                f.write(f"CONVERGENCE_THRESHOLD {config["data_generation"]["conf_search_settings"]["CONVERGENCE_THRESHOLD"]}\n")
                f.write(f"OUTCONFS_PER_SEARCH {config["data_generation"]["conf_search_settings"]["OUTCONFS_PER_SEARCH"]}\n")
                f.write(f"CONFSEARCH_STEPS {config["data_generation"]["conf_search_settings"]["CONFSEARCH_STEPS"]}\n")
                f.write(f"CONFSEARCH_STEPS_PER_ROTATABLE {config["data_generation"]["conf_search_settings"]["CONFSEARCH_STEPS_PER_ROTATABLE"]}\n")
                f.write(f"ENERGY_WINDOW {config["data_generation"]["conf_search_settings"]["ENERGY_WINDOW"]}\n")
                f.write(f"CONFSEARCH_TORSION_SAMPLING {config["data_generation"]["conf_search_settings"]["CONFSEARCH_TORSION_SAMPLING"]}\n")

            logger.info(f"Created {name}, running confsearch now...")
        if not os.path.exists(f"{name}-out.mae"):
            try:
                subprocess.run([
                    macro_model_path,
                    f"{name}"],check=True)
                waiting_for_file(working_dir, f"{name}-out.mae", wait_time)
                logger.info(f"Conformational search complete for {name}")
            except subprocess.CalledProcessError as e:
                logger.error(f"Error running MacroModel: {e}")
                raise RuntimeError(f"Failed to run conformational search for {name}. Please check the input file or the MacroModel installation.")

def mae_to_pdb(working_dir) -> None:
    """
    Converts a Maestro (.mae) file to a Protein Data Bank (.pdb) file using the
    Schrödinger structconvert utility. The function navigates to the specified
    working directory, performs the file conversion, and ensures the output file
    is created before returning to the original directory.

    :param working_dir: The working directory where the input file is located.
        It is also used for temporary operations during the file conversion.
    :type working_dir: str
    :return: None
    """

    os.chdir(working_dir)

    split_files = get_split_files()

    for file in split_files:
        if not os.path.exists(f"{file}-out.pdb"):
            try:
                subprocess.run([struct_convert_path, f"{file}-out.mae", f"{file}-out.pdb"])
                waiting_for_file(working_dir, f"{file}-out.pdb", 30)
                logger.info(f"Created {file}-out.pdb from conformational search")
                if not os.path.exists(f"{file}-out-template.pdb"):
                    create_template_model(f"{file}-out.pdb",f"{file}-out-template.pdb")
                    logger.info(f"Created {file}-out-template.pdb to use for feature calculations")
            except subprocess.CalledProcessError as e:
                logger.error(f"Error running structconvert: {e}")
                raise RuntimeError(f"Failed to convert MAE to PDB for {file}. Please check the input file or the structconvert installation.")


def pdb_to_xyz(working_dir) -> None:
    """
    Converts a PDB file to an XYZ file using Open Babel. This function changes the current
    working directory to the specified `working_dir` and triggers the conversion process.
    If the conversion output file is not found or its size changes during the process,
    it revalidates periodically until the file has stabilized or exists.

    :param working_dir: Directory where the PDB input file exists and where the XYZ
        output file will be created.
    :type working_dir: str
    :returns: None
    """
    os.chdir(working_dir)
    split_files = get_split_files()
    for f in split_files:
        if not os.path.exists(f"{f}-out.xyz"):
            try:
                subprocess.run(["obabel", "-ipdb", f"{f}-out.pdb", "-O", f"{f}-out.xyz"],check=True)           
                waiting_for_file(working_dir,f"{f}-out.xyz",15)
                logger.info(f"Created {f}-out.xyz from {f}-out.pdb")
            except subprocess.CalledProcessError as e:
                logger.error(f"Error running Open Babel: {e}")
                raise RuntimeError(f"Failed to convert PDB to XYZ for {f}. Please check the input file or the Open Babel installation.")


def xyz_to_individual_xyz(og_name,working_dir) -> None:
    """
    Converts XYZ file into its individual structure XYZ files if they exist. Creates {og_name}_Conformations folder, and combines all _model XYZ files into 1 folder (if multiple conf searches for the same molecule occur due to tautomer)

    :param og_name: The base name used to create the target directory and individual XYZ
        files.
    :type og_name: str
    :param working_dir: The current working directory where the operation begins
        and where the combined XYZ file is expected to be located.
    :type working_dir: str
    :return: None
    """
    os.chdir(working_dir)
    split_files = get_split_files()
    conf_counter = 1
    if not os.path.exists(f"{og_name}_Conformations"):  # create folder with all conforamtions
        os.mkdir(f"{og_name}_Conformations")
        temp_working_dir = working_dir + f"/{og_name}_Conformations"
        for f in split_files:
            os.chdir(temp_working_dir)
            num_confs = split_xyz(temp_working_dir,working_dir+f"/{f}-out.xyz",f,conf_counter)
            conf_counter+=num_confs
        os.chdir(working_dir)
        logger.info(f"Created {og_name}_Conformations folder with all conformations")


def extract_energies_to_csv(og_name,working_dir) -> None:
    """
    Extracts all energies from all logs from conf searches of each model and saves it into 1 .csv file, sorted from lowest energy to highest.

    :param og_name: The base name of the log file (without the file extension). The same
                 name is used as the base name for the output CSV file.
    :type og_name: str
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

    if not os.path.exists(f"{og_name}_total_energies.csv"):
        for name in names:
            log_file = f"{name}.log"

            # Extract energies directly from .log file
            with open(log_file, "r") as f:
                for line in f:
                    if "Conformation " in line:
                        try:
                            energy = float(line.split()[3])  # adjust if different column
                            all_energies.append((f"{og_name}_Conformation_{conf_counter}", energy))
                            conf_counter += 1
                        except ValueError:
                            logger.info(f"ERROR: Could not convert energy value in line: {line.strip()}")

        # Sort by energy
        all_energies.sort(key=lambda x: x[1])

        # Save to total_energies.csv
        df = pd.DataFrame(all_energies, columns=["Name", "Energy"])
        df.to_csv(f"{og_name}_total_energies.csv", index=False)
        logger.info(f"Created {og_name}_total_energies.csv")


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
    for ID in bfs_order:
        if mol.GetAtomWithIdx(ID).GetSymbol() == 'N':
            bfs_nitrogens.append(ID)
    return bfs_nitrogens


class AmideGroup:
    """
    Represents an amide group within a peptide, including its atomic components and relationships.

    This class identifies and stores the relevant atom IDs (Carbon, Oxygen, Nitrogen, and Hydrogen)
    of an amide group from a peptide structure. The primary purpose of this class is to facilitate
    access to these atom IDs, allowing further analysis or manipulation of the amide group within
    the peptide.
    """

    def __init__(self,ids_of_match, group_num, peptide,amide_groups):

        self.group_num = group_num
        self.N = ids_of_match[0] #nitrogen is always the first atom in the match group, since substructures start with "N"
        self.C =None
        self.O = None
        #residue associate with amide group, behind it if behind is towards the n-terminus
        self.Residue1 = None
        self.Residue2 = None #only the last amide group will have two residues, so this will be none for n-1 amide groups if there are n residues

        nitrogen = peptide.GetAtomWithIdx(self.N)
        hydrogen_id = None
        for neighbor in nitrogen.GetNeighbors(): #this looks for the hydrogen atom of the nitrogen atom if it exists.
                                                # Also looks for the Carbon and Oxygen in the amide group (this should be behind the nitrogen, or more towards the n-terminus)
            if neighbor.GetSymbol() == 'H':
                hydrogen_id = neighbor.GetIdx()
            if neighbor.GetSymbol() == 'C':
                for neighbor2 in neighbor.GetNeighbors():
                    if neighbor2.GetSymbol() == 'O':
                        self.C=neighbor.GetIdx()
                        self.O=neighbor2.GetIdx()


        self.H = hydrogen_id
        self.atom_IDs = (self.C,self.O,self.N)
        # ('NCC(=O)N')) normal case, 2 carbons
        #'NCCC(=O)N')) abnormal case, 3 carbons

        bond1,bond2 = None,None

        if None not in self.atom_IDs:
            if group_num == 1: #if it's the first amide, aka N-termini residue, we want to cut the bond that connects this residue to the rest of the molecule and save those ids
                atom1 = self.N
                atom2 = self.C
                for bond in peptide.GetBonds():
                    if {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()} == {atom1, atom2}:
                        bond1 = bond.GetIdx()
                mol_frag = FragmentOnBonds(peptide,[bond1],addDummies=False)
                frags_ids = Chem.GetMolFrags(mol_frag,asMols=False)
                for atom_ids in frags_ids:
                    if self.C in atom_ids:
                        self.Residue1 = atom_ids
                        break
            else: # if this is another amide, not the first one
                # the residue AFTER (towards C-terminus)  these atoms
                prev_amide = amide_groups[group_num-2]
                atom1 = prev_amide.getN()
                atom2 = prev_amide.getC()
                # the residue BEFORE (towards the N-terminus) these atoms (remember, the amide group you're currently on is the 0th index of atom_ids, or self.N)
                atom3 = self.N
                atom4 = self.C
                for bond in peptide.GetBonds():
                    if {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()} == {atom1, atom2}:
                        bond1 = bond.GetIdx()
                    if {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()} == {atom3, atom4}:
                        bond2 = bond.GetIdx()
                if bond1 == bond2: #trivial case, should never happen. necessary since we check all substructure matches
                    self.Residue1 = None
                    return
                mol_frag = FragmentOnBonds(peptide,[bond1,bond2],addDummies=False)
                frags_ids = Chem.GetMolFrags(mol_frag,asMols=False)
                for atom_ids in frags_ids:
                    if self.C in atom_ids:
                        self.Residue1 = atom_ids
                        break
            # if this is the last amide group, we want to cut the bond that connects this residue to the rest of the molecule and save those ids
            # this is the case where we have a C-terminus residue, which is the last residue in the peptide chain.
            if self.group_num == 6-1 and not self.getIDs() in [amide.getIDs() for amide in amide_groups]:
                atom1 = self.N
                atom2 = self.C
                for bond in peptide.GetBonds():
                    if {bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()} == {atom1, atom2}:
                        bond1 = bond.GetIdx()
                mol_frag = FragmentOnBonds(peptide, [bond1], addDummies=False)
                frags_ids = Chem.GetMolFrags(mol_frag, asMols=False)
                for atom_ids in frags_ids:
                    if self.N in atom_ids:
                        self.Residue2 = atom_ids
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


def find_n_terminus(input_peptide) -> int | None:
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
    n_terminus_methylated = input_peptide.GetSubstructMatches(Chem.MolFromSmarts("[NH1]([CH3])[C]"))

    try:
        if n_terminus_residue_normal:
            return n_terminus_residue_normal[0][0]
        elif n_terminus_methylated:
            return n_terminus_methylated[0][0]
        else:
            return n_terminus_residue_abnormal[0][0]
    except IndexError:
        return None


def add_amides(input_peptide) -> list[AmideGroup] | None:
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
    matches3 = input_peptide.GetSubstructMatches(Chem.MolFromSmiles('NCC(=O)O')) # C terminus
    matches = matches1 + matches2 + matches3

    if not matches:
        raise Exception("No amide groups found in the input peptide.")  



    n_terminus = find_n_terminus(input_peptide)

    bfs_order = bfs_traversal(input_peptide, n_terminus)

    i = 1
    for nID in bfs_order:
        for match in matches:
            if nID in match and n_terminus not in match:
                amide = AmideGroup(match, i, input_peptide,amide_groups)
                amide_ids = amide.getIDs()
                if amide_ids not in used_ids and None not in amide_ids:
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
    :param peptide: The input molecule, represented as an RDKit molecule object. It has the position coordinates necessary to be processed.
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
                    distance_matrix[i][j] = -1.0 #If the nitrogen is methylated, then the distance doesn't exist
                else:
                    amid2O = amid2.getO()
                    amid1H_pos = peptide.GetConformer().GetAtomPosition(amid1H)
                    amid2O_pos = peptide.GetConformer().GetAtomPosition(amid2O)
                    distance_matrix[i][j] = np.linalg.norm(np.array(amid1H_pos)-np.array(amid2O_pos))
    return distance_matrix



def add_double_bonds_to_pdb(input_pdb) -> Chem.Mol:
    """
    Takes an input peptide and adds double bonds to it based on the CONECT records in the PDB file.
    :param input_pdb:
    :return:
    """
    peptide = Chem.MolFromPDBFile(input_pdb,removeHs=False)
    peptide = Chem.RWMol(peptide)
    with open(input_pdb, 'r') as infile:
        lines = infile.readlines()
        lines = [line.split() for line in lines]
        conect_data = [line for line in lines if "CONECT"  in line]

    bond_counts = Counter()

    for row in conect_data:
        if row[0] == 'CONECT' and len(row) > 2:
            src = row[1]
            for tgt in row[2:]:
                pair = tuple(sorted((src, tgt)))  # ('1', '2') == ('2', '1')
                bond_counts[pair] += 1

    # Extract only duplicates
    duplicate_bonds = {bond: count for bond, count in bond_counts.items() if count > 1}

    # Output
    for atom_ids, num_bonds in duplicate_bonds.items():
        if num_bonds == 4: # meaning CONECT 1 2, CONECT 1 2, CONECT 2 1, CONECT 2 1 == double bond
            src, tgt = atom_ids
            bond = peptide.GetBondBetweenAtoms(int(src) - 1, int(tgt) - 1)
            bond.SetBondType(Chem.BondType.DOUBLE)

    peptide = peptide.GetMol()
    try:
        Chem.SanitizeMol(peptide)
    except ValueError as e:
        logger.error(f"Error sanitizing molecule: {e}")
        raise RuntimeError("Failed to sanitize the molecule after adding double bonds. Please check the input PDB file for correctness.")
    return peptide


def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def boltzmann_weight_distances(og_name,working_dir) -> None:
    """
    Calculate and store the Boltzmann-weighted energies based on molecular conformations.

    This function performs multiple operations including creating molecular
    representations, calculating amide group distances for conformations, computing
    Boltzmann-weighted matrices of distances, and saving results to a CSV file. The
    process involves reading input SMILES files, generating molecular conformations,
    and determining their Boltzmann-weighted contributions.

    :param og_name: the name to call the BWdistance matrix
    :type og_name: str
    :param working_dir: The file path to the working directory where input files
                        are located and output files are saved. Subdirectories for
                        intermediate data may also be created under this path.
    :type working_dir: str
    :return: None
    """

    os.chdir(working_dir)
    names = get_split_files()

    if not os.path.exists(f"{og_name}-BWdistances.csv") or config["rerun"]["distances"]:
        peptides = [add_double_bonds_to_pdb(f"{name}-out-template.pdb") for name in names]
        peptide = peptides[0]
        amideGroups = add_amides(peptides[0])
        temp_working_dir = working_dir + f"/{og_name}_Conformations"
        os.chdir(temp_working_dir) #working in conforamtions folder
        distances = []
        for conformation_xyz in natsorted(os.listdir(temp_working_dir)):
            model_num = conformation_xyz.split("_")[-3]
            if is_int(model_num):
                model_num = int(model_num)
                peptide = peptides[model_num-1]

            peptide = load_xyz_coords(peptide,f"{conformation_xyz}")

            distances.append(get_amide_distances(amideGroups,peptide))
        os.chdir(working_dir)

        boltzmann_matrix = boltzmann_weighted_average(distances, working_dir,og_name)
        print(len(boltzmann_matrix))
        df = pd.DataFrame(boltzmann_matrix)
        df.to_csv(working_dir+f'/{og_name}-BWdistances.csv', index=False, header=False)
        logger.info(f"The BW distance matrix for {og_name} is saved in {working_dir}/{og_name}-BWdistances.csv")
        index = random.randint(0, len(distances) - 1)
        element = distances[index]
        logger.info(f"Test matrix of distances is Conformation number {index}:\n{element}")






def create_template_model(input_pdb, output_pdb) -> None:
    """
    Creates a singular template model to be used for distances/feature calculations based on the input PDB file that originally contains all PDB structures.
    :param input_pdb: Path to the input `.pdb` file containing atomic coordinate data.
    :param output_pdb:
    :return:
    """
    with open(input_pdb, 'r') as infile:
        lines = infile.readlines()

    first_model = []
    conect_lines = []
    inside_model = False

    for line in lines:
        if line.startswith("MODEL"):
            inside_model = True
        if inside_model:
            first_model.append(line)
        if line.startswith("ENDMDL") and inside_model:
            break
    for line in lines:
        if line.startswith("CONECT"):
            conect_lines.append(line)


    # Append CONECT lines to the end of the first model
    first_model += conect_lines

    with open(output_pdb, 'w') as outfile:
        outfile.writelines(first_model)

    logger.info(f"The template model was created and saved in {output_pdb}.")


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
            try:
                conf.SetAtomPosition(i, (x, y, z))
            except ValueError as e:
                logger.error(f"Error setting atom position for line {i+1} in {xyz_path}: {e}")
                raise RuntimeError(f"Failed to set atom position for line {i+1} in {xyz_path}. Please check the file format.")

    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol


def get_dihedral_angle(p1s, p2s, p3s, p4s):
    """Calculate the dihedral angle between four points in 3D space."""
    # Convert to numpy arrays
    p1, p2, p3, p4 = map(np.array, (p1s, p2s, p3s, p4s))

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

def calculate_dihedrals(residue,mol) -> list[float]:
    """
    Calculates dihedral angles for a given residue in a mol object. For the given possible residues, this works.
    :param residue: List of IDs that correspond to the residue in the mol object
    :param mol: Chem.Mol object
    :return:
    """

    if len(residue) == 5: #n-terminus case (normal)  _NCC_
        temp_dihedrals = []
        p1= mol.GetConformer().GetAtomPosition(residue[0])
        p2 = mol.GetConformer().GetAtomPosition(residue[1])
        p3= mol.GetConformer().GetAtomPosition(residue[2])
        p4 = mol.GetConformer().GetAtomPosition(residue[3])
        temp_dihedrals.append(get_dihedral_angle(p1,p2,p3,p4))
        temp_dihedrals.append(5000)
        p1 = mol.GetConformer().GetAtomPosition(residue[1])
        p2 = mol.GetConformer().GetAtomPosition(residue[2])
        p3 = mol.GetConformer().GetAtomPosition(residue[3])
        p4 = mol.GetConformer().GetAtomPosition(residue[4])
        temp_dihedrals.append(get_dihedral_angle(p1, p2, p3, p4))
        return temp_dihedrals

    if len(residue) == 6: #n-terminus case (abnormal)  _NCCC_
        temp_dihedrals = []
        p1= mol.GetConformer().GetAtomPosition(residue[0])
        p2= mol.GetConformer().GetAtomPosition(residue[1])
        p3=mol.GetConformer().GetAtomPosition(residue[2])
        p4 = mol.GetConformer().GetAtomPosition(residue[3])
        temp_dihedrals.append(get_dihedral_angle(p1,p2,p3,p4))
        p1= mol.GetConformer().GetAtomPosition(residue[1])
        p2= mol.GetConformer().GetAtomPosition(residue[2])
        p3=mol.GetConformer().GetAtomPosition(residue[3])
        p4 = mol.GetConformer().GetAtomPosition(residue[4])
        temp_dihedrals.append(get_dihedral_angle(p1,p2,p3,p4))
        p1 = mol.GetConformer().GetAtomPosition(residue[2])
        p2 = mol.GetConformer().GetAtomPosition(residue[3])
        p3 = mol.GetConformer().GetAtomPosition(residue[4])
        p4 = mol.GetConformer().GetAtomPosition(residue[5])
        temp_dihedrals.append(get_dihedral_angle(p1, p2, p3, p4))
        return temp_dihedrals

    return None


def boltzmann_weight_dihedrals(name,working_dir) -> None:
    os.chdir(working_dir)
    ##if not os.path.exists(f"{name}-BWDihedralNormalized.csv") or config["rerun"]["dihedrals"] :
    if not os.path.exists(f"{name}-BWDihedralNormalized.csv") or config["rerun"]["dihedrals"]:
        peptide_normalized_dihedrals = []
        names = get_split_files()
        peptides = [add_double_bonds_to_pdb(f"{name}-out-template.pdb") for name in names]
        mol = peptides[0]
        #add embeding?
        n_terminus = find_n_terminus(mol)

        nitrogen_order = bfs_traversal(mol, n_terminus)
        possible_amides1 = mol.GetSubstructMatches(Chem.MolFromSmarts('[N]C[C](=O)'))
        possible_amides2 = mol.GetSubstructMatches(Chem.MolFromSmarts('[N]CC[C](=O)'))
        possible_amides= possible_amides1+possible_amides2
        possible_nitrogens = [amide[0] for amide in possible_amides]
        ordered_residues = [
            amide
            for n in nitrogen_order
            for amide in possible_amides
            if amide[0] == n and n in possible_nitrogens
        ]
        # this gets the amides and, the 4 atoms we found. just need to find out which atoms we take hte dihedral of again


            ##now to check cases. if any neighbor is a hydrogen, use the other neighbor





        """n_terminus_residue_normal = mol.GetSubstructMatches(Chem.MolFromSmarts('[NH2]C[C](=O)[N]'))
        n_terminus_residue_abnormal = mol.GetSubstructMatches(Chem.MolFromSmarts('[NH2]CC[C](=O)[N]'))
        normal_residues = mol.GetSubstructMatches(Chem.MolFromSmiles('C(=O)NCC(=O)N'))
        abnormal_residues = mol.GetSubstructMatches(Chem.MolFromSmiles('C(=O)NCCC(=O)N'))
        all_residues = normal_residues + abnormal_residues + n_terminus_residue_normal + n_terminus_residue_abnormal
        all_residues = list(all_residues)
        print(len(all_residues))
        ##########################
        # get rid of asparagine
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
                    all_residues.remove(residue)"""

        for conformation_xyz in natsorted(os.listdir(f"{name}_Conformations")):
            model_num = conformation_xyz.split("_")[-3]

            if is_int(model_num):
                model_num = int(model_num)
                mol = peptides[model_num - 1]
            mol.RemoveAllConformers()
            mol = load_xyz_coords(mol, f"{working_dir}/{name}_Conformations/{conformation_xyz}")
            conformation_dihedrals = [] #contains (phi,theta,psi) 6 times, 1 for each residue
            for residue in ordered_residues:
                c_residue = [a for a in residue]
                n_id = c_residue[0]
                c_id = c_residue[-2]
                new_atom_to_add_end = None
                ## we want to get the last carbons neighbor
                for neighbor in mol.GetAtomWithIdx(c_id).GetNeighbors():
                    if neighbor.GetIdx() not in c_residue:
                        new_atom_to_add_end = neighbor.GetIdx()
                        break
                # we dont use oxygen anyway, to just make it the carbon neighbor
                c_residue[-1] = new_atom_to_add_end
                new_atom_to_add_start = None
                n_neighbors = []
                for neighbor in mol.GetAtomWithIdx(n_id).GetNeighbors():
                    if neighbor.GetIdx() not in c_residue:
                        n_neighbors.append(mol.GetAtomWithIdx(neighbor.GetIdx()))

                # this will get the atom id that is not a hydrogen, however if there 2 hydrogens, we will use the one with the highest dihedral.
                if "H" in [atom.GetSymbol() for atom in n_neighbors]:
                    if n_neighbors[0].GetSymbol() == "H" and n_neighbors[1].GetSymbol() == "H":
                        hydrogen_dihedral_1 = calculate_dihedrals([n_neighbors[0].GetIdx()] + c_residue, mol)
                        hydrogen_dihedral_2 = calculate_dihedrals([n_neighbors[1].GetIdx()] + c_residue, mol)
                        if hydrogen_dihedral_1[0] < hydrogen_dihedral_2[0]:
                            new_atom_to_add_start = n_neighbors[0].GetIdx()
                        else:
                            new_atom_to_add_start = n_neighbors[1].GetIdx()
                    else:
                        new_atom_to_add_start = [atom.GetIdx() for atom in n_neighbors if atom.GetSymbol() != "H"][0]
                else: #if there is no hydrogen, pick non methyl carbon
                    for nitrogen_neighbor in n_neighbors:
                        for carbon_neighbor in nitrogen_neighbor.GetNeighbors():
                            if mol.GetAtomWithIdx(carbon_neighbor.GetIdx()).GetSymbol() != "H":
                                new_atom_to_add_start = nitrogen_neighbor.GetIdx()
                                break
                c_residue = [new_atom_to_add_start] + c_residue
                conformation_dihedrals.append(calculate_dihedrals(c_residue,mol))
            #convert each angle to sin/cos components, and add flag = 1 if row contains 5000
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
        df = pd.DataFrame(boltzmann_matrix)
        df.to_csv(working_dir+f'/{name}-BWDihedralNormalized.csv', index=False, header=False)
        logger.info(f"The BW dihedral matrix for {name} is saved in {working_dir}/{name}-BWDihedralNormalized.csv")
        index = random.randint(0, len(peptide_normalized_dihedrals) - 1)
        element = peptide_normalized_dihedrals[index]
        logger.info(f"Test matrix of dihedrals is Conformation number {index}:\n{element}")

def boltzmann_weighted_average(values, working_dir, name) -> np.ndarray:
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


def smart_parse_token(token: str) -> float | str:
    try:
        return float(token)
    except ValueError:
        return token


def remove_duplicates(smiles_lines_all,names_lines_all,percents_all) -> tuple[list,list,list]:
    smile_to_indices = defaultdict(list)
    for i, smile in enumerate(smiles_lines_all):
        smile_to_indices[smile].append(i)
    indices_to_keep = set()
    for indices in smile_to_indices.values():
        if len(indices) == 1:
            indices_to_keep.add(indices[0])
        else:
            best_index = max(indices, key=lambda x: float(percents_all[x])) #first index of the [6,12,18]
            indices_to_keep.add(best_index)

    indices_to_remove = [i for i in range(len(smiles_lines_all)) if i not in indices_to_keep]

    smiles_lines_all = [j for i,j in enumerate(smiles_lines_all) if i not in indices_to_remove]
    names_lines_all = [j for i,j in enumerate(names_lines_all) if i not in indices_to_remove]
    percents_all = [j for i,j in enumerate(percents_all) if i not in indices_to_remove]
    return smiles_lines_all, names_lines_all, percents_all



def create_new_descriptor(descriptor_name,og_name,working_dir) -> None:
    os.chdir(working_dir)
    #not if want to make new one
    if not os.path.exists(f"{og_name}_{descriptor_name}.csv") or config["rerun"]["new_descriptors"]: 

        working_dir = os.getcwd()
        peptide_descriptors = []


        peptides = [add_double_bonds_to_pdb(f"{name}-out-template.pdb") for name in get_split_files()]
        mol = peptides[0]
        amide_groups = add_amides(mol)

        for conformation_xyz in natsorted(os.listdir(f"{og_name}_Conformations")):
            model_num = conformation_xyz.split("_")[-3]
            if is_int(model_num):
                model_num = int(model_num)
                mol = peptides[model_num - 1]

            peptide = load_xyz_coords(mol, f"{working_dir}/{og_name}_Conformations/{conformation_xyz}")

            # here, put the function of what you want to calculate for each conformation
            # function of "peptide"
            peptide_descriptors.append(side_chain_descriptors(amide_groups,peptide))





        logger.info(f"The {descriptor_name} descriptors for {og_name} are saved in {working_dir}/{og_name}_{descriptor_name}.csv")
        peptide_boltzmann = boltzmann_weighted_average(peptide_descriptors,working_dir,og_name)
        peptide_boltzmann = peptide_boltzmann.reshape(len(peptide_boltzmann),-1)
        print(len(peptide_boltzmann[0]))
        df = pd.DataFrame(peptide_boltzmann)
        df.to_csv(f'{og_name}_{descriptor_name}.csv', index=False, header=False)


def read_file(file_name) -> list[str]:
    with open(file_name, 'r') as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines]
        if '\t' in lines[0]: #if its gonna be a list of lists
            lines = [[smart_parse_token(x) for x in line.strip().split('\t')] for line in lines]
        return lines


def make_submol_from_atom_ids(mol, atom_ids) -> Chem.Mol:
    """
    Create a submolecule from a list of atom indices, preserving coordinates and properties.

    Args:
        mol (Chem.Mol): RDKit molecule with 3D coordinates.
        atom_ids (list of int): Atom indices to include in the submolecule.

    Returns:
        Chem.Mol: Submolecule with copied coordinates and properties.
    """
    atom_ids_set = set(atom_ids)

    # Map old atom idx → new atom idx
    old_to_new = {}
    new_mol = Chem.RWMol()

    # Copy atoms
    for i, aid in enumerate(atom_ids):
        old_atom = mol.GetAtomWithIdx(aid)
        new_atom = Chem.Atom(old_atom.GetAtomicNum())
        new_atom.SetFormalCharge(old_atom.GetFormalCharge())
        new_atom.SetChiralTag(old_atom.GetChiralTag())
        new_atom.SetHybridization(old_atom.GetHybridization())
        new_atom.SetNumExplicitHs(old_atom.GetNumExplicitHs())
        new_idx = new_mol.AddAtom(new_atom)
        old_to_new[aid] = new_idx

    # Copy bonds only between included atoms
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if a1 in atom_ids_set and a2 in atom_ids_set:
            new_mol.AddBond(old_to_new[a1], old_to_new[a2], bond.GetBondType())

    # Add conformer with same coordinates
    conf = mol.GetConformer()
    new_conf = Chem.Conformer(len(atom_ids))
    for i, aid in enumerate(atom_ids):
        pos = conf.GetAtomPosition(aid)
        new_conf.SetAtomPosition(i, pos)

    new_mol = new_mol.GetMol()
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(new_conf, assignId=True)
    for atom in new_mol.GetAtoms():
        if atom.GetIsAromatic() and not atom.IsInRing():
            atom.SetIsAromatic(False)  # "Unmark" invalid aromatic atoms

    try:
        Chem.SanitizeMol(new_mol)
    except ValueError as e:
        logger.error(f"Error sanitizing submolecule: {e}")
        raise RuntimeError(f"Failed to sanitize submolecule in make_submol_from_atom_ids: {atom_ids}. Please check the input molecule for correctness.")
    return new_mol


def molecular_descriptors(peptide) -> list[list[float]]:
    descriptors_for_peptide = []
    descriptors_to_calculate = {

        "RadiusOfGyration": Descriptors3D.RadiusOfGyration,
        "Asphericity": Descriptors3D.Asphericity,
        "InertialShapeFactor": Descriptors3D.InertialShapeFactor,
        "Eccentricity": Descriptors3D.Eccentricity,
        "SpherocityIndex": Descriptors3D.SpherocityIndex,
    }
    results = {name: func(peptide) for name, func in descriptors_to_calculate.items()}
    descriptors_for_peptide.append(list(results.values()))
    return descriptors_for_peptide

def side_chain_descriptors(amidegroups,peptide):
    descriptors_for_peptide = []
    descriptors_to_calculate = {
        # Your existing descriptors
        "Radius": Descriptors3D.RadiusOfGyration, #1
        "Asphericity": Descriptors3D.Asphericity, #2
        "InertialShapeFactor": Descriptors3D.InertialShapeFactor, #3
        "Eccentricity": Descriptors3D.Eccentricity, #4
        "SpherocityIndex": Descriptors3D.SpherocityIndex, #5

        # Molecular properties
        "MolLogP": Descriptors.MolLogP,  # Partition coefficient #6
        "MolMR": Descriptors.MolMR,  # Molar refractivity #7
        "HeavyAtomCount": Descriptors.HeavyAtomCount, #8
        "NumHAcceptors": Descriptors.NumHAcceptors,  # H-bond acceptors #9
        "NumHDonors": Descriptors.NumHDonors,  # H-bond donors #10
        "NumRotatableBonds": Descriptors.NumRotatableBonds, #11
        "TPSA": Descriptors.TPSA, #13

        # Electronic descriptors
        "MaxEStateIndex": Descriptors.MaxEStateIndex, #14
        "MinEStateIndex": Descriptors.MinEStateIndex, #15
        "MaxAbsEStateIndex": Descriptors.MaxAbsEStateIndex, # 15
        "MinAbsEStateIndex": Descriptors.MinAbsEStateIndex, #16
        #17*6 descriptors

    }
    for amide_group in amidegroups:
        residue_ids = amide_group.getResidue1()
        residue = make_submol_from_atom_ids(peptide, residue_ids)
        results = {name: func(residue) for name, func in descriptors_to_calculate.items()}


        descriptors_for_peptide.append(list(results.values()))
        if amide_group.getResidue2() is not None:## CHECKING THE LAST AMIDE GROUP, IT IS THE ONLY ONE WITH 2 RESIDUES
            residue_ids = amide_group.getResidue2()
            residue = make_submol_from_atom_ids(peptide, residue_ids)
            for atom in residue.GetAtoms():
                if atom.GetIsAromatic() and not atom.IsInRing():
                    atom.SetIsAromatic(False)
            Chem.SanitizeMol(residue)


            results = {name: func(peptide) for name, func in descriptors_to_calculate.items()} ### I FOUND THE ISSUE HERE!!!!

            descriptors_for_peptide.append(list(results.values()))
    return descriptors_for_peptide







