import csv
import os
import subprocess
import time

from rdkit.Chem import AllChem
from rdkit import Chem
import pandas as pd
from collections import deque
import math
import numpy as np
from natsort import natsorted



main_dir = os.path.abspath("/Users/zaansaeed/Peptides")

def split_xyz(temp_working_dir,input_file,name):
    with open(input_file, "r") as f:
        lines = f.readlines()
    with open(input_file, "r") as f:
        num_coords = int(lines[0].split()[0])
        num_cords = num_coords + 2  # num of cords + 2 becuase of the length and name of peptide
        conformations = []
        while True:
            lines = [f.readline() for i in range(num_cords)]
            if not any(lines):
                break
            conformations.append(lines)
    for i, conf in enumerate(conformations):
        output_file = os.path.join(temp_working_dir, f"{name}_Conformation_{i + 1}.xyz")
        with open(output_file, "w") as f:
            f.writelines(conf)



def smile_to_mae(smile_string,name): #file with smiles, file wiht names, path to folder
    os.chdir(main_dir)
    # Process each SMILES string individually in Smiles_input_file and also add name

    folder_name = "Peptide_"+name

    if not os.path.exists(folder_name):     #if the folder already exists, skip - otherwise continue
        os.mkdir(folder_name)
    working_dir = main_dir+'/'+folder_name

    os.chdir(working_dir)
    if not os.path.exists(f"{name}.smi"):
        temp_smi = os.path.join(working_dir, f"{name}.smi")
        with open(temp_smi, "w") as temp_file:
            temp_file.write(f"{smile_string}\n")

    if not os.path.exists(f"{name}.mae"):

        #os.system("/opt/schrodinger/suites2024-3/ligprep -ismi " + name +".smi -omae " + name + ".mae")#converts smile to .mae
        subprocess.run(["/opt/schrodinger/suites2024-3/ligprep", "-ismi", f"{name}.smi", "-omae", f"{name}.mae"])
        while not os.path.exists(f"{name}.mae"):
            time.sleep(1)
        print("DONE CONVERTING TO MAE")

    os.chdir(main_dir)







def run_confsearch(name,working_dir):
    os.chdir(working_dir)
    #we now have a .mae file of the smile string we will run conf search on. we need to make confsearch file
    if not os.path.exists(name):
        conf_search_file = os.path.join(working_dir, f"{name}")
        with open(conf_search_file, "w") as f:
            f.write(f"INPUT_STRUCTURE_FILE {name}.mae\n")
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

    if not os.path.exists(f"{name}-out.mae"):
        #os.system(f"/opt/schrodinger/suites2024-3/macromodel {name}") #if there is no input, then run confsearch
        subprocess.run(["/opt/schrodinger/suites2024-3/macromodel", f"{name}"])
        prev_size = -1

        while not os.path.exists(f"{name}-out.mae") or os.path.getsize(f"{name}-out.mae") != prev_size:
            if os.path.exists(f"{name}-out.mae"):
                new_size = os.path.getsize(f"{name}-out.mae")
                if new_size == prev_size:
                    break
                prev_size = new_size
            time.sleep(60)

        print("DONE RUNNING CONFSEARCH")

    os.chdir(main_dir)



def mae_to_pdb(name,working_dir):
    #convert to pdb
    os.chdir(working_dir)
    if not os.path.exists(f"{name}-out.pdb"):
        #os.system(f"/opt/schrodinger/suites2024-3/utilities/structconvert {name}-out.mae {name}-out.pdb")
        subprocess.run(["/opt/schrodinger/suites2024-3/utilities/structconvert", f"{name}-out.mae",f"{name}-out.pdb"])
        prev_size = -1
        while not os.path.exists(f"{name}-out.pdb") or os.path.getsize(f"{name}-out.pdb") != prev_size:
            if os.path.exists(f"{name}-out.pdb"):
                new_size = os.path.getsize(f"{name}-out.pdb")
                if new_size == prev_size:
                    break
                prev_size = new_size
            time.sleep(30)

        print("DONE CONVERTING MAE TO PDB")

    os.chdir(main_dir)



def pdb_to_xyz(name,working_dir):
    #convert to xyz
    os.chdir(working_dir)
    if not os.path.exists(f"{name}-out.xyz"):
        os.system(f"obabel -ipdb {name}-out.pdb -O {name}-out.xyz")
        prev_size = -1
        while not os.path.exists(f"{name}-out.xyz") or os.path.getsize(f"{name}-out.xyz") != prev_size:
            if os.path.exists(f"{name}-out.xyz"):
                new_size = os.path.getsize(f"{name}-out.xyz")
                if new_size == prev_size:
                    break
                prev_size = new_size
            time.sleep(10)

        print("DONE CONVERTING PDB TO XYZ")

    os.chdir(main_dir)


def xyz_to_individual_xyz(name,working_dir):
    os.chdir(working_dir)
    if not os.path.exists(f"{name}_Conformations"):
        os.mkdir(f"{name}_Conformations")
        temp_working_dir = working_dir + f"/{name}_Conformations"
        os.chdir(temp_working_dir)
        split_xyz(temp_working_dir,working_dir+f"/{name}-out.xyz",name)

        print("DONE CONVERTING XYZ TO INDIVIDUAL XYZ")

    os.chdir(main_dir)


def extract_energies_to_csv(name,working_dir):
    os.chdir(working_dir)
    if not os.path.exists(f"{name}-energies.csv"):
        with open(f"{name}.log", "r") as f:
            conformations_list = []
            lines = f.readlines()
            for line in lines:
                if "Conformation " in line:
                    conformations_list.append(line.split()[3])
        df = pd.DataFrame(conformations_list, columns=["Energies"])
        num_energies = len(df)
        row_labels = [f"Conformation{i}" for i in range(1,num_energies+1)]
        df.index =row_labels
        df.to_csv(f"{name}-energies.csv")

        print("DONE EXTRACTING ENERGIES")

    os.chdir(main_dir)




#+====================================================================================

def bfs_traversal(mol, startingID):
    visited = set()
    queue = deque([startingID])
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

    def __init__(self, atom_IDs, group_num, Peptide):
        self.group_num = group_num

        self.N= atom_IDs[0]

        nitrogen = Peptide.GetAtomWithIdx(self.N)
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
        self.H = hydrogen_id
        self.atom_IDs = (self.C,self.O,self.N)


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
    def getGroupNum(self):
        return self.group_num

def find_n_terminus(input_peptide):
    n_terminus_residue_normal = input_peptide.GetSubstructMatches(Chem.MolFromSmarts('[NH2]C[C](=O)'))
    n_terminus_residue_abnormal = input_peptide.GetSubstructMatches(Chem.MolFromSmarts('[NH2]CC[C](=O)'))
    if n_terminus_residue_normal:
        return n_terminus_residue_normal[0][0]
    else:
        return n_terminus_residue_abnormal[0][0]

def addAmides(input_peptide):
    amideGroups = []
    matches1 = input_peptide.GetSubstructMatches(Chem.MolFromSmiles('NCC(=O)N'))# N[C;X4][C;X3](=[O])N
    matches2 = input_peptide.GetSubstructMatches(Chem.MolFromSmiles('NCCC(=O)N'))#N[C;X4][C;X4][C;X3](=[O])N
    matches = matches1 + matches2

    n_terminus = find_n_terminus(input_peptide)

    bfs_order = bfs_traversal(input_peptide, n_terminus)
    i = 1
    used_IDS = []
    for nID in bfs_order:
        for match in matches:
            if nID in match and n_terminus not in match:
                amide = AmideGroup(match, i, input_peptide)
                amide_IDS = amide.getIDs()
                if amide_IDS not in used_IDS:
                    amideGroups.append(AmideGroup(match, i, input_peptide))
                    used_IDS.append(amide_IDS)
                    i += 1

    return amideGroups


def getAmideDistances(amideGroups,atom_coordinates):
    distance_matrix = [[0.0 for _ in range(len(amideGroups))] for _ in range(len(amideGroups))]
    for i,amid1 in enumerate(amideGroups):
        for j,amid2 in enumerate(amideGroups):
            if i == j and amid1.getH() is not None:
                distance_matrix[i][j] = 0.0
            else:
                amid1H = amid1.getH()
                if amid1H is None:
                    distance_matrix[i][j] = -1.0 #prolines have no nitrogen Hydrogen, so the distance is non existant
                else:
                    amid2O = amid2.getO()
                    amid1H_pos = atom_coordinates[amid1H]
                    amid2O_pos = atom_coordinates[amid2O]
                    distance = ((amid1H_pos[0] - amid2O_pos[0])**2 + (amid1H_pos[1] - amid2O_pos[1])**2 + (amid1H_pos[2] - amid2O_pos[2])**2)**0.5
                    distance_matrix[i][j] = distance
    return distance_matrix

def xyz_to_array(xyz_file):
    #oth index is num atoms
    #ith index is ith atom ID
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






def boltzmann_weight_energies(name,working_dir, update_matrices):
    os.chdir(working_dir)
    if not os.path.exists(f"{name}-BWdistances.csv") or update_matrices: #hcange to not
        with open(f"{name}.smi", "r") as f:
            smiles_string = f.readlines()[0]

        peptide = Chem.AddHs(Chem.MolFromSmiles(smiles_string))
        AllChem.EmbedMolecule(peptide)
        amideGroups = addAmides(peptide)

        temp_working_dir = working_dir + f"/{name}_Conformations"
        os.chdir(temp_working_dir)
        distances = []
        for conformation_xyz in natsorted(os.listdir(temp_working_dir)):
            if conformation_xyz.endswith('.xyz'):
                atom_coordinates = xyz_to_array(f"{temp_working_dir}/{conformation_xyz}")
                distances.append(getAmideDistances(amideGroups,atom_coordinates))
        os.chdir(working_dir)
        boltzmann_matrix = boltzmann(distances, working_dir,name)


        df = pd.DataFrame(boltzmann_matrix)
        df.to_csv(working_dir+f'/{name}-BWdistances.csv', index=False, header=False)

    if update_matrices:
        print("UPDATED BOLTZMANN MATRIX")


    os.chdir(main_dir)






################################################
################################################

def load_xyz_coords(mol, xyz_path):
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



def extract_boltzmann_weighted_dihedrals_normalized():
    os.chdir(main_dir)
    for folder in natsorted(os.listdir(main_dir)):
        if os.path.isdir(folder):
            os.chdir(folder)
            working_dir = os.getcwd()
            name = folder.split("_")[1]
            print(name)
            if not os.path.exists(f"{name}-BWdihedrals.csv") or not os.path.exists(f"{name}-BWDihedralNormalized.csv"):
                smiles_string = open(f"{name}.smi").read().strip() #generate the smiles string, currently working in Peptide _XXXX folder
                peptide_normalized_dihedrals = []
                mol = Chem.MolFromSmiles(smiles_string)
                mol = Chem.AddHs(mol)
                mol = load_xyz_coords(mol, f"{working_dir}/{name}_Conformations/{conformation_xyz}")
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
                boltzmann_matrix = boltzmann(peptide_normalized_dihedrals, working_dir,name)


                boltzmann_matrix = boltzmann_matrix.reshape(len(boltzmann_matrix),-1)
                print(boltzmann_matrix)
                df = pd.DataFrame(boltzmann_matrix)
                df.to_csv(working_dir+f'/{name}-BWDihedralNormalized.csv', index=False, header=False)
                print("dihedral calculation done for " +name )
                os.chdir(main_dir)

            os.chdir(main_dir)



def boltzmann(values, working_dir,name):
    energies = pd.read_csv(os.path.join(working_dir, f'{name}-energies.csv'))
    energy_vals = energies['Energies'].values
    print(energy_vals[0])

    R = 8.314e-3  # kJ/mol·K
    T = 298  # Kelvin

    weights = np.exp(-energy_vals / (R * T))
    weights = weights / np.sum(weights)  # Normalize weights

    weighted_sum = np.zeros_like(values[0], dtype=float)
    for i, arr in enumerate(values):
        weighted_sum += weights[i] * np.array(arr)

    return weighted_sum














