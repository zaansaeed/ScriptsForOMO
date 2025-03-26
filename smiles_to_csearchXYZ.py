import csv
import os
import subprocess
import time

from rdkit.Chem import AllChem
from rdkit import Chem
import pandas as pd
from collections import deque
import math



def split_xyz(input_file):
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


schrodinger_path ="/opt/schrodinger/suites2024-3/"
# Define the working directory (where results will be stored)
main_dir = os.path.abspath("/Users/zaansaeed/Peptides")
os.chdir(main_dir)
smiles_input_file = "all_peptides.smi"
names_input_file = "all_names.txt"


with open(smiles_input_file, "r") as f:
    smiles_lines = f.readlines()
with open(names_input_file, "r") as f:
    names_lines = f.readlines()

# Process each SMILES string individually in Smiles_input_file and also add name
for i, line in enumerate(smiles_lines):
    parts = line.strip().split()
    smile_string = parts[0]
    name = names_lines[i].strip()
    folder_name = "Peptide_"+name

    if not os.path.exists(folder_name):     #if the folder already exists, skip - otherwise continue
        os.mkdir(folder_name)
    working_dir = main_dir+'/'+folder_name

    os.chdir(working_dir)

    if not os.path.exists(f"{name}.smi"):
        temp_smi = os.path.join(working_dir, f"{name}.smi")

        with open(temp_smi, "w") as temp_file:
            temp_file.write(f"{smile_string}\n")

    if not os.path.exists(name+".mae"):
        #os.system("/opt/schrodinger/suites2024-3/ligprep -ismi " + name +".smi -omae " + name + ".mae")#converts smile to .mae
        subprocess.run(["/opt/schrodinger/suites2024-3/ligprep", "-ismi", f"{name}.smi", "-omae", f"{name}.mae"])
        while not os.path.exists(f"{name}.mae"):
            print("waiting for smile to .mae")
            time.sleep(1)








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
                    print("Conf search completed")
                    break
                prev_size = new_size
            time.sleep(60)
            print(f"Waiting for conf search for {name}...")




    #convert to pdb
    if not os.path.exists(f"{name}-out.pdb"):
        #os.system(f"/opt/schrodinger/suites2024-3/utilities/structconvert {name}-out.mae {name}-out.pdb")
        subprocess.run(["/opt/schrodinger/suites2024-3/utilities/structconvert", f"{name}-out.mae",f"{name}-out.pdb"])

        while not os.path.exists(f"{name}-out.pdb") or os.path.getsize(f"{name}-out.pdb") != prev_size:
            if os.path.exists(f"{name}-out.pdb"):
                new_size = os.path.getsize(f"{name}-out.pdb")
                if new_size == prev_size:
                    print("converted to pdb")
                    break
                prev_size = new_size
            time.sleep(30)
            print(f"waiting to convert output to pdb {name}...")



    #convert to xyz
    if not os.path.exists(f"{name}-out.xyz"):
        os.system(f"obabel -ipdb {name}-out.pdb -O {name}-out.xyz")

        while not os.path.exists(f"{name}-out.xyz") or os.path.getsize(f"{name}-out.xyz") != prev_size:
            if os.path.exists(f"{name}-out.xyz"):
                new_size = os.path.getsize(f"{name}-out.xyz")
                if new_size == prev_size:
                    print("converted to xyz")
                    break
                prev_size = new_size
            time.sleep(10)
            print(f"waiting to convert pdb to xyz {name}...")




    if not os.path.exists(f"{name}_Conformations"):
        os.mkdir(f"{name}_Conformations")
        temp_working_dir = working_dir + f"/{name}_Conformations"
        os.chdir(temp_working_dir)
        split_xyz(working_dir+f"/{name}-out.xyz")

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
        self.atom_IDs = atom_IDs
        self.C = self.atom_IDs[0]
        self.O = self.atom_IDs[1]
        self.N = self.atom_IDs[2]

        nitrogen = Peptide.GetAtomWithIdx(self.N)
        neighbors = nitrogen.GetNeighbors()
        hydrogen_id = None

        for neighbor in neighbors:
            if neighbor.GetSymbol() == 'H':
                hydrogen_id = neighbor.GetIdx()
                break
        self.H = hydrogen_id

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

def addAmides(input_peptide):
    amideGroups = []
    n_terminus = Chem.MolFromSmarts('[N;H2]')
    matches = input_peptide.GetSubstructMatches(Chem.MolFromSmarts('[N][C;X4][C;X3](=[O])-[N]'))
    matches = [tpl[2:] for tpl in matches]
    n_terminus_match = input_peptide.GetSubstructMatch(n_terminus)
    n_terminus_id = n_terminus_match[0]
    bfs_order = bfs_traversal(input_peptide, n_terminus_id)
    i = 1
    for nID in bfs_order:
        for match in matches:
            if nID in match:
                amideGroups.append(AmideGroup(match, i, input_peptide))
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


def boltzmann(values, properties_of_each_conformer):
    numerator = 0
    denominator = 0
    boltzmann_results = []
    new_array = []
    for amide_array in range(len(values)):
        new_array = []
        for amide_row in range(len(values[0])):
            new_row = []
            for amide_col in range(len(values[0][0])):
                e_term = 0
                denominator = 0
                numerator = 0
                answer = 0
                for k in range(len(properties_of_each_conformer)):
                    e_term = math.exp(
                        -(properties_of_each_conformer[k]['Energies']) / (298 * 8.314 * 10 ** -3))
                    denominator += e_term

                    numerator += e_term * values[k][amide_row][amide_col]
                    answer = numerator / denominator
                new_row.append(answer)
            new_array.append(new_row)
    boltzmann_results.append(new_array)

    return boltzmann_results

for item in os.listdir(main_dir):
    if item.startswith('Peptide_'):

        name = item[8:]
        working_dir = main_dir+'/'+item
        os.chdir(working_dir)

        if not os.path.exists(f"{name}-BWdistances.csv"):
            print(name)
            with open(f"{name}.smi", "r") as f:
                smiles_string = f.readlines()[0]

            peptide = Chem.AddHs(Chem.MolFromSmiles(smiles_string))
            AllChem.EmbedMolecule(peptide)
            amideGroups = addAmides(peptide)
            temp_working_dir = working_dir + f"/{name}_Conformations"
            os.chdir(temp_working_dir)
            distances = []
            for conformation_xyz in os.listdir(temp_working_dir):
                if conformation_xyz.endswith('.xyz'):
                    atom_coordinates = xyz_to_array(f"{temp_working_dir}/{conformation_xyz}")
                    distances.append(getAmideDistances(amideGroups,atom_coordinates))
            os.chdir(working_dir)
            boltzmann_matrix = boltzmann(distances, pd.read_csv(working_dir+f'/{name}-energies.csv').to_dict(orient="records"))
            with open(f'{name}-BWdistances.csv', 'w',newline='') as file:
                writer = csv.writer(file)
                writer.writerows(boltzmann_matrix)


    os.chdir(main_dir)

















