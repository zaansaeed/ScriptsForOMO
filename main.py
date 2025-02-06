from rdkit.Chem import AllChem
from rdkit import Chem
import py3Dmol
import math
import numpy as np
import os
import pandas as pd
import glob
import re


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

def sort_files(s):
    return [int(text) if text.isdigit() else text for text in re.split(r'(\d+)', os.path.basename(s))]



class Peptide:

    peptide_registry ={} # each peptide will store itself in the library

    def __init__(self, smiles_code, PercentCyclization, name, csv_of_properties):
        self.molecule = Chem.AddHs(Chem.MolFromSmiles(smiles_code))
        AllChem.EmbedMolecule(self.molecule)
        self.atom_IDs = [atom.GetIdx() for atom in self.molecule.GetAtoms()]
        self.conformers = []
        self.amideGroups = []
        self.PercentCyclization = PercentCyclization
        self.Name = name
        self.properties = pd.read_csv(csv_of_properties).to_dict(orient="records")
        Peptide.peptide_registry[name] = self #the key is the name, the oobject is the peptide itself

    def getName(self):
        return self.Name

    def getPercentCyclization(self):
        return self.PercentCyclization

    def csearch(self):
        pass # = Conformer(self, xyz)

    def addConformer(self, conformer):
        self.conformers.append(conformer)

    def getConformers(self):
        return self.conformers

    def addAmides(self):
        group = Chem.MolFromSmarts('[C](=[O])-[N]')
        matches = self.molecule.GetSubstructMatches(group)
        if len(matches) == 0:
            return "There are no amide groups in this peptide."
        else:
            for i, IDs in enumerate(matches):
                self.amideGroups.append(AmideGroup(IDs,i,self))
            return "Added " + str(len(matches)) + " Amide groups to GroupIDs."


class Conformer(Peptide):
    def __init__(self, parent_peptide: Peptide, xyz_file):
        self.atom_coordinates = xyz_to_array(xyz_file)
        #self.relative_energy = relative_energy

    def getAmideDistances(self):
        atom_coordinates = xyz_to_array("R1C1_conformation_001.xyz")
        distance_matrix = [[0.0 for _ in range(len(self.amideGroups))] for _ in range(len(self.amideGroups))]
        for i,amid1 in enumerate(self.amideGroups):
            for j,amid2 in enumerate(self.amideGroups):
                if i == j:
                    distance_matrix[i][j] = 0.0
                else:
                    amid1H =amid1.getH()
                    amid2O = amid2.getO()
                    amid1H_pos = atom_coordinates[amid1H]
                    amid2O_pos = atom_coordinates[amid2O]
                    distance = ((amid1H_pos[0] - amid2O_pos[0])**2 + (amid1H_pos[1] - amid2O_pos[1])**2 + (amid1H_pos[2] - amid2O_pos[2])**2)**0.5
                    distance_matrix[i][j] = distance
        return distance_matrix




class AmideGroup:

    def __init__(self, atom_IDs, group_num, Peptide):
        self.group_num = group_num
        self.atom_IDs = atom_IDs
        self.C = self.atom_IDs[0]
        self.O = self.atom_IDs[1]
        self.N = self.atom_IDs[2]

        nitrogen = Peptide.molecule.GetAtomWithIdx(self.N)
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

peptides = []
df = pd.read_csv("PeptideCyclizationSmiles.csv") #readas all peptide data
df = df[['Address','Smiles','Percent cyclization']].dropna() # drops other columns
peptide_data = df.set_index('Address').to_dict(orient="index") #makes data a dictionary accessible by the name of the peptide

for folder_name in os.listdir("/Users/zaansaeed/Google Drive/My Drive/OMO Lab - Peptide Cyclization - Zaan Saeed/Data/Peptide Library"):
    files = sorted(glob.glob(f"{"/Users/zaansaeed/Google Drive/My Drive/OMO Lab - Peptide Cyclization - Zaan Saeed/Data/Peptide Library/" + folder_name}/**/*",recursive=True),key=sort_files)
    peptide_name = folder_name[:4]
    peptide_csv = ""
    for file in files:
        if file.endswith(".csv"):
            peptide_csv = file

    if peptide_name in peptide_data.keys():
        smiles = peptide_data[peptide_name]['Smiles']
        percent = peptide_data[peptide_name]['Percent cyclization']
        peptide = Peptide(smiles, percent, peptide_name,peptide_csv)
        for file in files:
            if file.endswith(".xyz"):
                peptide.addConformer(Conformer(peptide, os.path.abspath(file)))

        peptides.append(peptide)


amide_distances= []
for peptide in peptides:
    peptide.addAmides()
    for conformer in peptide.getConformers():
        conformer.getAmideDistances()



