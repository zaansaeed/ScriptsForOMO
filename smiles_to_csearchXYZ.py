import os
import subprocess

from ScriptsForOMO.main import smiles

schrodinger_path ="/opt/schrodinger/suites2024-3/"
# Define the working directory (where results will be stored)
working_dir = os.path.abspath("/Users/zaansaeed/Peptides")

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
    os.mkdir(folder_name)
    os.chdir(working_dir+folder_name)

    temp_smi = os.path.join(working_dir, f"{name}.smi")

    with open(temp_smi, "w") as temp_file:
        temp_file.write(f"{smile_string} {name}\n")
