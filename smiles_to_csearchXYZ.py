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

    if not os.path.exists(folder_name):     #if the folder already exists, skip - otherwise continue
        os.mkdir(folder_name)
        os.chdir(working_dir+folder_name)

        temp_smi = os.path.join(working_dir, f"{name}.smi")

        with open(temp_smi, "w") as temp_file:
            temp_file.write(f"{smile_string}\n")

        if not os.path.exists(name+".mae"):
            os.system("/opt/schrodinger/suites2024-3/ligprep -ismi " + name +".smi -omae " + name + ".mae")#converts smile to .mae

        #we now have a .mae file of the smile string we will run conf search on. we need to make confsearch file
        if not os.path.exists(name):
            conf_search_file = os.path.join(working_dir, f"{name}")
            with open(conf_search_file, "r") as f:
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







        os.chdir(working_dir)