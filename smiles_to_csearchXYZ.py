import os
import subprocess
import time


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
        os.system("/opt/schrodinger/suites2024-3/ligprep -ismi " + name +".smi -omae " + name + ".mae")#converts smile to .mae
        time.sleep(30)

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

    if not os.path.exists(f"{name}-out.mae"): #if there is no input, then run confsearch
        os.system(f"/opt/schrodinger/suites2024-3/macromodel {name}")
        time.sleep(60)

    #convert to pdb
    if not os.path.exists(f"{name}-out.pdb"):
        os.system(f"/opt/schrodinger/suites2024-3/utilities/structconvert {name}-out.mae {name}-out.pdb")
    #convert to xyz
    if not os.path.exists(f"{name}-out.xyz"):
        os.system(f"obabel -ipdb {name}-out.pdb -O {name}-out.xyz")

    def split_xyz(input_file):
        with open(input_file, "r") as f:
            lines = f.readlines()
        with open(input_file, "r") as f:
            num_coords = int(lines[0].split()[0])
            num_cords = num_coords+2 #num of cords + 2 becuase of the length and name of peptide
            conformations = []
            while True:
                lines = [f.readline() for i in range(num_cords)]
                if not any(lines):
                    break
                conformations.append(lines)
        for i, conf in enumerate(conformations):
            output_file = os.path.join(temp_working_dir, f"{name}_Conformation_{i+1}.xyz")
            with open(output_file, "w") as f:
                f.writelines(conf)



    if not os.path.exists(f"{name}_Conformations"):
        os.mkdir(f"{name}_Conformations")
        temp_working_dir = working_dir + f"/{name}_Conformations"
        os.chdir(temp_working_dir)
        split_xyz(working_dir+f"/{name}-out.xyz")
    os.chdir(main_dir)






