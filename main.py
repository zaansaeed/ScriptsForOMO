import os
import subprocess
# Set paths
schrodinger_path ="/opt/schrodinger/suites2024-3/"
# Define the working directory (where results will be stored)
working_dir = os.path.abspath("/Users/zaansaeed/maestro_files") #set working directory
os.chdir(working_dir)
# Input SMILES file
print(working_dir)
input_file = "input.smi"

# Read the input SMILES file
with open(input_file, "r") as f:
    lines = f.readlines()

# Process each SMILES string individually
for i, line in enumerate(lines):
    parts = line.strip().split()


    smiles = parts[0]  # First part is the SMILES string
    name = parts[1] if len(parts) > 1 else f"molecule_{i+1}"  # Use name if available, else auto-name

    temp_smi = os.path.join(working_dir, f"{name}.smi")

    # Write SMILES to a temporary .smi file
    with open(temp_smi, "w") as temp_file:
        temp_file.write(f"{smiles} {name}\n")

    # Output .mae file
    output_mae = os.path.join(working_dir, f"{name}.mae")

    # Run LigPrep for the individual molecule
    ligprep_command = ["/opt/schrodinger/suites2024-3/ligprep", "-ismi", name+".smi", "-omae", name+".mae"]

    try:
        print(f"Processing {name}...")
        subprocess.run(ligprep_command,check=True,shell=True)
        print(f"Generated: {output_mae}")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {name}: {e}")

print(f"All .mae files are saved in {working_dir}/")