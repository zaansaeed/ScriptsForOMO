from ScriptsForOMO.functions import add_double_bonds_to_pdb
from ScriptsForOMO.functions import add_amides
from ScriptsForOMO.functions import boltzmann_weight_dihedrals
from natsort import natsorted

import pandas as pd
main_dir = "/Users/zaansaeed/Desktop/PermeabilityDataset"
import os

os.chdir(main_dir+"/"+"Peptide_7459")
working_dir = os.getcwd()
mol = add_double_bonds_to_pdb("7459_model-out-template.pdb")
amide_groups = add_amides(mol)

os.chdir(main_dir)
boltzmann_weight_dihedrals("7459",working_dir)

"""for folder in natsorted(os.listdir(main_dir)):
    if folder.startswith("Peptide_"):
        os.chdir(main_dir+"/"+folder)
        name = folder.split("_")[1]
        data = pd.read_csv(f"{name}-BWDihedralNormalized.csv",header=None)
        print(data.shape,name)"""