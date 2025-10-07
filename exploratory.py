from ScriptsForOMO.functions import add_double_bonds_to_pdb
from ScriptsForOMO.functions import add_amides
from ScriptsForOMO.functions import boltzmann_weight_dihedrals
from natsort import natsorted

import pandas as pd
main_dir = "/Users/zaansaeed/Desktop/PermeabilityDataset"
import os


for folder in natsorted(os.listdir(main_dir)):
    if folder.startswith("Peptide_"):
        os.chdir(main_dir+"/"+folder)
        name = folder.split("_")[1]
        data = pd.read_csv(f"{name}-BWDihedralNormalized.csv",header=None)
        if not data.shape[0] == 6:
            print(name)