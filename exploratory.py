from ScriptsForOMO.functions import add_double_bonds_to_pdb
from ScriptsForOMO.functions import add_amides
from ScriptsForOMO.functions import boltzmann_weight_dihedrals
main_dir = "/Users/zaansaeed/Desktop/PermeabilityDataset"
import os

os.chdir(main_dir+"/"+"Peptide_2")
working_dir = os.getcwd()
mol = add_double_bonds_to_pdb("2_model-out-template.pdb")
amide_groups = add_amides(mol)

os.chdir(main_dir)
boltzmann_weight_dihedrals("2",working_dir)

