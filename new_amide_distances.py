from array import array
from tkinter.font import names

from rdkit import Chem
from ML_functions import *
from natsort import natsorted
from collections import defaultdict
from rdkit.Chem import Descriptors3D, Descriptors, Crippen, Lipinski
import inspect
from ScriptsForOMO.functions import get_amide_distances
from functions import boltzmann
from functions import add_amides

def distance_between_two_atoms(mol,ID1,ID2):
    conf = mol.GetConformer()
    atom_1_pos = conf.GetAtomPosition(ID1)
    atom_2_pos = conf.GetAtomPosition(ID2)
    return np.linalg.norm(atom_1_pos - atom_2_pos)

def getAmideDistances(amideGroups,mol):
    distance_matrix = []
    for i,amid1 in enumerate(amideGroups):
        distances_for_one_amid = []
        for j,amid2 in enumerate(amideGroups):
            if i == j:
                distances_for_one_amid.append([0,0,0])
            else:
                amid1C=amid1.getC()
                amid2C=amid2.getC()
                amid1N=amid1.getN()
                amid2N=amid2.getN()
                amid1O=amid1.getO()
                amid2O=amid2.getO()
                distances_for_one_amid.append([distance_between_two_atoms(mol,amid1C,amid2C),distance_between_two_atoms(mol,amid1O,amid2O),distance_between_two_atoms(mol,amid1N,amid2N)])
        distance_matrix.append(distances_for_one_amid)
    return distance_matrix

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


def create_new_descriptor(descriptor_name,directory_of_peptides):
    main_dir = os.path.abspath(directory_of_peptides)
    os.chdir(main_dir)
    for folder in natsorted(os.listdir(main_dir)):
        if os.path.isdir(folder):
            os.chdir(folder) #currenlty working in Peptide_{name}
            peptide_name  = folder.split('_')[1]
            if not os.path.exists(f"{peptide_name}_{descriptor_name}.csv"): #change to/from not
                working_dir = os.getcwd()
                name = folder.split("_")[1]
                smiles_string = open(f"{name}.smi").read().strip()
                print(name)
                peptide_descriptors = []
                mol = Chem.MolFromSmiles(smiles_string)
                mol = Chem.AddHs(mol)

                for conformation_xyz in natsorted(os.listdir(f"{name}_Conformations")):
                    if conformation_xyz.endswith('.xyz'):  # working within 1 conformer
                        mol.RemoveAllConformers()
                        mol = load_xyz_coords(mol, f"{working_dir}/{name}_Conformations/{conformation_xyz}")
                        #here, put the function of what you want to calcualte for each conformation
                        amideGroups = add_amides(mol)
                        peptide_descriptors.append(side_chain_descriptors(amideGroups))

                peptide_boltzmann = boltzmann(peptide_descriptors,working_dir,name)
                print(peptide_boltzmann.shape)
                peptide_boltzmann = peptide_boltzmann.reshape(len(peptide_boltzmann),-1)
                print(peptide_boltzmann)
                df = pd.DataFrame(peptide_boltzmann)
                df.to_csv(f'{peptide_name}_{descriptor_name}.csv', index=False, header=False)
            os.chdir(main_dir)

def side_chain_descriptors(amidegroups):
    descriptors_for_peptide = []

    for amide_group in amidegroups:
        mol = amide_group.getResidue()[0] #getresidue is (frag_structure, ids of atoms)
        for atom in mol.GetAtoms():
            if atom.GetIsAromatic() and not atom.IsInRing():
                atom.SetIsAromatic(False)
        Chem.SanitizeMol(mol)
        results = {
            "TPSA": Descriptors.TPSA(mol),  # Topological Polar Surface Area
            "Radius" :Descriptors3D.RadiusOfGyration(mol)

        }


        descriptors_for_peptide.append(list(results.values()))

    return descriptors_for_peptide




def smart_parse_token(token: str):
    try:
        return float(token)
    except ValueError:
        return token

def read_file(file_name) -> list[str]:
    with open(file_name, 'r') as f:
        lines = f.readlines()
        lines = [line.strip() for line in lines]
        if '\t' in lines[0]: #if its gonna be a list of lists
            lines = [[smart_parse_token(x) for x in line.strip().split('\t')] for line in lines]
        return lines



def remove_duplicates(smiles_lines_all,names_lines_all,percents_all) -> tuple[list,list,list]:
    smile_to_indices = defaultdict(list)
    for i, smile in enumerate(smiles_lines_all):
        smile_to_indices[smile].append(i)


    indices_to_keep = set()
    for indices in smile_to_indices.values():
        if len(indices) == 1:
            indices_to_keep.add(indices[0])
        else:
            best_index = max(indices, key=lambda x: percents_all[x][0]/100) #first index of the [6,12,18]
            indices_to_keep.add(best_index)

    indices_to_remove = [i for i in range(len(smiles_lines_all)) if i not in indices_to_keep]

    smiles_lines_all = [j for i,j in enumerate(smiles_lines_all) if i not in indices_to_remove]
    names_lines_all = [j for i,j in enumerate(names_lines_all) if i not in indices_to_remove]
    percents_all = [j for i,j in enumerate(percents_all) if i not in indices_to_remove]
    return smiles_lines_all, names_lines_all, percents_all

#computes your specified descriptor for each conformation
def compute_global_descriptors(mol):
    descriptor_funcs = {
        "MolWt": Descriptors.MolWt(mol),
        "HeavyAtomCount": Descriptors.HeavyAtomCount(mol),
        "BertzCT": Descriptors.BertzCT(mol),
        "MolLogP": Descriptors.MolLogP(mol),
        "FractionCSP3": Descriptors.FractionCSP3(mol),
        "TPSA": Descriptors.TPSA(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol),
        "NOCount": Descriptors.NOCount(mol),
        "NumRotatableBonds": Descriptors.NumRotatableBonds(mol),
        "Kappa1": Descriptors.Kappa1(mol),
        "Chi1v": Descriptors.Chi1v(mol),
        "RingCount": Descriptors.RingCount(mol),
        "NumAliphaticRings": Descriptors.NumAliphaticRings(mol),
        "NumAromaticRings": Descriptors.NumAromaticRings(mol),
        "NumSaturatedRings": Descriptors.NumSaturatedRings(mol),
        "BalabanJ": Descriptors.BalabanJ(mol),
        "HallKierAlpha": Descriptors.HallKierAlpha(mol),
    }
    return list(descriptor_funcs.values())

def main():
    main_dir = "/Users/zaansaeed/Peptides"
    os.chdir(main_dir)

    create_new_descriptor('side_chain_descriptors', main_dir)



    smiles_lines_all = read_file("all_peptides.smi") #list of smiles
    names_lines_all = read_file("all_names.txt") # list of names
    percents_all = read_file("percent6-12-18.txt")# list of lists of percents
    smiles_lines_all = sort_by_names_alphabetically(names_lines_all,smiles_lines_all) #sorts by names
    percents_all = sort_by_names_alphabetically(names_lines_all,percents_all) #sorts by names
    names_lines_all = sort_by_names_alphabetically(names_lines_all,names_lines_all) #changes order of names_lines, must be last


    smiles_final, names_final, percents_final = remove_duplicates(smiles_lines_all,names_lines_all,percents_all) #sorted, with no duplicates

    X = create_X(main_dir,names_final,["BWDihedralNormalized"])


    Y = create_Y(percents_final)

    #print(X.shape)
    ### removal of 1s and 0s

    X = X.tolist()
    Y =Y.tolist()

    for i in reversed(range(len(X))):
        if Y[i] == 0 or Y[i] == 1:
            del X[i]
            del Y[i]
    X = np.array(X)
    Y = np.array(Y)
    plot_Y_distribution(Y)
    run_RFR(X, Y, 0.3, 5)
    dummy_RFR(X,Y,0.25)



if __name__ == "__main__":
    main()









