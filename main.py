from rdkit import Chem
from ML_functions import *
from natsort import natsorted
from collections import defaultdict
from rdkit.Chem import Descriptors3D, Descriptors
from functions import add_amides
from visualization import visualize
from functions import boltzmann_weighted_average
from visualization import visualize_molecule_with_highlights
from rdkit.Chem import AllChem



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
            x, y, z = map(float, tokens[1:4])
            conf.SetAtomPosition(i, (x, y, z))

    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)
    return mol


def create_new_descriptor(descriptor_name,directory_of_peptides,descriptor_funcs):
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
                peptide = Chem.MolFromSmiles(smiles_string)
                peptide = Chem.AddHs(peptide)
                AllChem.EmbedMolecule(peptide)
                Chem.MolToPDBFile(peptide, f"/Users/zaansaeed/Desktop/{name}.pdb")
                for conformation_xyz in natsorted(os.listdir(f"{name}_Conformations")):
                    if conformation_xyz.endswith('.xyz'):  # working within 1 conformer
                        mol = load_xyz_coords(peptide, f"{working_dir}/{name}_Conformations/{conformation_xyz}")
                        #here, put the function of what you want to calcualte for each conformation
                        amide_groups = add_amides(mol)
                        peptide_descriptors.append(side_chain_descriptors(amide_groups,descriptor_funcs,mol))


                peptide_boltzmann = boltzmann_weighted_average(peptide_descriptors,working_dir,name)
                print(peptide_boltzmann.shape)
                peptide_boltzmann = peptide_boltzmann.reshape(len(peptide_boltzmann),-1)
                print(peptide_boltzmann)
                df = pd.DataFrame(peptide_boltzmann)
                df.to_csv(f'{peptide_name}_{descriptor_name}.csv', index=False, header=False)
            os.chdir(main_dir)
def make_submol_from_atom_ids(mol, atom_ids):
    """
    Create a submolecule from a list of atom indices, preserving coordinates and properties.

    Args:
        mol (Chem.Mol): RDKit molecule with 3D coordinates.
        atom_ids (list of int): Atom indices to include in the submolecule.

    Returns:
        Chem.Mol: Submolecule with copied coordinates and properties.
    """
    atom_ids_set = set(atom_ids)

    # Map old atom idx → new atom idx
    old_to_new = {}
    new_mol = Chem.RWMol()

    # Copy atoms
    for i, aid in enumerate(atom_ids):
        old_atom = mol.GetAtomWithIdx(aid)
        new_atom = Chem.Atom(old_atom.GetAtomicNum())
        new_atom.SetFormalCharge(old_atom.GetFormalCharge())
        new_atom.SetChiralTag(old_atom.GetChiralTag())
        new_atom.SetHybridization(old_atom.GetHybridization())
        new_atom.SetNumExplicitHs(old_atom.GetNumExplicitHs())
        new_idx = new_mol.AddAtom(new_atom)
        old_to_new[aid] = new_idx

    # Copy bonds only between included atoms
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtomIdx()
        a2 = bond.GetEndAtomIdx()
        if a1 in atom_ids_set and a2 in atom_ids_set:
            new_mol.AddBond(old_to_new[a1], old_to_new[a2], bond.GetBondType())

    # Add conformer with same coordinates
    conf = mol.GetConformer()
    new_conf = Chem.Conformer(len(atom_ids))
    for i, aid in enumerate(atom_ids):
        pos = conf.GetAtomPosition(aid)
        new_conf.SetAtomPosition(i, pos)

    new_mol = new_mol.GetMol()
    new_mol.RemoveAllConformers()
    new_mol.AddConformer(new_conf, assignId=True)
    for atom in new_mol.GetAtoms():
        if atom.GetIsAromatic() and not atom.IsInRing():
            atom.SetIsAromatic(False)  # "Unmark" invalid aromatic atoms

    Chem.SanitizeMol(new_mol)
    return new_mol

def side_chain_descriptors(amidegroups,descriptors_to_calculate,mol):
    descriptors_for_peptide = []

    for amide_group in amidegroups:
        residue_ids = amide_group.getResidue1()
        residue = make_submol_from_atom_ids(mol, residue_ids)

        results = {name: func(residue) for name, func in descriptors_to_calculate.items()}


        descriptors_for_peptide.append(list(results.values()))
        if amide_group.getResidue2() is not None:## CHECKING THE LAST AMIDE GROUP, IT IS THE ONLY ONE WITH 2 RESIDUES
            residue_ids = amide_group.getResidue2()
            residue = make_submol_from_atom_ids(mol, residue_ids)
            for atom in residue.GetAtoms():
                if atom.GetIsAromatic() and not atom.IsInRing():
                    atom.SetIsAromatic(False)
            Chem.SanitizeMol(residue)


            results = {name: func(mol) for name, func in descriptors_to_calculate.items()}

            descriptors_for_peptide.append(list(results.values()))
    print(len(descriptors_for_peptide))
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
    print(smile_to_indices)
    indices_to_keep = set()
    for indices in smile_to_indices.values():
        if len(indices) == 1:
            indices_to_keep.add(indices[0])
        else:
            best_index = max(indices, key=lambda x: float(percents_all[x])) #first index of the [6,12,18]
            indices_to_keep.add(best_index)

    indices_to_remove = [i for i in range(len(smiles_lines_all)) if i not in indices_to_keep]

    smiles_lines_all = [j for i,j in enumerate(smiles_lines_all) if i not in indices_to_remove]
    names_lines_all = [j for i,j in enumerate(names_lines_all) if i not in indices_to_remove]
    percents_all = [j for i,j in enumerate(percents_all) if i not in indices_to_remove]
    return smiles_lines_all, names_lines_all, percents_all

#computes your specified descriptor for each conformation
def compute_global_descriptors(mol):
    descriptor_funcs = {
        'Eccentricity': Descriptors3D.Eccentricity(mol)
    }
    return list(descriptor_funcs.values())

def main():
    #
    main_dir = "/Users/zaansaeed/Desktop/NewPeptides"
    os.chdir(main_dir)
    descriptor_funcs = {
        # Your original descriptors
        "Radius": Descriptors3D.RadiusOfGyration,
        "BertzCT": Descriptors.BertzCT,
        "Kappa1": Descriptors.Kappa1,
        "MolLogP": Descriptors.MolLogP,
        "TPSA": Descriptors.TPSA,
        "Eccentricity": Descriptors3D.Eccentricity,
        "Asphericity": Descriptors3D.Asphericity,
        "SpherocityIndex": Descriptors3D.SpherocityIndex,
        "MolWt": Descriptors.MolWt,
        "NumRotatableBonds": Descriptors.NumRotatableBonds,
    }


    create_new_descriptor('side_chain_descriptors', main_dir,descriptor_funcs)



    smiles_lines_all = read_file("all_peptides.smi") #list of smiles
    names_lines_all = read_file("all_names.txt") # list of names
    percents_all = read_file("percents.txt")# list of lists of percents
    smiles_lines_all = sort_by_names_alphabetically(names_lines_all,smiles_lines_all) #sorts by names
    percents_all = sort_by_names_alphabetically(names_lines_all,percents_all) #sorts by names
    names_lines_all = sort_by_names_alphabetically(names_lines_all,names_lines_all) #changes order of names_lines, must be last

    smiles_final, names_final, percents_final = remove_duplicates(smiles_lines_all,names_lines_all,percents_all) #sorted, with no duplicates
    features = ["BWDihedralNormalized","BWdistances"]

    X = create_X(main_dir,names_final,features)
    #BWdistances , BWDihedralNormalized, side_chain_descriptors
    os.chdir(main_dir)
    Y = create_Y(percents_final)
    print(X.shape,Y.shape)
    os.chdir("/Users/zaansaeed/PycharmProjects/pythonProject/ScriptsForOMO")
    plot_Y_distribution(Y)

    #run_elasticnet(X,Y,5,0.2)
    #run_RFR(X,Y,5,0.20)
    #run_GBR(X,Y,0.2,5)

    #run_SVR(X,Y,5,0.2)
    #run_NN(X,Y,0.2,5)

    #visualize("elasticnet_model.joblib","X.csv","y.csv",descriptor_funcs)

if __name__ == "__main__":
    main()








