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








