from matplotlib import pyplot as plt
from sklearn.tree import plot_tree
import pandas as pd
from natsort import natsorted
from rdkit.Chem import AllChem, FragmentOnBonds
from rdkit import Chem
from functions import add_amides
import py3Dmol
import tempfile
import webbrowser
import os
import re
from joblib import load
import csv

def visualize_model(pipeline, X, Y):
    model = pipeline.named_steps['model']

    # Check if the pipeline has polynomial features
    has_poly = 'poly' in pipeline.named_steps

    if has_poly:
        # Handle polynomial case
        poly = pipeline.named_steps['poly']
        original_feature_names = [f"Feature {i + 1}" for i in range(X.shape[1])]

        # Get polynomial feature names
        if hasattr(poly, 'get_feature_names_out'):
            feature_names = poly.get_feature_names_out(original_feature_names)
        else:
            # Fallback for older scikit-learn versions
            X_transformed = poly.fit_transform(X)
            feature_names = [f"Poly_Feature_{i + 1}" for i in range(X_transformed.shape[1])]
    else:
        # Handle non-polynomial case
        feature_names = [f"Feature {i + 1}" for i in range(X.shape[1])]

    if hasattr(model, 'feature_importances_'):  # For tree-based models
        importances = model.feature_importances_
    elif hasattr(model, 'coef_'):  # For linear models
        importances = abs(model.coef_)

    else:
        raise ValueError("Model doesn't support feature importance visualization")

    # Create DataFrame with features and their importances
    feat_importances = pd.DataFrame({
        'feature': feature_names,
        'importance': importances
    })

    # Sort by importance in descending order and get top 20 features
    feat_importances = feat_importances.sort_values('importance', ascending=False).head(20)

    plt.figure(figsize=(15, 8))
    plt.bar(range(len(feat_importances)), feat_importances['importance'])

    if hasattr(model, 'feature_importances_'):
        plt.title('Top 20 Feature Importances')
    else:
        plt.title('Top 20 Feature Coefficients (Absolute Values)')

    plt.xlabel('Features')
    plt.ylabel('Importance')

    # Add feature names as x-axis labels
    plt.xticks(range(len(feat_importances)),
               [f"{feature }" for feature in feat_importances['feature']],
               rotation=90,
               ha='right')

    # Adjust layout to prevent label cutoff
    plt.tight_layout()

    # Show the plot
    plt.show()

    # Print the detailed mapping of feature numbers to actual features
    print("\nFeature Number Mapping:")
    print("=" * 80)
    for i, (feature, importance) in enumerate(zip(feat_importances['feature'], feat_importances['importance'])):
        print(f"F{i + 1}: {feature:<60} {importance:.4f}")


def analyze_feature_ranges(model, X):
    from sklearn.inspection import partial_dependence, permutation_importance
    import numpy as np

    if hasattr(model, 'named_steps'):
        # For pipeline (e.g., ElasticNet with transformers)
        X_transformed = X
        for _, transformer in model.named_steps.items():
            if transformer != model.named_steps['model']:
                X_transformed = transformer.transform(X_transformed)
        model = model.named_steps['model']

    feature_names = [f"Feature {i + 1}" for i in range(X.shape[1])]

    # Get feature importance using permutation importance
    if hasattr(model, 'feature_importances_'):
        importances = model.feature_importances_
    elif hasattr(model, 'coef_'):
        importances = abs(model.coef_)
    else:
        raise ValueError("Model doesn't support feature importance visualization")

    # Create sorted indices based on importance
    sorted_idx = np.argsort(importances)[::-1]  # Reverse to get descending order
    sorted_features = [feature_names[i] for i in sorted_idx]

    # Get top 10 features
    top_features_idx = sorted_idx[:10]  # Take first 10 from sorted list
    feature_ranges = {}


    for idx in range(len(feature_names)):
        plt.figure(figsize=(10, 6))

        # Calculate partial dependence
        feature_values = np.linspace(X[:, idx].min(), X[:, idx].max(), 50)
        # Reshape feature values for partial_dependence input
        X_temp = X.copy()
        pdp_values = []

        # Manual calculation of partial dependence
        for val in feature_values:
            X_temp[:, idx] = val
            if hasattr(model, 'predict'):
                predictions = model.predict(X_temp)
            else:
                predictions = model.predict_proba(X_temp)[:, 1]
            pdp_values.append(np.mean(predictions))

        pdp_values = np.array(pdp_values)

        # Plot partial dependence
        plt.plot(feature_values, pdp_values, 'b-', label='Partial Dependence')

        # Find the optimal range (where target value is in top 25%)
        threshold = np.percentile(pdp_values, 75)
        optimal_mask = pdp_values >= threshold
        optimal_ranges = feature_values[optimal_mask]

        if len(optimal_ranges) > 0:
            plt.axvspan(optimal_ranges.min(), optimal_ranges.max(),
                        alpha=0.2, color='green', label='Optimal Range')

        #bottom 25%
        threshold = np.percentile(pdp_values, 25)
        least_optimal_mask = pdp_values <= threshold
        least_optimal_ranges = feature_values[least_optimal_mask]
        if len(least_optimal_ranges) > 0:
            plt.axvspan(least_optimal_ranges.min(), least_optimal_ranges.max(),
                        alpha=0.2, color='red', label='Least Optimal Range')
        ##plot
        plt.xlabel(feature_names[idx])
        plt.ylabel('Partial dependence')
        plt.title(f'{feature_names[idx]} (Importance: {importances[idx]:.4f})')
        plt.grid(True)
        plt.legend()
        plt.show()

        if idx in top_features_idx:
            print(f"\n{feature_names[idx]}:")
            print(f"Optimal range: [{optimal_ranges.min():.3f}, {optimal_ranges.max():.3f}]")
            print(f"Least Optimal Range: [{least_optimal_ranges.min():.3f}, {least_optimal_ranges.max():.3f}]")
            print(f"Current data range: [{X[:, idx].min():.3f}, {X[:, idx].max():.3f}]")
            print(f"Mean target value in optimal range: {pdp_values[optimal_mask].mean():.3f}")
            print(f"Mean target value in least optimal range: {pdp_values[least_optimal_mask].mean():.3f}")

        feature_ranges[feature_names[idx]] = {
            "optimal_range": [float(optimal_ranges.min()), float(optimal_ranges.max())],
            "least_optimal_range": [float(least_optimal_ranges.min()), float(least_optimal_ranges.max())],
            "current_range": [float(X[:, idx].min()), float(X[:, idx].max())],
            "mean_target_in_optimal_range": float(pdp_values[optimal_mask].mean()),
            "mean_target_in_least_optimal_range": float(pdp_values[least_optimal_mask].mean()),
            "importance": float(importances[idx])

        }

    feature_ranges = dict((k, feature_ranges[k]) for k in natsorted(feature_ranges.keys()))
    return importances, sorted_features, feature_ranges







def generate_feature_map(atom1,atom2,feature_blocks,feature_ranges):
    """
    feature_blocks: list of tuples
        Each tuple = (feature_type: str, shape: int or (int, int))

    Returns:
        feature_index_map: dict with keys = feature names, values = column indices
    """
    feature_index_map = {}
    idx = 0

    for feature_type, shape in feature_blocks:
        if isinstance(shape, int):  # 1D array
            for i in range(shape):
                key = f"{feature_type}_{i+1}"
                feature_index_map[f"Feature_{idx+1}_{key}"] = feature_ranges[f"Feature {idx+1}"]
                idx += 1
        elif isinstance(shape, tuple) and len(shape) == 2:  # 2D matrix
            for i in range(shape[0]):
                for j in range(shape[1]):
                    key = f"{feature_type}_{atom1}_{i+1}_to_{atom2}_{j+1}"
                    feature_index_map[f"Feature_{idx + 1}_{key}"] = feature_ranges[f"Feature {idx + 1}"]
                    idx += 1
        else:
            raise ValueError(f"Unsupported shape: {shape}")

    return feature_index_map


def visualize_molecule_with_highlights(mol, highlight_atoms):
    """
    Visualizes a 3D molecule and highlights specified atoms by coloring spheres around them.

    Args:
        mol: RDKit Mol object (with or without 3D conformer).
        highlight_atoms: List of atom indices (int) to highlight.
    """
    mol = Chem.AddHs(mol)

    # Generate 3D coords if needed
    if mol.GetNumConformers() == 0:
        if AllChem.EmbedMolecule(mol) != 0:
            raise RuntimeError("3D embedding failed.")
        AllChem.UFFOptimizeMolecule(mol)

    mol_block = Chem.MolToMolBlock(mol)
    viewer = py3Dmol.view(width=2000, height=1000)
    viewer.addModel(mol_block, 'mol')

    # Style the whole molecule as sticks
    viewer.setStyle({'stick': {}})

    # Highlight atoms with spheres (red color and some radius)
    for atom_id in highlight_atoms:
        viewer.addStyle({'serial': atom_id},  # atom serial numbers start at 1 in py3Dmol
                        {'sphere': {'radius': 0.7, 'color': 'red', 'opacity': 0.7}})

    if len(highlight_atoms) == 2:
        atom1 = highlight_atoms[0]
        atom2 = highlight_atoms[1]
        # Add bonds between atoms
        conf = mol.GetConformer()
        p1 = conf.GetAtomPosition(atom1)
        p2 = conf.GetAtomPosition(atom2)
        viewer.addLine({
            'start': {'x': p1.x, 'y': p1.y, 'z': p1.z},
            'end': {'x': p2.x, 'y': p2.y, 'z': p2.z},
            'dashed': True,
            'color': 'black',
            'linewidth': 50
        })

    viewer.zoomTo()


    viewer.setBackgroundColor('white')

    # Save to temp file and open in browser
    html = viewer._make_html()
    with tempfile.NamedTemporaryFile(delete=False, suffix=".html", mode='w', encoding='utf-8') as f:
        f.write(html)
        temp_path = f.name

    print(f"Opening viewer in browser: {temp_path}")
    webbrowser.open('file://' + os.path.realpath(temp_path))


def visualize_peptide_and_save_features(feature_map,arbitrary_peptide_smiles,feature_to_examine):
    peptide = Chem.MolFromSmiles(arbitrary_peptide_smiles)
    peptide = Chem.AddHs(peptide)
    amide_groups = add_amides(peptide)
    pb = Chem.MolToMolBlock(peptide)
    feature = [key for key in feature_map if feature_to_examine in key]
    feature_in_map = feature[0]
    if "side_chain" in feature_in_map:
        side_chain_num = int(feature_in_map.split("_")[4])
        if side_chain_num == 6:
            group = amide_groups[-1]
            residue_ids = group.getResidue2()[1]
            visualize_molecule_with_highlights(peptide,residue_ids)
        else:
            group = amide_groups[side_chain_num-1]
            residue_ids = group.getResidue1()[1]
            visualize_molecule_with_highlights(peptide,residue_ids)
    if "distance" in feature_in_map:
        match = re.search(r'distance_hydrogen_(\d+)_to_oxygen_(\d+)', feature_in_map)
        if match:
            h_group_idx = int(match.group(1))
            o_group_idx = int(match.group(2))
            hydrogen_id = amide_groups[h_group_idx-1].getH()
            oxygen_id = amide_groups[o_group_idx-1].getO()
            visualize_molecule_with_highlights(peptide,[hydrogen_id,oxygen_id])

    all_keys = set(k for sub in feature_map.values() for k in sub)
    fieldnames = ['Feature'] + sorted(all_keys)

    sorted_items = sorted(feature_map.items(), key=lambda x: x[1].get('importance', 0), reverse=True)

    with open('features_sorted.csv', 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for key, subdict in sorted_items:
            row = {'Feature': key}
            # Convert list values to strings
            for k, v in subdict.items():
                row[k] = str(v) if isinstance(v, list) else v
            writer.writerow(row)

def visualize(path_to_model,path_to_X,path_to_y):
    model = load(path_to_model)
    X = pd.read_csv(path_to_X).values
    y = pd.read_csv(path_to_y).values.ravel()
    visualize_model(model,X,y)
    _, _, feature_ranges = analyze_feature_ranges(model, X)
    feature_blocks = [
        ("side_chain", 6),
        ("distance", (5, 5)),
    ]
    # r1c1
    smiles = "CC(C)[C@H](C(N[C@H](CCC(OCC=C)=O)C(NC(C)(C)C(N(c(cc(cc1)C(NC)=O)c1N1C)C1=O)=O)=O)=O)N(C)C([C@H](Cc1c[nH]c2c1cccc2)NC(CN(C)C(CCN)=O)=O)=O"
    feature_map = generate_feature_map("hydrogen", "oxygen", feature_blocks, feature_ranges)
    print(feature_map)
    visualize_peptide_and_save_features(feature_map, smiles, "distance_hydrogen_2_to_oxygen_3")

if __name__ == "__main__":
    visualize("random_forest.joblib","X.csv","y.csv")