from matplotlib import pyplot as plt
from sklearn.tree import plot_tree
import pandas as pd
from natsort import natsorted

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
        plt.xlabel(feature_names[idx])
        plt.ylabel('Partial dependence')
        plt.title(f'{feature_names[idx]} (Importance: {importances[idx]:.4f})')
        plt.grid(True)
        plt.legend()
        plt.show()

        if idx in top_features_idx:
            print(f"\n{feature_names[idx]}:")
            print(f"Optimal range: [{optimal_ranges.min():.3f}, {optimal_ranges.max():.3f}]")
            print(f"Current data range: [{X[:, idx].min():.3f}, {X[:, idx].max():.3f}]")
            print(f"Mean target value in optimal range: {pdp_values[optimal_mask].mean():.3f}")

        feature_ranges[feature_names[idx]] = {
            "optimal_range": [float(optimal_ranges.min()), float(optimal_ranges.max())],
            "current_range": [float(X[:, idx].min()), float(X[:, idx].max())],
            "mean_target": float(pdp_values[optimal_mask].mean()),
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
                feature_index_map[f"Feature {idx+1} - {key}"] = feature_ranges[f"Feature {idx+1}"]
                idx += 1
        elif isinstance(shape, tuple) and len(shape) == 2:  # 2D matrix
            for i in range(shape[0]):
                for j in range(shape[1]):
                    key = f"{feature_type} of {atom1} on Amide {i+1} to {atom2} on Amide {j+1}"
                    feature_index_map[f"Feature {idx + 1} - {key}"] = feature_ranges[f"Feature {idx + 1}"]
                    idx += 1
        else:
            raise ValueError(f"Unsupported shape: {shape}")

    return feature_index_map

feature_blocks = [
    ("Side_Chain", 6),
    ("Distance", (5, 5)),
]

