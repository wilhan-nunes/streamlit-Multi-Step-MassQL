def extract_all_paths(tree, current_path=[]):
    """
    Recursively extract all possible paths from root to leaves in the classification tree.

    Args:
        tree (dict): The classification tree structure
        current_path (list): Current path being built

    Returns:
        list: All possible paths from root to leaves
    """
    paths = []

    for key, value in tree.items():
        new_path = current_path + [key]

        if isinstance(value, dict) and value:
            # Continue down the tree
            paths.extend(extract_all_paths(value, new_path))
        else:
            # Leaf node - add the complete path
            paths.append(new_path)

    return paths


def check_classification_paths(matches, classification_tree):
    """
    Check if all combinations of different levels in the classification tree
    are subsets of the matches for a specific compound.

    Args:
        matches (list): List of classification matches for a compound
        classification_tree (dict): Hierarchical classification structure

    Returns:
        dict: Results showing which paths are satisfied and which are not
    """

    # Extract all paths from the classification tree
    all_paths = extract_all_paths(classification_tree)

    # Check each path
    results = {
        'satisfied_paths': [],
        'unsatisfied_paths': [],
        'summary': {}
    }

    for path in all_paths:
        if set(path[1:]).issubset(set(matches)):
            results['satisfied_paths'].append(path)
        else:
            results['unsatisfied_paths'].append(path)

    # Generate summary
    total_paths = len(all_paths)
    satisfied_count = len(results['satisfied_paths'])

    results['summary'] = {
        'total_paths': total_paths,
        'satisfied_paths': satisfied_count,
        'unsatisfied_paths': total_paths - satisfied_count,
        'satisfaction_rate': satisfied_count / total_paths if total_paths > 0 else 0
    }

    return results
