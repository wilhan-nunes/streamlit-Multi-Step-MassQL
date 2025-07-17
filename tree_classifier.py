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
    Find the single most specific path that matches the given classifications.
    Returns the most complete match possible, including the bile acid category name.

    Args:
        matches (list): List of classification matches for a compound
        classification_tree (dict): Hierarchical classification structure

    Returns:
        dict: Results showing the most specific satisfied path with bile acid category
    """

    # Extract all paths from the classification tree
    all_paths = extract_all_paths(classification_tree)

    matches_set = set(matches)  # Convert to set once for efficiency

    best_match = None
    all_matches = []
    best_match_length = 0
    best_match_bile_acid_category = None

    for path in all_paths:
        # Keep the first element (bile acid category name) for the result
        bile_acid_category = path[0]
        # Resets best_match to look for other possible classifications
        if bile_acid_category != best_match_bile_acid_category:
            if best_match:
                all_matches.append([best_match_bile_acid_category] + best_match)
            best_match = None
            best_match_length = 0
        # Remove the first element for matching logic
        path_without_category = path[1:]

        # Find the most complete match by checking from full path backwards
        for i in range(1, len(path_without_category)):
            current_subset = path_without_category[:i]

            # Check if this subset of the path has all elements in matches
            if set(current_subset).issubset(matches_set) and best_match != current_subset:
                # If this is a better (longer) match, update our best match
                if len(current_subset) > best_match_length:
                    best_match = current_subset
                    best_match_length = len(current_subset)
                    best_match_bile_acid_category = bile_acid_category
                # break

    # Append the final best_match if it exists
    if best_match:
        all_matches.append([best_match_bile_acid_category] + best_match)

    # Prepare results
    if all_matches:
        if len(all_matches) >= 2 and all([lst[-1].endswith("stage2") for lst in all_matches]):
            pass
        elif len(all_matches) >= 2:
            all_matches = [max(all_matches, key=len)]
        results = {
            'satisfied_paths': all_matches,
            'most_specific_path': max([lst for lst in all_matches], key=len),
            'bile_acid_category': best_match_bile_acid_category,
            'path_length': best_match_length,
        }
    else:
        results = {
            'satisfied_paths': [],
            'most_specific_path': None,
            'bile_acid_category': None,
            'path_length': 0,
        }

    return results
if __name__ == '__main__':
    from utils import bile_acid_tree
    # classification = ['Monohydroxy', "", 'Monohydroxy_stage1', 'Monohydroxy_stage2', 'Mono-7b-OH']
    # classification = ['Dihydroxy_stage1', 'Di-3,6-OH', 'Mono-3a-OH', 'Di-3,7-OH']
    all_paths = extract_all_paths(bile_acid_tree)
    # query_validation = """Di-3,12a-OH|Di-7,12a-OH;Tri-3,6b,7a-OH|Tri-3,6a,7b-OH;Di-3,6-OH;Dihydroxy_stage1"""
    query_validation = """Di-3,12a-OH;Dihydroxy_stage2;Trihydroxy_stage2;Di-3,6-OH;Dihydroxy_stage1;Trihydroxy_stage1"""
    # query_validation = """Dihydroxy_stage2;Di-3,12a-OH|Di-7,12a-OH;Tri-3,6b,7a-OH|Tri-3,6a,7b-OH;Di-3,6-OH;Dihydroxy_stage1"""
    # query_validation = """Monohydroxy_stage1;Dihydroxy_stage1;Dihydroxy_stage2;Monohydroxy_stage2"""


    results = check_classification_paths(query_validation.split(";"), bile_acid_tree)
    print(results)
