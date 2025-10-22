import json
import re

def is_categorical(arr, unique_threshold = 20):
    """
    Determine if a list of values behaves like a categorical variable.

    Parameters:
    arr (list): Input list.
    unique_threshold (int or float): Max number or ratio of unique values to consider the list categorical.

    Returns:
    bool: True if the list is considered categorical, False otherwise.
    """

    if len(arr) == 0:
        return True

    # Unique values
    unique_values = len(set(arr))
    if unique_threshold < 1:
        is_few_uniques = unique_values / len(arr) <= unique_threshold
    else:
        is_few_uniques = unique_values <= unique_threshold

    return is_few_uniques


def dumps_compact_lists(obj, indent=4):
    pretty = json.dumps(obj, indent=indent)
    
    # Match any list that spans multiple lines
    def compact_list(match):
        # Extract the list content and remove newlines/extra spaces
        items = match.group(1).strip().splitlines()
        compacted_items = [item.strip().rstrip(',') for item in items if item.strip()]
        return '[' + ','.join(compacted_items) + ']'

    # Replace multi-line lists with compact single-line ones
    compacted = re.sub(r'\[\s*((?:.|\n)*?)\s*\]', compact_list, pretty)
    
    return compacted


