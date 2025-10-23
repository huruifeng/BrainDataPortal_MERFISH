# Load required libraries
import scanpy as sc
import anndata as ad
import sys
import numpy as np


## Base on the current h5ad file, create a new h5ad file with necessary modifications
## current h5ad file: data/20251020_s80n_s80p_demo.h5ad, it has 2 samples: s80n and s80p
## The new h5ad file will be saved to: data/20251020_s80n_s80p_modified.h5ad
## The modifications include:
## 1. Add a new column 'condition' to adata.obs, with values 'normal' for sample 's80n' and 'pathology' for sample 's80p
## 2. split s80n into 5 samples: s80n_1, s80n_2, ..., s80n_5 randomlyly, and add a new column 'sample_new_id' to adata.obs
## 3. split s80p into 5 samples: s80p_1, s80p_2, ..., s80p_5 randomlyly, and add a new column 'sample_new_id' to adata.obs
## 4. save the new h5ad file to data/20251020_s80n_s80p_modified.h5ad

def make_modified_adata(input_file, output_file):
    print(f"Reading h5ad file from {input_file}...")
    adata = ad.read_h5ad(input_file, backed=None)
    
    print("Modifying adata.obs...")
    np.random.seed(42)  # for reproducibility
    
    # Add 'condition' column
    adata.obs['condition'] = adata.obs['sample'].map({'s80n': 'normal', 's80p': 'pathology'})
    
    # Create new sample IDs
    def assign_new_sample_ids(original_sample, n_splits=5):
        cell_indices = np.where(adata.obs['sample'] == original_sample)[0]
        np.random.shuffle(cell_indices)
        split_size = len(cell_indices) // n_splits
        new_ids = []
        for i in range(n_splits):
            start_idx = i * split_size
            end_idx = (i + 1) * split_size if i < n_splits - 1 else len(cell_indices)
            new_id = f"{original_sample}_{i + 1}"
            new_ids.extend([new_id] * (end_idx - start_idx))
        return new_ids
    
    new_sample_ids = []
    for sample in ['s80n', 's80p']:
        new_sample_ids.extend(assign_new_sample_ids(sample))
    
    adata.obs['sample_new_id'] = new_sample_ids
    
    print(f"Saving modified h5ad file to {output_file}...")
    adata.write_h5ad(output_file)
    print("Done.")

    return adata
# Specify input and output files
input_h5ad_file = "data/20251020_s80n_s80p_demo.h5ad"
output_h5ad_file = "data/20251020_s80n_s80p_modified.h5ad"

# Call the function
make_modified_adata(input_h5ad_file, output_h5ad_file)

