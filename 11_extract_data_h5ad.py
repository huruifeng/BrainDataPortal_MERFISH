# Load required libraries
import scanpy as sc
import anndata as ad
import os

import pandas as pd

# Function to read and extract data from h5ad file
data_file = "data/20251020_s80n_s80p_demo.h5ad"  # Replace with your .h5ad file path
dataset_name = "MERFISH_Demo"  # Replace with your dataset name


output_folder = "datasets/" + dataset_name + "/"
if not os.path.exists(output_folder):
    os.makedirs(output_folder)

adata = ad.read_h5ad(data_file, backed=None)

## extract expression matrix, save as csv
print("===========================================")
print(f"Extracting data from {data_file} for dataset {dataset_name}")
expr_matrix = adata.X
## Convert expression matrix to DataFrame as long format, with columns: cs_id, Gene, Expression
expr_df = pd.DataFrame.sparse.from_spmatrix(expr_matrix, index=adata.obs_names, columns=adata.var_names)
print(f"Expression check: {expr_df.loc[["s80n_2.0"], ["ACKR2","AQP4","CCN1","COL5A1", "CXCL12"]]}")
expr_df = expr_df.stack().reset_index()
expr_df.columns = ['cs_id', 'Gene', 'Expression']

##  remove zero expression rows
expr_df = expr_df[expr_df['Expression'] > 0]

print(f"Saving expression matrix with shape {expr_df.shape} to {output_folder + 'raw_expr_matrix.csv'}")
expr_df.to_csv(output_folder + "raw_expr_matrix.csv", index=False)
print("Example data:")
print(expr_df.iloc[:5, :5])

## Extract and save obs dataframe
print("===========================================")
print(f"Saving cell metadata with shape {adata.obs.shape} to {output_folder + 'raw_cell_metadata.csv'}")
obs_df = adata.obs
obs_df.to_csv(output_folder + "raw_cell_metadata.csv")

## Extract and save var dataframe
print("===========================================")
print(f"Saving gene metadata with shape {adata.var.shape} to {output_folder + 'raw_genes_metadata.csv'}")
var_df = adata.var
var_df.to_csv(output_folder + "raw_genes_metadata.csv")

## Extract the cell coordinates if available
print("===========================================")
print("Extracting cell coordinates if available...")
if 'spatial' in adata.obsm.keys():
    spatial_matrix = adata.obsm['spatial']
    spatial_df = pd.DataFrame.sparse.from_spmatrix(spatial_matrix, index=adata.obs_names)
    spatial_df.to_csv(output_folder + "raw_cell_coordinates.csv")
    print(f"Saved cell coordinates to {output_folder + 'raw_cell_coordinates.csv'}")
elif 'X_spatial' in adata.obsm.keys():
    spatial_matrix = adata.obsm['X_spatial']
    spatial_df = pd.DataFrame.sparse.from_spmatrix(spatial_matrix, index=adata.obs_names)
    spatial_df.to_csv(output_folder + "raw_cell_coordinates.csv")
    print(f"Saved cell coordinates to {output_folder + 'raw_cell_coordinates.csv'}")
elif 'x' in adata.obs.keys() and 'y' in adata.obs.keys():
    spatial_df = adata.obs[['x', 'y']]
    spatial_df.to_csv(output_folder + "raw_cell_coordinates.csv")
    print(f"Saved cell coordinates to {output_folder + 'raw_cell_coordinates.csv'}")
else:
    print("No cell coordinates found.")

## Extract and save umap coordinates
print("===========================================")
print("Extracting UMAP coordinates if available...")
if 'X_umap' in adata.obsm.keys():
    umap_matrix = adata.obsm['X_umap']
    umap_df = pd.DataFrame(umap_matrix, index=adata.obs_names, columns=['UMAP_1', 'UMAP_2'])
    umap_df.to_csv(output_folder + "raw_umap_embeddings.csv")
    print(f"Saved UMAP coordinates to {output_folder + 'raw_umap_embeddings.csv'}")

## Extract and save obsm data
# print("Extracting obsm data...")
# for key in adata.obsm.keys():
#     obsm_matrix = adata.obsm[key]
#     obsm_df = pd.DataFrame(obsm_matrix, index=adata.obs_names)
#     obsm_df.to_csv(output_folder + f"raw_obsm_{key}.csv")
#     print(f"Saved obsm data for key '{key}' to {output_folder + f'raw_obsm_{key}.csv'}")
