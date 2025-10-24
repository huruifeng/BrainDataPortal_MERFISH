# %%
## Calculate the  markers in each cluster
import scanpy as sc
import anndata as ad
import os
import pandas as pd

from utils.degs import create_pseudobulk, deg_within_clusters

# %%
## ========== Parameters ==========
adata_file = "data/20251020_s80n_s80p_modified.h5ad"
output_dir = "datasets/MERFISH_Demo/clustermarkers/"
cluster_col = "celltype"
condition_col = "condition"
sample_col = "sample_new_id"

min_cells_per_gene = 10
min_cells_per_group = 5
logfc_threshold = 0.20
adjusted_pval_threshold = 0.05

if not os.path.exists(output_dir):
    os.makedirs(output_dir)

# %%
## ========== Cluster markers (FindAllMarkers) ==========
print("============================================")
print(f"Calculating cluster markers for {adata_file} using column '{cluster_col}'")
adata = ad.read_h5ad(adata_file, backed=None)
sc.tl.rank_genes_groups(adata, groupby=cluster_col, use_raw=True, method='wilcoxon')

markers_df = sc.get.rank_genes_groups_df(adata, group=None)
markers_df.rename(columns={'group': 'cluster', 'names': 'gene', 'logfoldchanges': 'avg_log2FC', 'pvals': 'p_val', 'pvals_adj': 'p_val_adj'}, inplace=True)
markers_df.to_csv(os.path.join(output_dir, "cluster_FindAllMarkers.csv"), index=False)
print(f"✅ Cluster markers (FindAllMarkers) saved to {output_dir}")

# %%
## ========== calculate DEGs within each cell type/cluster ==========
print("============================================")
# calcaulate differential expression within each cell type between the conditions
if condition_col is not None:
    print(f"Calculating DEGs within each cluster using condition column '{condition_col}'")
    deg_results = deg_within_clusters(adata, cluster_col=cluster_col, group_col=condition_col)
    
    ## merge and save results
    print("Merging DEGs from all clusters...")
    all_deg_dfs = []
    topN = 10
    topN_deg_dfs = []
    for cluster, comparisons in deg_results.items():
        for comparison, res in comparisons.items():
            deg_df = res['deg_df']
            deg_df['cluster_DE'] = f"{cluster}.{comparison}"
            deg_df.rename(columns={'names': 'gene', 'logfoldchanges': 'avg_log2FC', 'pvals': 'p_val', 'pvals_adj': 'p_val_adj'}, inplace=True)
            deg_df.drop(columns=['scores'], inplace=True)
            deg_df = deg_df[['cluster_DE', 'gene', 'avg_log2FC', 'p_val', 'p_val_adj']]
            deg_df['avg_log2FC'] = deg_df['avg_log2FC'].round(2)
            deg_df['p_val'] = deg_df['p_val'].apply(lambda x: f"{x:.2e}").astype(float)
            deg_df['p_val_adj'] = deg_df['p_val_adj'].apply(lambda x: f"{x:.2e}").astype(float)
            all_deg_dfs.append(deg_df)
            
            # Get top N upregulated and N downregulated genes by avg_log2FC with p_val_adj < 0.05,
            # topN_deg = deg_df.nlargest(topN, 'avg_log2FC')
            deg_df_filtered = deg_df[deg_df['p_val_adj'] < adjusted_pval_threshold]
            topN_up = deg_df_filtered.nlargest(topN, 'avg_log2FC')
            topN_down = deg_df_filtered.nsmallest(topN, 'avg_log2FC')
            topN_deg = pd.concat([topN_up, topN_down])
            topN_deg_dfs.append(topN_deg)
            
    all_deg_df = pd.concat(all_deg_dfs, ignore_index=True)
    topN_deg_df = pd.concat(topN_deg_dfs, ignore_index=True)
    all_deg_df.to_csv(os.path.join(output_dir, "cluster_DEGs.csv"), index=False)
    topN_deg_df.to_csv(os.path.join(output_dir, "cluster_DEGs_topN.csv"), index=False)
    print(f"✅ Cluster DEGs saved to {output_dir}")


# %%
## ========== pseudo-bulk expression matrix ==========
print("============================================")
# create pseudo-bulk expression matrix
print(f"Creating pseudo-bulk expression matrix from {adata_file} using cluster column '{cluster_col}' and sample column '{sample_col}'")
pb_adata = create_pseudobulk(adata, cluster_col=cluster_col, sample_col=sample_col,condition_col=condition_col, min_cells=1)
pb_expr_file = os.path.join(output_dir, "pb_expr_matrix.csv")
pb_expr_df = pd.DataFrame(pb_adata.X, index=pb_adata.obs_names, columns=pb_adata.var_names)
## round to 2 decimal places
pb_expr_df = pb_expr_df.round(2)
pb_expr_df.to_csv(pb_expr_file)
print(f"✅ Pseudo-bulk expression matrix saved to {pb_expr_file}")

# # %%
# ## validate pseudo-bulk values by comparing with manual aggregation, use gene 'ACADL' as example
# print("Validating pseudo-bulk values by manual aggregation...")
# gene_of_interest = "ACADL"
# gene_df = pb_expr_df[pb_expr_df.index.str.startswith("s80n_5")][gene_of_interest]
# gene_sum = gene_df.sum()
# print(f"Sum of pseudo-bulk expression for gene {gene_of_interest} in sample s80n_5: {gene_sum}")

# %%
## ========== pseudo-bulk DE analysis in each cell type/cluster ==========
print("============================================")
# pseudo-bulk DE analysis in each cell type/cluster
print("Calculating pseudo-bulk differential expression within each cell type...")
pb_deg_results = deg_within_clusters(pb_adata, cluster_col="cluster", group_col=condition_col)

## merge and save results
print("Merging DEGs from all clusters...")
all_pb_deg_dfs = []
topN = 10
topN_pb_deg_dfs = []
for cluster, comparisons in pb_deg_results.items():
    for comparison, res in comparisons.items():
        deg_df = res['deg_df']
        deg_df['cluster_DE'] = f"{cluster}.{comparison}"
        deg_df.rename(columns={'names': 'gene', 'logfoldchanges': 'avg_log2FC', 'pvals': 'p_val', 'pvals_adj': 'p_val_adj'}, inplace=True)
        deg_df.drop(columns=['scores'], inplace=True)
        deg_df = deg_df[['cluster_DE', 'gene', 'avg_log2FC', 'p_val', 'p_val_adj']]
        deg_df['avg_log2FC'] = deg_df['avg_log2FC'].round(2)
        
        ## replace inf values with maximum/minimum values for better visualization
        max_log2fc = deg_df.loc[deg_df['avg_log2FC'] != float('inf'), 'avg_log2FC'].max()
        min_log2fc = deg_df.loc[deg_df['avg_log2FC'] != float('-inf'), 'avg_log2FC'].min()
        deg_df['avg_log2FC'] = deg_df['avg_log2FC'].replace(float('inf'), max_log2fc)
        deg_df['avg_log2FC'] = deg_df['avg_log2FC'].replace(float('-inf'), min_log2fc)

        deg_df['p_val'] = deg_df['p_val'].apply(lambda x: f"{x:.2e}").astype(float)
        deg_df['p_val_adj'] = deg_df['p_val_adj'].apply(lambda x: f"{x:.2e}").astype(float)
        all_pb_deg_dfs.append(deg_df)
        
        # Get top N upregulated and N downregulated genes by avg_log2FC with p_val_adj < 0.05,
        # topN_deg = deg_df.nlargest(topN, 'avg_log2FC')
        deg_df_filtered = deg_df[deg_df['p_val_adj'] < adjusted_pval_threshold]
        topN_up = deg_df_filtered.nlargest(topN, 'avg_log2FC')
        topN_down = deg_df_filtered.nsmallest(topN, 'avg_log2FC')
        topN_pb_deg = pd.concat([topN_up, topN_down])
        topN_pb_deg_dfs.append(topN_pb_deg)

all_pb_deg_df = pd.concat(all_pb_deg_dfs, ignore_index=True)
topN_pb_deg_df = pd.concat(topN_pb_deg_dfs, ignore_index=True)
all_pb_deg_df.to_csv(os.path.join(output_dir, "cluster_pb_DEGs.csv"), index=False)
topN_pb_deg_df.to_csv(os.path.join(output_dir, "cluster_pb_DEGs_topN.csv"), index=False)
print(f"✅ Cluster pseudo-bulk DEGs saved to {output_dir}")


print("All done!")

# %%
