# %% ==============================
import pandas as pd
import json
from collections import defaultdict
import os
import sys
import numpy as np

# %% ============================================================================
## top 10 marker genes in each cell type
# dataset_folder = sys.argv[1]  # e.g. "datasets/PD5D_MTG_snRNAseq"
# cluster_col = sys.argv[2]  # e.g. "SubCellTypes"
# sex_col = sys.argv[3]  # e.g. "sex"

dataset_folder = "datasets/MERFISH_Demo"
cluster_col = "celltype"
sex_col = None

output_folder = dataset_folder + "/clustermarkers"

marker_gene_file = output_folder + "/cluster_FindAllMarkers.csv"
marker_genes = pd.read_csv(marker_gene_file, index_col=None, header=0)

# Filter for significant genes
filtered_df = marker_genes[marker_genes['p_val_adj'] < 0.05]

# Rank by absolute log2FC and select top 10 per cluster
top_genes = (
    filtered_df
    .assign(abs_log2FC = filtered_df['avg_log2FC'])  ## chenge to abs(filtered_df['avg_log2FC']) will include negative log2FC genes
    .sort_values(['cluster', 'abs_log2FC'], ascending=[True, False])
    .groupby('cluster')
    .head(10)
    .drop(columns='abs_log2FC')  # Optional: remove helper column
)

# Save or inspect the result
top_genes.loc[:,["cluster","gene","avg_log2FC","p_val_adj"]].to_csv(output_folder+'/cluster_markergenes_topN.csv', index=False)

## save marker genes of each cell type into a dictionary
marker_genes_dict = defaultdict(list)
for index, row in top_genes.iterrows():
    marker_genes_dict[row["cluster"]].append([row["gene"], row["avg_log2FC"], row["p_val_adj"]])
# print("marker genes dictionary:")
# print(marker_genes_dict)
# Save the dictionary to a JSON file
with open(output_folder + "/cluster_markergenes_topN.json", "w") as f:
    json.dump(marker_genes_dict, f, indent=4)

# %%============================================================================
pool_genes = top_genes["gene"].unique().tolist()
cluster_list = top_genes["cluster"].unique().tolist()
metadata = pd.read_csv(dataset_folder + "/cellspot_metadata_original.csv", index_col=0, header=0)

#%% ============================================================================
## calculate the percentage of cells where the gene is detected in that cell type
condition_col = "Condition"
pct_detected = {}

marker_genes_df = pd.DataFrame()
for cluster in cluster_list:
    print("==========================")
    print("Processing cluster: ", cluster)
    pct_detected[cluster] = []
    cells_in_cluster = metadata[metadata[cluster_col] == cluster]
    num_cells = len(cells_in_cluster)
    

    conditions = metadata[condition_col].unique().tolist()
    if sex_col is not None and sex_col in metadata.columns: 
        sex = metadata[sex_col].unique().tolist()
        for condition in conditions:
            for s in sex:
                subgroup_counts = {}
                print("Processing condition: ", condition, " and sex: ", s)
                diagnosis_sex_group = cells_in_cluster[(cells_in_cluster[condition_col] == condition) & (cells_in_cluster[sex_col] == s)]
                n_cells = len(diagnosis_sex_group)
                subgroup_counts["condition"] = condition
                subgroup_counts["sex"] = s
                subgroup_counts["count"] = n_cells
                pct_detected[cluster].append(subgroup_counts)  # add the subgroup counts to the list
    else:
        for condition in conditions:
            subgroup_counts = {}
            print("Processing condition: ", condition)
            diagnosis_group = cells_in_cluster[cells_in_cluster[condition_col] == condition]
            n_cells = len(diagnosis_group)
            subgroup_counts["condition"] = condition
            subgroup_counts["count"] = n_cells
            pct_detected[cluster].append(subgroup_counts)  # add the subgroup counts to the list
                
    marker_genes = {}
    for gene in pool_genes:
        marker_genes[gene] = {}
        gene_expr = json.load(open(dataset_folder + "/gene_jsons/"+gene+".json"))
        gene_expr_in_cells =  list(gene_expr.keys())

        cells_with_gene_expr = [cell for cell in gene_expr_in_cells if cell in cells_in_cluster.index]
        num_cells_with_gene_expr = len(cells_with_gene_expr)

        if num_cells_with_gene_expr == 0:
            avg_expr = 0.00
        else:
            # Calculate the average expression value for the cells with gene expression and convert to float
            # avg_expr = np.mean([float(gene_expr[cell]) for cell in cells_with_gene_expr])

            # Use np.nanmean to ignore NaN values in the calculation
            # This will return NaN if all values are NaN, so we need to handle that case
            avg_expr = np.nanmean([float(gene_expr[cell]) for cell in cells_with_gene_expr])
            # If all values are NaN, set avg_expr to 0.00
            if np.isnan(avg_expr):
                avg_expr = 0.00
            else:
                avg_expr = round(avg_expr, 2)
        marker_genes[gene]["avg_expr"] = avg_expr
        marker_genes[gene]["is_marker"] = gene in [g[0] for g in marker_genes_dict[cluster]] 
        marker_genes[gene]["n_expr_cells"] = num_cells_with_gene_expr
    
    marker_df = pd.DataFrame.from_dict(marker_genes, orient='index')
    marker_df[cluster_col] = cluster
    marker_df["cluster_n_cells"] = num_cells

    marker_genes_df = pd.concat([marker_genes_df, marker_df], axis=0)
marker_genes_df["gene"] = marker_genes_df.index
marker_genes_df = marker_genes_df.reset_index(drop=True)
marker_genes_df = marker_genes_df[["gene", cluster_col, "cluster_n_cells"] + [col for col in marker_genes_df.columns if col not in ["gene", cluster_col, "cluster_n_cells"]]]  
marker_genes_df.to_csv(output_folder + "/cluster_markergenes_cellcounts.csv", index=False)

# Save the dictionary to a JSON file
with open(output_folder + "/cluster_cellcounts.json", "w") as f:
    json.dump(pct_detected, f, indent=4)    



# %%
