import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
from itertools import combinations
import anndata as ad

def deg_within_clusters(adata, cluster_col='cluster', group_col='group'):
    """Calculate DEGs for all pairwise group comparisons within each cluster"""
    results = {}
    
    clusters = adata.obs[cluster_col].unique()
    
    for cluster in clusters:
        print(f"Processing cluster {cluster}...")
        
        # Subset to current cluster
        cluster_mask = adata.obs[cluster_col] == cluster
        adata_cluster = adata[cluster_mask].copy()
        
        groups = adata_cluster.obs[group_col].unique()
        
        if len(groups) < 2:
            continue
            
        results[cluster] = {}
        
        # Get all pairwise combinations of groups
        group_pairs = list(combinations(groups, 2))
        
        for group1, group2 in group_pairs:
            print(f"  Comparing {group1} vs {group2}")
            
            # Subset to these two groups
            group_mask = adata_cluster.obs[group_col].isin([group1, group2])
            adata_pair = adata_cluster[group_mask].copy()
            
            # Calculate DEGs
            sc.tl.rank_genes_groups(adata_pair, group_col, 
                                  groups=[group1], reference=group2,
                                  method='wilcoxon')
            
            # Extract results
            deg_df = sc.get.rank_genes_groups_df(adata_pair, group=group1)
            
            results[cluster][f"{group1}_vs_{group2}"] = {
                'deg_df': deg_df,
                'adata': adata_pair
            }
    
    return results

def create_pseudobulk(adata, cluster_col, sample_col, condition_col=None, min_cells=10):
    """
    Create pseudo-bulk profiles by aggregating counts per sample within each cluster
    """
    # Get unique clusters and samples
    clusters = adata.obs[cluster_col].unique()
    samples = adata.obs[sample_col].unique()
    
    if condition_col is not None:
        sample_to_condition = adata.obs.set_index(sample_col)[condition_col].to_dict()
    
    pseudobulk_data = []
    pseudobulk_obs = []
    pseudobulk_vars = adata.var.copy()
    
    for cluster in clusters:
        cluster_mask = adata.obs[cluster_col] == cluster
        
        for sample in samples:
            sample_mask = adata.obs[sample_col] == sample
            cell_mask = cluster_mask & sample_mask
            
            n_cells = np.sum(cell_mask)
            
            if n_cells >= min_cells:
                # Aggregate counts for this sample-cluster combination
                if sparse.issparse(adata.X):
                    sample_counts = adata[cell_mask].X.sum(axis=0).A1
                else:
                    sample_counts = adata[cell_mask].X.sum(axis=0)
                
                pseudobulk_data.append(sample_counts)
                
                # Create observation metadata
                obs_row = {
                    'sample': sample,
                    'cluster': cluster,
                    'n_cells': n_cells,
                    'condition': sample_to_condition[sample] if condition_col is not None else "NA"
                }
                
                pseudobulk_obs.append(obs_row)
    
    # Create pseudo-bulk AnnData object
    pseudobulk_data = np.array(pseudobulk_data)
    pseudobulk_obs = pd.DataFrame(pseudobulk_obs)
    ## combine sample, cluster, condition to form unique index
    pseudobulk_obs.index = pseudobulk_obs.apply(lambda row: f"{row['sample']}_{row['cluster']}_{row['condition']}", axis=1)
    
    pseudobulk_adata = sc.AnnData(
        X=pseudobulk_data,
        obs=pseudobulk_obs,
        var=pseudobulk_vars
    )
    
    # Add sample and cluster as categorical variables
    pseudobulk_adata.obs['sample'] = pseudobulk_adata.obs['sample'].astype('category')
    pseudobulk_adata.obs['cluster'] = pseudobulk_adata.obs['cluster'].astype('category')
    if condition_col is not None:
        pseudobulk_adata.obs['condition'] = pseudobulk_adata.obs['condition'].astype('category')
    
    return pseudobulk_adata

def pseudobulk_deg_pairwise(adata, cluster_col, sample_col, condition_col, 
                           min_cells=10, layer=None):
    """
    Perform pseudo-bulk DEG for all pairwise condition comparisons within each cluster
    
    # Run pairwise analysis
    pairwise_results = pseudobulk_deg_pairwise(
        adata,
        cluster_col='cluster',
        sample_col='sample', 
        condition_col='condition'
    )
    
    # Access results for a specific cluster
    cluster_results = results['c1']
    pseudobulk_data = cluster_results['pseudobulk_adata']

    # Get DEG dataframe
    deg_df = sc.get.rank_genes_groups_df(pseudobulk_data, group=None)

    # Filter significant genes (adjust threshold as needed)
    significant_genes = deg_df[deg_df['pvals_adj'] < 0.05]

    # Get top upregulated genes
    top_upregulated = deg_df[deg_df['logfoldchanges'] > 0].head(10)

    # Get top downregulated genes  
    top_downregulated = deg_df[deg_df['logfoldchanges'] < 0].head(10)

    print(f"Significant DEGs in cluster c1: {len(significant_genes)}")
    print("Top upregulated:", top_upregulated['names'].tolist())
    print("Top downregulated:", top_downregulated['names'].tolist())
    """
    results = {}
    
    clusters = adata.obs[cluster_col].unique()
    
    for cluster in clusters:
        print(f"\n=== Processing Cluster {cluster} ===")
        
        # Subset to current cluster
        cluster_mask = adata.obs[cluster_col] == cluster
        adata_cluster = adata[cluster_mask].copy()
        
        # Create pseudo-bulk
        try:
            pseudobulk = sc.get.aggregate(
                adata_cluster,
                by=sample_col,
                func='sum',
                layer=layer,
                min_count=min_cells
            )
            
            # Add metadata
            sample_to_condition = adata_cluster.obs.groupby(sample_col)[condition_col].first()
            pseudobulk.obs[condition_col] = pseudobulk.obs_names.map(sample_to_condition)
            pseudobulk.obs[condition_col] = pseudobulk.obs[condition_col].astype('category')
            pseudobulk.obs['cluster'] = cluster
            
            # Get conditions and their sample counts
            condition_counts = pseudobulk.obs[condition_col].value_counts()
            print(f"  Condition counts: {dict(condition_counts)}")
            
            # Skip if not enough samples per condition
            if condition_counts.min() < 2:
                print(f"  Not enough replicates for DEG in cluster {cluster}")
                continue
            
            # Normalize
            sc.pp.normalize_total(pseudobulk, target_sum=1e4)
            sc.pp.log1p(pseudobulk)
            
            # Perform DEG for all pairwise comparisons
            conditions = pseudobulk.obs[condition_col].cat.categories
            cluster_results = {}
            
            for i, cond1 in enumerate(conditions):
                for cond2 in conditions[i+1:]:
                    comparison_name = f"{cond1}_vs_{cond2}"
                    print(f"  Comparing {comparison_name}")
                    
                    # Subset to these two conditions
                    condition_mask = pseudobulk.obs[condition_col].isin([cond1, cond2])
                    pseudobulk_subset = pseudobulk[condition_mask].copy()
                    
                    # Ensure we still have both conditions after subsetting
                    if len(pseudobulk_subset.obs[condition_col].unique()) < 2:
                        print(f"    Warning: Missing conditions after subsetting")
                        continue
                    
                    # Perform DEG
                    sc.tl.rank_genes_groups(
                        pseudobulk_subset, 
                        condition_col, 
                        groups=[cond1], 
                        reference=cond2,
                        method='t-test'
                    )
                    
                    # Extract results
                    deg_df = sc.get.rank_genes_groups_df(pseudobulk_subset, group=cond1)
                    
                    cluster_results[comparison_name] = {
                        'deg_df': deg_df,
                        'pseudobulk_adata': pseudobulk_subset,
                        'condition1': cond1,
                        'condition2': cond2
                    }
            
            results[cluster] = cluster_results
            
        except Exception as e:
            print(f"  Error in cluster {cluster}: {e}")
            continue
    
    return results
