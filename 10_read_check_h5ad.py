# Load required libraries
import scanpy as sc
import anndata as ad
import sys

# Function to read and check h5ad file
def read_and_check_h5ad(file_path):
    try:
        print(f"Reading h5ad file from {file_path}...")
        adata = ad.read_h5ad(file_path, backed=None)
        return adata
    except Exception as e:
        print(f"Error reading h5ad file: {e}")
        return None

def save_h5ad_structure(adata, output_file=None):
    """
    Print or save the structure summary of an AnnData (.h5ad) object,
    similar to `str(seurat_obj)` in R.

    Parameters
    ----------
    adata : anndata.AnnData
        The AnnData object to summarize.
    output_file : str or None, optional
        If provided, save output to this file. If None, print to console.
    """
    # Decide output stream
    if output_file is None:
        f = sys.stdout
    else:
        f = open(output_file, "w")

    print("AnnData object summary", file=f)
    print("-" * 60, file=f)
    print(f"Cells (obs): {adata.n_obs}", file=f)
    print(f"Genes (var): {adata.n_vars}", file=f)
    print(f"Main matrix (X): {type(adata.X).__name__} with shape {adata.X.shape}", file=f)
    print("", file=f)

    # Layers
    if adata.layers:
        print(f"Layers ({len(adata.layers)}): {', '.join(adata.layers.keys())}", file=f)
    else:
        print("Layers: None", file=f)

    # obs
    print("", file=f)
    print(f"obs: {len(adata.obs.columns)} columns", file=f)
    for col in adata.obs.columns:
        print(f"  • {col} ({adata.obs[col].dtype})", file=f)

    # var
    print("", file=f)
    print(f"var: {len(adata.var.columns)} columns", file=f)
    for col in adata.var.columns:
        print(f"  • {col} ({adata.var[col].dtype})", file=f)

    # obsm
    print("", file=f)
    if adata.obsm:
        print("obsm:", file=f)
        for k, v in adata.obsm.items():
            print(f"  • {k}: shape {v.shape}", file=f)
    else:
        print("obsm: None", file=f)

    # obsp
    print("", file=f)
    if adata.obsp:
        print("obsp:", file=f)
        for k, v in adata.obsp.items():
            print(f"  • {k}: shape {v.shape}", file=f)
    else:
        print("obsp: None", file=f)

    # uns
    print("", file=f)
    if adata.uns:
        print(f"uns keys ({len(adata.uns.keys())}): {list(adata.uns.keys())}", file=f)
    else:
        print("uns: None", file=f)

    print("-" * 60, file=f)

    # Close file if needed
    if output_file is not None:
        f.close()
        print(f"✅ Saved summary to {output_file}")


# Main function
if __name__ == "__main__":
    h5ad_file = "data/20251020_s80n_s80p_demo.h5ad"  # Replace with your .h5ad file path
    adata = read_and_check_h5ad(h5ad_file)

    if adata is not None:
        print("h5ad file read successfully.")
    else:
        print("Error reading h5ad file.")

    # Save or print structure summary
    save_h5ad_structure(adata, output_file="h5ad_structure.txt")

