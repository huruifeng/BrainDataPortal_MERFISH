# Load required libraries
library(Seurat)
library(SeuratDisk)

# Set input and output file paths
h5ad_file <- "data/MERFISH_Demo_s73.h5ad"      # Replace with your .h5ad file path
seurat_file <- "data/MERFISH_Demo_s73.rds"    # Replace with desired Seurat .rds output path

# Convert h5ad to h5seurat (intermediate format)
Convert(h5ad_file, dest = "h5seurat", overwrite = TRUE, )

# Load h5seurat file as Seurat object
seurat_obj <- LoadH5Seurat(sub("\\.h5ad$", ".h5seurat", h5ad_file))

# Save Seurat object as .rds
saveRDS(seurat_obj, seurat_file)