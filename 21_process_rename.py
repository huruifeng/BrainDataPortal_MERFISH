# %% ==============================
import pandas as pd
import json
import os
import sys

from utils.funcs import is_categorical, dumps_compact_lists

import functools
print = functools.partial(print, flush=True)

print("============================================")
# %% ==============================
# Get the arguments
# dataset_path = sys.argv[1]
# kept_features = sys.argv[2].split(",")
# sample_col = sys.argv[3]
# cluster_col = sys.argv[4]
# condition_col = sys.argv[5]

dataset_path = "Visium_MTG_10samples"
kept_features =["nCount_Spatial","nFeature_Spatial","sample_name","sex","diagnosis","last_mmse_test_score","motor_updrs_score","smoothed_label_s5"]
sample_col = "sample_name"
cluster_col = "smoothed_label_s5"
condition_col = "diagnosis"

print("============================================")
print("Dataset path: ", dataset_path)
print("Kept features: ", kept_features)
print("Sample column: ", sample_col)
print("Cluster column: ", cluster_col)
print("Condition column: ", condition_col)

print("Checking inputs...")
if sample_col not in kept_features:
    kept_features.append(sample_col)
if cluster_col not in kept_features:
    kept_features.append(cluster_col)
if condition_col not in kept_features:
    kept_features.append(condition_col)


print("Loading metadata...")
metadata = pd.read_csv(dataset_path + "/raw_metadata.csv", index_col=0, header=0)
metadata = metadata.loc[:, kept_features]


metadata = metadata.rename(columns={condition_col: "Condition"})
kept_features.remove(condition_col)
kept_features.append("Condition")

## if the column data is float, keep 2 digits after the decimal point
# Round only float columns to 2 decimal places
metadata[metadata.select_dtypes(include=['float']).columns] = metadata.select_dtypes(include=['float']).round(2)

## check if "sample_id" column exists
if "sample_id" != sample_col:
    print("Renaming sample id...")
    metadata.drop("sample_id", axis=1, inplace=True, errors="ignore")
    metadata = metadata.rename(columns={sample_col: "sample_id"})
    kept_features.remove(sample_col)
    kept_features.append("sample_id")

## Rename the spot id as: SampleID_spotSerialNumber
print("Renaming spot id...")
new_ids = []
sample_cellspot_n = {}
for index, row in metadata.iterrows():
    sample_id = row["sample_id"]
    if sample_id not in sample_cellspot_n:
        sample_cellspot_n[sample_id] = 0
    sample_cellspot_n[sample_id] += 1
    c_id = sample_id + "_s" + str(sample_cellspot_n[sample_id])
    new_ids.append(c_id)
metadata["cs_id"] = new_ids

barcode_to_csid = metadata["cs_id"].to_dict()

metadata["barcode"] = metadata.index.tolist()
metadata = metadata.set_index("cs_id")

all_samples = metadata["sample_id"].unique().tolist()
with open(dataset_path + "/sample_list.json", "w") as f:
    json.dump(sorted(all_samples), f)

spot_to_sample = metadata["sample_id"].to_dict()
# Save spot_to_sample mapping
with open(f"{dataset_path}/cellspot_to_sample.json", "w") as f:
    json.dump(spot_to_sample, f, indent=2)

# %% ==============================================
## Process spot metadata
print("Processing spot metadata...")
## save original metadata
metadata.loc[:,kept_features].to_csv(dataset_path + "/cellspot_metadata_original.csv")

sample_level_features = []
spot_level_features = []
sample_groups = metadata.groupby("sample_id")
for feature in kept_features:
    is_sample_level = all(group[feature].nunique() <= 2 for _, group in sample_groups)
    if is_sample_level:
        sample_level_features.append(feature)
    else:
        spot_level_features.append(feature)

spot_meta_list = spot_level_features
metadata_lite = metadata.loc[:, spot_meta_list]


spot_meta_mapping = {}
for spot_meta in spot_meta_list:
    # Check if the column is categorical
    if is_categorical(metadata_lite[spot_meta], unique_threshold=0.2):
        # Convert to categorical
        cat_series = metadata_lite[spot_meta].astype("category")

        cat_counts = cat_series.value_counts().to_dict()

        # Replace original column with codes
        metadata_lite[spot_meta] = cat_series.cat.codes

        # Store mapping, and calculate the number of spots in each category
        mapping = {i: [cat, cat_counts[cat]] for i, cat in enumerate(cat_series.cat.categories)} ## with counts
        spot_meta_mapping[spot_meta] = mapping

# Save mapping to JSON
with open(dataset_path + "/cellspot_meta_mapping.json", "w") as f:
    f.write(dumps_compact_lists(spot_meta_mapping, indent=4))

metadata_lite.to_csv(dataset_path + "/cellspot_metadata.csv")

print("Processing sample metadata...")
sample_meta_list = sample_level_features
sample_meta = metadata.loc[:, sample_meta_list]
sample_meta = sample_meta.drop_duplicates()
sample_meta = sample_meta.set_index("sample_id")
sample_meta.fillna("", inplace=True)
sample_meta.to_csv(dataset_path + "/sample_metadata.csv")

with open(dataset_path + "/meta_list.json", "w") as f:
    json.dump(sorted(spot_meta_list + sample_meta_list), f)

# %% ==============================================
## Process embedding data
print("Loading embedding ....")
embeddings_data = pd.read_csv(dataset_path + "/raw_umap_embeddings.csv", index_col=0, header=0)

embeddings_data["UMAP_1"] = embeddings_data["UMAP_1"].round(2)
embeddings_data["UMAP_2"] = embeddings_data["UMAP_2"].round(2)

## reset index use barcode_cid map
# Reset index and rename using the mapping
print("Renaming embeddings....")
embeddings_data = embeddings_data.reset_index()  # Move index to a column
embeddings_data["index"] = embeddings_data["index"].map(barcode_to_csid)  # Rename using mapping
embeddings_data = embeddings_data.set_index("index")  # Set the renamed column as index
embeddings_data.to_csv(dataset_path + "/umap_embeddings.csv", index_label="cs_id")

## sampling umap, get 100k spots
print("Sampling umap...")
n_rows = embeddings_data.shape[0]

sample_rows = 100000 if n_rows > 100000 else n_rows
embeddings_data_nk = embeddings_data.sample(n=sample_rows, random_state=42)
embeddings_data_nk.to_csv(dataset_path + "/umap_embeddings_100k.csv", index_label="cs_id")

sample_rows = 50000 if n_rows > 50000 else n_rows
embeddings_data_nk = embeddings_data.sample(n=sample_rows, random_state=42)
embeddings_data_nk.to_csv(dataset_path + "/umap_embeddings_50k.csv", index_label="cs_id")

# %% ============================================================================
print("Processing coordinates...[renaming barcode to cs_id]")
files = os.listdir(dataset_path + "/coordinates")
for file in files:
    if file.endswith(".csv"):
        df = pd.read_csv(dataset_path + "/coordinates/" + file, index_col=0, header=0)
        df.rename(index=barcode_to_csid, inplace=True)
        df.to_csv(dataset_path + "/coordinates/" + file, index_label="cs_id")

# %% ============================================================================
## Process expression data
print("Loading expression data...(Takes a while...be patient...)")
## rename_expression_data
expression_data = pd.read_csv(dataset_path + "/raw_normalized_counts.csv", index_col=0, header=0)
## rename "Spot" column use barcode_cid map
print("Renaming expression....")
expression_data["cs_id"] = expression_data["Spot"].map(barcode_to_csid)
expression_data.drop("Spot", axis=1, inplace=True)

## "Expression" column keep 4 digits after the decimal point
expression_data["Expression"] = expression_data["Expression"].round(2)

# %% ============================================================================
## Save gene jsons
print("Saving gene jsons...")
expression_data["sample_id"] = expression_data["cs_id"].map(spot_to_sample)

# Group data by gene
print("Grouping by gene... be patient...")
grouped_by_gene = expression_data.groupby("Gene")

all_genes = grouped_by_gene.groups.keys()
all_genes = [gene_i.replace("/", "_") for gene_i in list(set(all_genes))]
with open(dataset_path + "/gene_list.json", "w") as f:
    json.dump(sorted(all_genes), f)

# Create directory for genes
os.makedirs(dataset_path + "/gene_jsons", exist_ok=True)

# Save each gene as a separate JSON file
print("Saving gene... be patient...")
i = 0
total_n = len(all_genes)
for gene, df in grouped_by_gene:
    try:
        i += 1
        if i % 1000 == 0:
            print(f"{i}/{total_n}")

        gene_dict = dict(zip(df["cs_id"], df["Expression"]))

        safe_gene_name = gene.replace("/", "_")

        # Create JSON file
        file_name = f"{dataset_path}/gene_jsons/{safe_gene_name}.json"
        with open(file_name, "w") as f:
            json.dump(gene_dict, f, indent=4)

    except Exception as e:
        print(f"Error in processing {gene} !!! Check the error_gene.txt")
        with open(dataset_path + "/error_gene_json.txt", "a") as f_err:
            f_err.write(gene + "\n")

# %% ============================================================================
## calculate psuedo count of each gene in each sample
print("=============================================")
print("Calculating pseudo count...")
os.makedirs(dataset_path + "/gene_pseudobulk", exist_ok=True)

print("Grouping by gene... be patient...")
# Compute pseudo-bulk counts by summing expression values per (sample, gene)
pseudo_bulk = expression_data.groupby(["sample_id", "Gene"])["Expression"].sum().reset_index()
# Rename the expression value column
pseudo_bulk.rename(columns={"Expression": "pseudobulk_expr"}, inplace=True)

# Save each gene's data to a separate JSON file
i = 0
for gene, df_gene in pseudo_bulk.groupby("Gene"):
    i += 1
    if i % 1000 == 0:
        print(f"{i}/{total_n}")
    gene_dict = df_gene.set_index("sample_id")["pseudobulk_expr"].to_dict()
    safe_gene_name = gene.replace("/", "_")
    with open(f"{dataset_path}/gene_pseudobulk/{safe_gene_name}.json", "w") as f:
        json.dump(gene_dict, f, indent=4)

print("Done! Feature/Gene data processed and saved.")

