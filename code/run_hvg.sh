#!/bin/bash
#SBATCH --job-name=hvg
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=hvg.log

source ~/miniconda3/bin/activate covid_brain

python - << 'PYEOF'
import scanpy as sc
import subprocess
import sys

adata = sc.read_h5ad("/scratch/users/dagnyr/lung_merged_filtered.h5ad")

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    batch_key="batch",
    flavor="seurat_v3",
    layer="counts",
)

genes_of_interest = ["CX3CL1", "CX3CR1", "ADAM10", "ADAM17", "MMP2"]

for gene in genes_of_interest:
    if gene in adata.var_names:
        if not adata.var.loc[gene, "highly_variable"]:
            adata.var.loc[gene, "highly_variable"] = True
        else:
            print(f"{gene} already in list")
    else:
        print(f"{gene} not found in var_names at all")

adata = adata[:, adata.var.highly_variable].copy()

sc.pp.scale(adata, max_value=10)

adata.write_h5ad("/scratch/users/dagnyr/lung_merged_filtered.h5ad", compression="gzip")

PYEOF
