#!/bin/bash
#SBATCH --job-name=normalize
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --output=normalize.log

source ~/miniconda3/bin/activate covid_brain

python - << 'PYEOF'
import scanpy as sc

adata = sc.read_h5ad("/scratch/users/dagnyr/lung_merged_filtered.h5ad")
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata
adata.write_h5ad("/scratch/users/dagnyr/lung_merged_filtered.h5ad", compression="gzip")

PYEOF
