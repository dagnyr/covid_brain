#!/bin/bash
#SBATCH --job-name=harmony
#SBATCH --mem=64G
#SBATCH --time=04:00:00
#SBATCH --output=harmony.log

source ~/miniconda3/bin/activate covid_brain

python - << 'PYEOF'
import scanpy as sc
import harmonypy as hm
import numpy as np

adata = sc.read_h5ad("/scratch/users/dagnyr/lung_merged_filtered.h5ad")

sc.tl.pca(adata, svd_solver="arpack", n_comps=50)
print("PCA done.")

ho = hm.run_harmony(adata.obsm["X_pca"], adata.obs, "batch", max_iter_harmony=20)
adata.obsm["X_pca_harmony"] = np.array(ho.Z_corr)
use_rep = "X_pca_harmony"

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep=use_rep)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.3, key_added="leiden_0.3")
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_0.5")
sc.tl.leiden(adata, resolution=1.0, key_added="leiden_1.0")

adata.write_h5ad("/scratch/users/dagnyr/lung_merged_filtered.h5ad", compression="gzip")
PYEOF
