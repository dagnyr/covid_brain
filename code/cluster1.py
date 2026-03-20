import os
import re
import scanpy as sc
import anndata as ad
import numpy as np
import pandas as pd

# ----------- settings
sc.settings.verbosity = 3
sc.settings.set_figure_params(dpi=120, facecolor="white")

DATA_DIR = "/scratch/users/dagnyr/brain/data_raw/choroid_plexus"          # folder containing all the GSM files
OUT_H5AD = "choroid_merged.h5ad"

# ----------- find + group samples
files = os.listdir(DATA_DIR)

prefix_re = re.compile(r"^(GSM\d+_[^_]+(?:_[^_]+)?)_(barcodes\.tsv|features\.tsv|matrix\.mtx)$")

samples: dict[str, dict] = {}
for f in files:
    m = prefix_re.match(f)
    if m:
        prefix, kind = m.group(1), m.group(2)
        samples.setdefault(prefix, {})[kind] = os.path.join(DATA_DIR, f)

print(f"Found {len(samples)} samples:\n  " + "\n  ".join(sorted(samples)))

# ----------- load each sample as adata object
def load_sample(prefix: str, paths: dict) -> ad.AnnData:

    import scipy.io as sio

    mat   = sio.mmread(paths["matrix.mtx"]).T.tocsr()
    
    barcodes = pd.read_csv(paths["barcodes.tsv"], header=None, sep="\t")[0].values
    features = pd.read_csv(paths["features.tsv"], header=None, sep="\t")

    gene_ids   = features[0].values # features file columns: gene_id, gene_name[, type]
    gene_names = features[1].values if features.shape[1] > 1 else gene_ids 

    adata = ad.AnnData(X=mat)
    adata.obs_names = [f"{prefix}_{bc}" for bc in barcodes]
    adata.var_names = gene_names
    adata.var["gene_ids"] = gene_ids
    adata.obs["sample"]   = prefix

    batch = re.sub(r"^GSM\d+_", "", prefix) # create unique batch id for each sample before merge
    adata.obs["batch"] = batch

    adata.var_names_make_unique()
    return adata

adatas = []
for prefix in sorted(samples):
    
    paths = samples[prefix]
    if len(paths) != 3:
        print(f"  ⚠  Skipping {prefix} — missing files: {paths}")
        continue
    adata = load_sample(prefix, paths)
    adatas.append(adata)

# ----------- concat samples into larger file
adata = ad.concat(adatas, join="outer", fill_value=0, label="sample", index_unique=None)
adata.obs_names_make_unique()

# quick print to verify merge worked
print(f"merged -> {adata.n_obs} cells × {adata.n_vars} genes")
print(adata.obs[["sample", "batch"]].drop_duplicates())

# ----------- compute qc metrics
adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

# ----------- filter cells

min_genes  = 200
#max_genes  = 10000
max_pct_mt = 20

print(f"before cell filter: {adata.n_obs}") # quick print to verify

sc.pp.filter_cells(adata, min_genes=min_genes)
adata = adata[adata.obs.n_genes_by_counts < max_genes].copy()
adata = adata[adata.obs.pct_counts_mt  < max_pct_mt ].copy()

print(f"after  cell filter: {adata.n_obs}")  # quick print to verify filtering step

# ----------- filter genes
print(f"before  gene filter: {adata.n_vars} genes") # quick print to verify filtering step
sc.pp.filter_genes(adata, min_cells=3)
print(f"after  gene filter: {adata.n_vars} genes") # quick print to verify filtering step

# ----------- normalize + log transform data
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

# ----------- find highly variable genes
sc.pp.highly_variable_genes(
    adata,
    n_top_genes=3000,
    batch_key="batch",
    flavor="seurat_v3",
    layer="counts",
)
print(f"No. HVG: {adata.var.highly_variable.sum()}") # quick print to check

adata = adata[:, adata.var.highly_variable].copy() # subset to HVG only

# ----------- scale
sc.pp.scale(adata, max_value=10)

# ----------- pca
sc.tl.pca(adata, svd_solver="arpack", n_comps=50)
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False, save="_variance_ratio.png")

# ----------- harmony
try:
    sc.external.pp.harmony_integrate(adata, key="batch", max_iter_harmony=20)
    use_rep = "X_pca_harmony" # use for KNN and UMAP if avail

except Exception as e:
    use_rep = "X_pca"

# ----------- KNN, UMAP, Leiden

sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30, use_rep=use_rep)
sc.tl.umap(adata)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_0.5")
sc.tl.leiden(adata, resolution=1.0, key_added="leiden_1.0")

# ----------- marker genes per cluster

sc.tl.rank_genes_groups(adata, groupby="leiden_0.5", method="wilcoxon", use_raw=True, key_added="rank_genes_leiden_0.5",)

adata.write_h5ad(OUT_H5AD, compression="gzip")
print("saved h5ad")
print(adata)

