import os
import gseapy as gp
import numpy as np
import pandas as pd
import scanpy as sc

adata = sc.read_h5ad("choroid_merged.h5ad")

out_dir = "/scratch/users/dagnyr/brain/results/choroid_plexus/gsea"
os.makedirs(out_dir, exist_ok=True)

gene_sets = [
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023",
    "KEGG_2021_Human",
    "Reactome_2022",
]

def run_gsea(adata_subset, label, out_dir):

    sc.tl.rank_genes_groups(
        adata_subset,
        groupby="disease",
        groups=["COVID-19"],
        reference="Control",
        method="wilcoxon",
        use_raw=True,
        key_added="rank_genes_disease",
    )

    de = sc.get.rank_genes_groups_df(
        adata_subset, group="COVID-19", key="rank_genes_disease"
    )
    de = de[~de["names"].str.startswith(("MT-", "MT"))]
    de = de.sort_values("scores", ascending=False)

    ranked = de.set_index("names")["scores"].dropna()

    results = {}

    for gs in gene_sets:
        try:
            pre_res = gp.prerank(
                rnk=ranked,
                gene_sets=gs,
                threads=4,
                min_size=15,
                max_size=500,
                permutation_num=1000,
                outdir=f"{out_dir}/{label}/{gs}",
                seed=42,
                verbose=False,
            )
            sig = pre_res.res2d[pre_res.res2d["FDR q-val"] < 0.05]

            results[gs] = pre_res.res2d

        except Exception as e:
            print(f"  {gs} failed: {e}")

    return results

# ----------- overall
overall_results = run_gsea(adata, "overall", out_dir)

for gs, df in overall_results.items():
    df_sorted = df.sort_values("FDR q-val")

    df_sorted.to_csv(f"{out_dir}/overall_{gs}.csv", index=False)

    print(f"\ntop 10 pathways — {gs}:")

# ----------- by cell type
for cell_type in adata.obs["cell_type"].unique():
    subset = adata[adata.obs["cell_type"] == cell_type].copy()

    if subset.n_obs < 50:
        continue

    if subset.obs["disease"].nunique() < 2:
        continue

    if subset.obs["disease"].value_counts().min() < 10:
        continue

    fname = cell_type.replace("/", "_").replace(" ", "_").replace("+", "plus")

    ct_results = run_gsea(subset, fname, out_dir)

    for gs, df in ct_results.items():
        df_sorted = df.sort_values("FDR q-val")

        df_sorted.to_csv(f"{out_dir}/{fname}_{gs}.csv", index=False)

