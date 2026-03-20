import scanpy as sc
import gseapy as gp
import pandas as pd
import numpy as np
import os

adata = sc.read_h5ad("cortex_merged.h5ad")

out_dir = "/scratch/users/dagnyr/brain/results/cortex/gsea"
os.makedirs(out_dir, exist_ok=True)

# gene sets to test against — good options for brain/immune data
gene_sets = [
    "GO_Biological_Process_2023",
    "GO_Molecular_Function_2023", 
    "KEGG_2021_Human",
    "Reactome_2022",
]

def run_gsea(adata_subset, label, out_dir):
    """
    Run GSEA comparing COVID-19 vs Control using preranked approach.
    Genes are ranked by log2FC * -log10(pval) from wilcoxon test.
    """
    print(f"\n── GSEA: {label} ──")
    
    # run differential expression COVID vs Control
    sc.tl.rank_genes_groups(
        adata_subset,
        groupby="disease",
        groups=["COVID-19"],
        reference="Control",
        method="wilcoxon",
        use_raw=True,
        key_added="rank_genes_disease",
    )
    
    # extract results into a ranked gene list
    de = sc.get.rank_genes_groups_df(
        adata_subset,
        group="COVID-19",
        key="rank_genes_disease",
    )
    
    # remove MT genes
    de = de[~de["names"].str.startswith(("MT-", "MT"))]
    
    # rank by score
    de = de.sort_values("scores", ascending=False)
    ranked = de.set_index("names")["scores"].dropna()
    
    print(f"  Ranked {len(ranked)} genes")
    
    # run preranked GSEA
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
            print(f"  {gs}: {len(sig)} significant pathways")
            results[gs] = pre_res.res2d
        except Exception as e:
            print(f"  {gs} failed: {e}")
    
    return results


# ── Overall GSEA ──────────────────────────────────────────────────────────────
print("=== Overall GSEA ===")
overall_results = run_gsea(adata, "overall", out_dir)

# save overall results
for gs, df in overall_results.items():
    df_sorted = df.sort_values("FDR q-val")
    df_sorted.to_csv(f"{out_dir}/overall_{gs}.csv", index=False)
    print(f"\nTop 10 pathways — {gs}:")
    print(df_sorted[["Term", "ES", "NES", "NOM p-val", "FDR q-val"]].head(10).to_string(index=False))

# ── Per cell type GSEA ────────────────────────────────────────────────────────
print("\n=== Per Cell Type GSEA ===")
for cell_type in adata.obs["cell_type"].unique():
    subset = adata[adata.obs["cell_type"] == cell_type].copy()
    
    # skip if too few cells or only one disease group
    if subset.n_obs < 50:
        print(f"Skipping {cell_type} — too few cells ({subset.n_obs})")
        continue
    if subset.obs["disease"].nunique() < 2:
        print(f"Skipping {cell_type} — only one disease group")
        continue
    if subset.obs["disease"].value_counts().min() < 10:
        print(f"Skipping {cell_type} — too few cells in one group")
        continue

    fname = cell_type.replace("/", "_").replace(" ", "_").replace("+", "plus")
    ct_results = run_gsea(subset, fname, out_dir)
    
    # save per cell type results
    for gs, df in ct_results.items():
        df_sorted = df.sort_values("FDR q-val")
        df_sorted.to_csv(f"{out_dir}/{fname}_{gs}.csv", index=False)
