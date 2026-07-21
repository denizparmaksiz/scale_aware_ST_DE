from __future__ import annotations

import os
import re
from pathlib import Path

import anndata as ad
import pandas as pd
import scanpy as sc
from scipy.io import mmread

IN_DIR = Path("/gpfs/Home/dpp5572/simulation_results/")
OUT_DIR = Path("/gpfs/Home/dpp5572/simulation_results_scanpy/")
OUT_DIR.mkdir(parents=True, exist_ok=True)

target_sim_raw = os.environ.get("TARGET_SIM")
target_sim = int(target_sim_raw) if target_sim_raw is not None else None
manifest_pattern = re.compile(r"^manifest_(.+)_sim_(\d{3})_N_(\d{3})\.rds$")
datasets: list[tuple[str, int, int]] = []

for path in IN_DIR.glob("manifest_*.rds"):
    match = manifest_pattern.match(path.name)
    if match is None:
        continue
    ct_region = match.group(1)
    sim_id = int(match.group(2))
    n_batches = int(match.group(3))
    if target_sim is None or sim_id == target_sim:
        datasets.append((ct_region, sim_id, n_batches))

datasets.sort(key=lambda item: (item[1], item[2], item[0]))
if not datasets:
    raise RuntimeError("No simulation datasets were found.")

for index, (ct_region, sim_id, n_batches) in enumerate(datasets, start=1):
    rep_id = f"{ct_region}_sim_{sim_id:03d}_N_{n_batches:03d}"
    print(f"[{index}/{len(datasets)}] {rep_id}", flush=True)
    counts = mmread(IN_DIR / f"counts_{rep_id}_observed.mtx").T.tocsr()
    obs = pd.read_csv(IN_DIR / f"obs_{rep_id}.csv")
    var = pd.read_csv(IN_DIR / f"var_{rep_id}.csv")
    adata = ad.AnnData(X=counts, obs=obs.copy(), var=var.copy())
    adata.obs_names = adata.obs["cell_id"].astype(str)
    adata.var_names = adata.var["gene"].astype(str)
    adata.var["keep_gene"] = adata.var["keep_gene"].astype(bool)
    adata = adata[:, adata.var["keep_gene"]].copy()
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.tl.rank_genes_groups(
        adata,
        groupby="age_perm",
        method="wilcoxon",
        reference="Yng",
    )
    results = sc.get.rank_genes_groups_df(adata, group="Old")
    truth = adata.var.loc[:, ["gene", "is_de", "true_log2fc"]].reset_index(drop=True)
    results = (
        results.rename(
            columns={
                "names": "gene",
                "scores": "age_perm:score",
                "logfoldchanges": "age_perm:logfc",
                "pvals": "age_perm:pval",
                "pvals_adj": "age_perm:pval.adj",
            }
        )
        .merge(truth, on="gene", how="left")
    )
    results.to_csv(OUT_DIR / f"scanpy_{rep_id}_observed.csv", index=False)


# Combine all currently available per-dataset Scanpy results.
result_pattern = re.compile(
    r"^scanpy_(.+)_sim_(\d+)_N_(\d+)_observed\.csv$"
)
combined = []

for path in OUT_DIR.glob("scanpy_*_observed.csv"):
    match = result_pattern.match(path.name)
    if match is None:
        continue

    result = pd.read_csv(path)
    result["ct_region"] = match.group(1)
    result["sim_id"] = int(match.group(2))
    result["N_batches"] = int(match.group(3))
    result["condition"] = "observed"
    result["source_file"] = path.name
    combined.append(result)

if combined:
    pd.concat(combined, ignore_index=True).to_csv(
        OUT_DIR / "scanpy_simulation_results_observed.csv",
        index=False,
    )
