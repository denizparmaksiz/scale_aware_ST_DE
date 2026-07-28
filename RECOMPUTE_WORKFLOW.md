# Published and recompute workflow contract

## What the mode flag controls

`ALDEX_RESULT_MODE` selects the source of ALDEx-derived result tables used by a
notebook:

- `published` reads the pooled, frozen manuscript tables from Additional Files
  4–6. This is the recommended route for reproducing downstream figures and
  tables without refitting the models.
- `recompute` uses individual result workbooks produced by the separate ALDEx R
  workflow and, where applicable, regenerates the pooled Python outputs.

The flag does not mean that every notebook can reconstruct all of its inputs in
one uninterrupted execution. The analysis has three stages:

```text
Python notebook: prepare count/metadata inputs
                  |
                  v
R/ALDEx workflow: fit each model and write per-stratum workbooks
                  |
                  v
Python notebook: pool/postprocess results and create downstream analyses
```

The R step is documented under
`documentation_code/aldex_main_crossstudy_repo/`. It may be run locally or on
HPC, but its individual output workbooks must be copied to the notebook's
documented results directory before Python postprocessing continues.

## Why switching modes is asymmetric

A completed recompute workflow can always be followed by a downstream notebook
in `published` mode because the frozen pooled tables are independent of local
intermediates. The reverse is not generally true: running an upstream notebook
in `published` mode deliberately skips creation of local intermediates, so a
later notebook set to `recompute` may require an earlier notebook and/or the R
workflow to have been run first.

The mode may be selected independently in each notebook, but the user is
responsible for satisfying the recompute prerequisites listed below. The
notebooks now stop at each external boundary with an explanation instead of
passing an empty file list into postprocessing.

## Notebook-by-notebook contract

| Notebook | `published` source | What `recompute` does | Recompute prerequisites/boundary |
|---|---|---|---|
| 03 | Additional File 4 | Exports primary MERFISH count/metadata inputs and postprocesses primary ALDEx workbooks | Requires the deposited primary AnnData. After export, run the primary R workflows and place the individual result workbooks in the named `results/DEG/` subdirectories before continuing. |
| 04 | Additional File 4 | Reads locally postprocessed primary ALDEx tables | Requires the postprocessed cell-type-by-region workbook and gene pivot produced from the notebook 03/R workflow. It does not rerun notebook 03. |
| 05 | Additional File 4 | Reads a locally postprocessed primary ALDEx workbook and frozen annotation resources | Requires the recomputed workbook at the exact path printed by the preflight check. It does not fit ALDEx. |
| 06 | Additional Files 4–5 | Exports Sun binary-age inputs, postprocesses Sun workbooks, and recomputes rank-based comparisons | Requires the external Sun AnnData and deposited primary AnnData. The Sun R workflow must run after input export; primary recompute comparisons also require the relevant notebook 03 output. |
| 07 | Additional File 5 | Exports Sun continuous-age inputs, postprocesses individual continuous-age workbooks, and recalculates expression summaries | Requires the external Sun AnnData. Run the continuous cell-type R workflow after export, then continue the notebook. |
| 08 | Additional File 5 | Recomputes Sun subregion correlations/inputs and postprocesses continuous subregion ALDEx workbooks | Requires the external Sun AnnData and published Sun comparison tables. Run the continuous subregion R workflow after export, then continue the notebook. |
| 09 | Additional Files 4 and 6 | Exports CosMx inputs, postprocesses CosMx workbooks, and recomputes rank-based comparisons | Requires the deposited CosMx and primary MERFISH AnnData objects. Run the CosMx R workflow after export; the cross-platform ALDEx comparison also requires the relevant primary notebook 03 output. |

Notebook 02 consumes final simulation summaries from Additional File 3 and does
not use `ALDEX_RESULT_MODE`. Full simulation regeneration is a separate staged
workflow under `documentation_code/simulation_benchmarking/`.

## Historical BLMM sensitivity outputs

Complete per-stratum BLMM sensitivity workbooks are not redistributed. The
published engine comparison is in Additional File 4, sheet
`ALDEx_lme_vs_blmm_engine`. Recomputing that sensitivity analysis requires the
pinned development implementation described in
`documentation_code/aldex_main_crossstudy_repo/README.md`.
