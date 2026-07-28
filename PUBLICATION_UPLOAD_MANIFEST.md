# Authoritative publication upload manifest

Zenodo record: [10.5281/zenodo.21420404](https://doi.org/10.5281/zenodo.21420404)

Audit date: 2026-07-20

This manifest gives explicit destinations for the repository in its current
analysis form. It supersedes the provisional upload language in earlier planning
notes. Decisions are based on direct notebook/script consumption, verified file
contents, and whether the data originate from this study or an external source.

## 1. Upload to Zenodo

The current analysis requires the following 35 files, totaling approximately
9.48 GiB. Preserve the directory/group labels below in the Zenodo README or in
a single archive so notebook configuration remains understandable.

### A. Primary and derived analysis-ready objects

Upload all six:

| Zenodo group | Current file | Size | Why required |
|---|---|---:|---|
| `primary_merfish/` | `adata_glia_aldex_pp.h5ad` | 1.60 GB | Direct input to notebooks 03, 04, 06, and 09; source for primary ALDEx input export and rank-based comparisons |
| `primary_merfish/` | `adata_ct_centered_dec25.hdf5` | 3.25 GB | Direct input to notebook 04; no producer for this exact object exists in the repository |
| `figure1/` | `adata_500_forfig1.hdf5` | 1.07 GB | Direct input to notebook 01 |
| `figure1/` | `adata_140_forfig1.hdf5` | 181.75 MB | Direct 140-gene input to notebook 01 |
| `figure1/` | `adata_1000_forfig1.hdf5` | 712.25 MB | Direct 1000-gene input to notebook 01 |
| `cosmx/` | `cosmx_data_ct_final.h5ad` | 1.42 GB | Direct input to notebooks 09 and 10; preserves the exact consensus labels produced by the multi-stage MapMyCells/scANVI workflow |

Clear consequence: if the full 140- or 1000-gene objects cannot be released,
they cannot simply be omitted. Before publication, notebook 01 must instead be
rewritten to consume deposited figure-specific derived tables, and those tables
must be generated and validated against the full objects. Until that replacement
exists, both full objects are required Zenodo files.

### B. Raw Vizgen inputs for spatial panels

Upload the three files listed below from each of the five unique section folders.
Do not upload `detected_transcripts.csv`, image pyramids, result tiles, `.vzg`
archives, segmentation parquet files, or other folder contents because the
numbered notebooks do not read them.

Required within every listed section:

- `cell_by_gene.csv`
- `cell_metadata.csv`
- `images/micron_to_mosaic_pixel_transform.csv`

Required sections:

1. `young_anterior_section_example_1` — notebooks 01 and 04; required subset 193.34 MB.
2. `old_anterior_section_example_1` — notebook 04; required subset 143.39 MB.
3. `young_anterior_section_example_2` — notebook 01; required subset 132.04 MB.
4. `old_posterior_section_example_1` — notebook 04; required subset 185.30 MB.
5. `young_posterior_section_example_1` — notebook 04; required subset 200.91 MB.

Also upload these five study-generated CCF annotation tables used for the
representative spatial figures in notebooks 01 and 04:

- `cell_metadata_wCCF_regions_old_anterior_section_example_1.csv` — 23,861,200 bytes.
- `cell_metadata_wCCF_regions_old_posterior_section_example_1.csv` — 31,951,500 bytes.
- `cell_metadata_wCCF_regions_young_anterior_section_example_1.csv` — 31,947,073 bytes.
- `cell_metadata_wCCF_regions_young_anterior_section_example_2.csv` — 22,463,802 bytes.
- `cell_metadata_wCCF_regions_young_posterior_section_example_1.csv` — 33,704,080 bytes.

### C. Simulation source inputs

Upload all three:

- `adata_nn_raw_counts_transposed.csv` — 240.76 MB.
- `adata_nn_raw_counts_transposed_MD.xlsx` — 41.03 MB.
- `updated_sim_design.rds` — 7.78 KB.

These are direct inputs to `01_generate_simulations.R`. The RDS must be retained
even though Additional File 3 contains a human-readable `Simulation_design`
sheet, because the released simulation script reads the R object directly.

Do **not** upload per-replicate simulated matrices, per-method intermediate
outputs, SLURM logs, or standalone combined result workbooks. The source inputs
and stored design/seeds support full recomputation; Additional File 3 supports
the fast figure-reproduction route.

### D. Final manuscript Additional Files

Upload the final versions of all six to Zenodo as well as submitting them to the
journal:

- `Additional_File_1.xlsx`
- `Additional_File_2.pdf`
- `Additional_File_3.xlsx`
- `Additional_File_4.xlsx`
- `Additional_File_5.xlsx`
- `Additional_File_6.xlsx`

They total approximately 21.43 MB and provide the stable published-results route
for simulation, primary MERFISH, Sun, and CosMx downstream analyses. They should
not be committed to GitHub.

## 2. Upload to GitHub

### A. Main notebooks

Upload exactly the ten numbered notebooks after path cleanup and testing:

1. `01_fig1_panel_composition.ipynb`
2. `02_simulation_benchmarking_plots.ipynb`
3. `03_main_merfish_aldex_prepost.ipynb`
4. `04_main_MERFISH_bio_figs.ipynb`
5. `05_primaryMERFISH_functional_annotations.ipynb`
6. `06_Sun_main_binary_DE.ipynb`
7. `07_sun_continuous_raw_validation.ipynb`
8. `08_sun_supplementary_analyses.ipynb`
9. `09_cosmx_figs_DE.ipynb`
10. `10_cosmx_merfish_module_analysis.ipynb`

Do not upload `aging_studies_venn.ipynb`; it is a superseded source notebook now
integrated into notebook 01.

### B. Analysis scripts and manifests

Upload:

- `documentation_code/simulation_benchmarking/R/`
- `documentation_code/simulation_benchmarking/python/`
- `documentation_code/simulation_benchmarking/hpc/`
- `documentation_code/simulation_benchmarking/environments/`
- `documentation_code/aldex_main_crossstudy_repo/R/`
- `documentation_code/aldex_main_crossstudy_repo/hpc/`
- `documentation_code/aldex_main_crossstudy_repo/config/`
- `documentation_code/aldex_main_crossstudy_repo/environments/`
- `documentation_code/aldex_main_crossstudy_repo/tutorial/`
- `documentation_code/cosmx_python_scripts/`
- `documentation_code/go_enrichment/`

### C. Small deterministic resources

Move these out of `to_upload/` into a documented `resources/` hierarchy and
commit them to GitHub:

- `CosMx Data.xlsx`
- `cosmx_cell_type_markers.json`
- `cosmx_neuron_markers.json`
- `mygene_annotations.xlsx`
- `filtered_mygene_annotations_crossstudy.xlsx`
- `module_dictionary.json`
- `cosmx_module_dictionary.json`
- `updated_merfish_module_dictionary.json`
- `shared_crossstudy_module_dictionary.json`
- `recovered_updated_merfish_assignments.csv`
- `unassigned_updated_merfish_terms.csv`
- `tutorial_counts.csv`
- `tutorial_metadata.csv`

Also commit `MERFISH_Master-List-03_SM.xlsx` under a clearly named primary
MERFISH resource folder. It is only 25,452 bytes and is directly read during
primary ALDEx postprocessing.

### D. Environments, shared code, tests, and documentation

Upload:

- The historical provenance exports and the portable core/optional environment
  specifications in `environment_exports/`.
- `resources/allen_ccf/CCFv3OntologyStructure_u16.xlsx`, the exact small
  ontology workbook read by notebook 04. Do not substitute
  `atlas_info_KimRef_v2_segmentation.xlsx`.
- The final minimal `src/scale_aware_st/` package.
- The tested minimal subset of `kimlabspatial`; do not publish the current full
  historical directory unchanged.
- `tests/` without its `__pycache__` directory.
- `pyproject.toml`.
- Main README, data-download guide, notebook execution guide, Zenodo manifest,
  license, citation metadata, and contribution/contact information.

## 3. External download only — do not upload copies

Document exact source links, citations, expected filenames, versions/releases,
and preparation steps for:

- Sun et al. `aging_coronal.h5ad` source data.
- Sun et al. Supplementary Tables S7, S8, and S11.
- The three Tabula Muris Senis, Ximerakis et al., and Jin et al. workbooks used
  in Figure 1D.
- Five 10x Xenium `outs` directories obtained from three dataset pages and used
  by notebook 01: three fresh-frozen replicate directories, one wild-type
  time-course directory (the notebook variable is currently named
  `xenium_350_dir`), and one 5K Xenium Prime directory.
- SenNet `SNT638.MWRV.378` raw CosMx expression and metadata files.
- Allen MapMyCells marker JSON and precomputed-statistics HDF5.
- Allen ABC Atlas metadata and WMB-10Xv2 expression partitions used by scANVI.

These are third-party inputs. The repository should validate expected files and
provide download instructions; it should not silently substitute locally renamed
copies.

## 4. Do not upload anywhere

- `data_for_deeper_test/` as a folder. Individual files from it are assigned to
  Zenodo above; the folder itself is local staging.
- `Supplementary tables/old/`.
- `repo_planning.md` after the public README and manifest are complete.
- `repository_inventory.json` and `asset_inventory.json` after final QA.
- `__pycache__/`, `.pyc`, `.ipynb_checkpoints/`, local caches, temporary files,
  intermediate figures, debug logs, memory logs, and trained model caches.
- Full Allen reference downloads, MapMyCells extended HDF5 outputs, scANVI
  caches/models, and regenerated public raw datasets.
- The full contents of the five Vizgen result folders beyond the 15 explicitly
  listed input files.
- Individual ALDEx workbooks when their pooled rows are already present in
  Additional Files 4–6.
- Per-stratum ALDEx count/metadata exports: these are regenerated by notebooks
  03, 06, 07/08, and 09 from the deposited analysis-ready object or the external
  Sun object.
- Per-replicate simulation data and method intermediates.
- `HYBRID_benchmark_ready.xlsx`.
- `per_gene_results_with_effect_size_correlations.xlsx`.

## 5. Verified simulation-workbook equivalence

The following comparisons were performed cell-by-cell, treating blank cells
consistently and using a numeric tolerance of 1e-12:

| Standalone source | Additional File 3 sheet | Shape compared | Differences |
|---|---|---:|---:|
| `HYBRID_benchmark_ready.xlsx:signal_plot_ready` after removing constant `ct_region` | `Performance_summary` | 901 × 13 | 0 |
| `HYBRID_benchmark_ready.xlsx:signal_avg_by_N` | `Performance_summary_by_N` | 31 × 8 | 0 |
| `per_gene_results_with_effect_size_correlations.xlsx:per_gene` | `Gene_level_results` | 119,341 × 10 | 0 |
| `per_gene_results_with_effect_size_correlations.xlsx:effect_size_correlations` | `Effect_size_correlations` | 901 × 10 | 0 |

The removed `ct_region` column contained only `Astro_Isocortex` in all 900 data
rows. The other method-specific sheets in the standalone workbooks are component
views of the combined tables and are not read by notebook 02. Therefore the two
standalone workbooks are not publication inputs.

## 6. Notebook-by-notebook source and upload routing

| Notebook | Zenodo inputs | GitHub inputs | External-download inputs | Additional File route |
|---|---|---|---|---|
| 01 | Three panel HDF5s; five section input triplets; study-generated CCF metadata CSVs for its displayed sections | Figure 1D loader/config code | Xenium, raw SenNet CosMx, three published aging workbooks | None |
| 02 | Simulation source inputs only for full recomputation | Plotting notebook and simulation scripts | None | Additional File 3 for all plotting tables |
| 03 | `adata_glia_aldex_pp.h5ad` | MERFISH master list, ALDEx scripts/manifest | None | Additional File 4 for pooled results |
| 04 | Two primary MERFISH objects; four section input triplets; study-generated CCF metadata CSVs for its displayed sections | Plotting code and `CCFv3OntologyStructure_u16.xlsx` | None | Additional File 4 for primary results |
| 05 | None beyond Additional File 4 | Frozen annotations and module dictionary | Optional live MyGene service only | Additional File 4 |
| 06 | Primary MERFISH H5AD | Comparison code | Sun AnnData | Additional Files 4 and 5 |
| 07 | None | Validation code | Sun AnnData and Table S7 | Additional File 5 |
| 08 | None | Validation and GO code | Sun AnnData and Tables S7/S8/S11 | Additional File 5 |
| 09 | Final CosMx and primary MERFISH H5ADs | Marker JSONs and comparison code | Raw SenNet data only for reconstruction | Additional Files 4 and 6 |
| 10 | Final CosMx H5AD | Frozen annotations and four module dictionaries | Optional live MyGene service only | Additional Files 4 and 6 |

## 7. Remaining work before uploading

The file destinations above are final for the current analysis contract. The
following work changes code quality and usability but does not postpone the
collection of the listed Zenodo files:

1. Replace absolute paths with repository configuration.
2. Make notebook 02 read Additional File 3 directly.
3. Add explicit recompute-versus-published-result modes to ALDEx consumers.
4. Compare nonidentical repeated functions semantically before consolidation.
5. Reduce and test `kimlabspatial`.
6. Execute notebook smoke and end-to-end tests.
7. Produce checksums and a final Zenodo archive layout.
