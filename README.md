# SpatialProteomics

SpatialProteomics documents an end-to-end workflow for tandem mass tag (TMT) fractionation studies that detect covarying protein networks in wild-type (WT) and Parkinson’s disease (MUT) mouse models. The repository bundles the scripts, notebooks, and reference data needed to denoise proteomic measurements, cluster proteins with Leiden community detection, and quantify WT–MUT contrasts with linear mixed models.

## Background

Each biological replicate is fractionated into seven cellular compartments, labeled with TMT barcodes, and analyzed via mass spectrometry. Intensities are normalized, filtered, and converted into adjacency matrices describing protein co-fractionation. After network enhancement, Leiden clustering identifies modules that capture shared behavior across fractions. Linear mixed models then test whether the MUT condition shifts protein abundance within each module. Significant modules (p ≤ 0.05 after FDR correction) are summarized, visualized, and paired with pathway-level enrichment analyses.

## Repository layout

| Path | Purpose |
| --- | --- |
| `src/` | Primary analysis pipeline; subdirectories `1_`–`7_` mirror the processing stages (QC, normalization, clustering, module modeling, plotting, Shiny assets). |
| `postprocess/` | WT vs MUT plotting notebooks, KEGG lookups, and publication-ready figure templates (formerly `compare_MutvsWTfractions/`). |
| `data/`, `rawdata/`, `transformeddata/` | Representative input files, intermediate tidy tables, and transformed adjacency matrices. |
| `tables/` | Aggregated Excel/CSV outputs (module memberships, MSstats results, enrichment summaries). |
| `example_data/` | Lightweight datasets covering every pipeline stage; useful for smoke tests or training. |
| `ShinyGo/`, `geneontologies/` | Gene-ontology helpers including ShinyGO automation notebooks and custom KEGG-based scripts. |
| `rdata/`, `rdata.tar` | Serialized R objects for quick reuse across machines. |

## Example data

Any script under `src/` can be driven by the contents of `example_data/`:

1. Copy the matching example files into `data/` (or point scripts directly to `example_data/`).
2. Run the stages listed below to confirm the environment is configured before using full datasets.

## Workflow overview

1. **Acquisition QC** – `src/1_WASH-iBioID` scripts inspect BioID runs, annotate accessions, and create tidy inputs.
2. **TMT normalization** – `src/2_SWIP-TMT` executes MSstatsTMT normalization, removes control channels, and exports contrasts.
3. **Network building** – `src/3_Clustering` generates adjacency matrices, applies network enhancement, and clusters proteins with Leiden (`2_leidenalg-clustering.py`).
4. **Module-level modeling** – `src/4_Module-Analysis/1_module-lmerTest-analysis.R` fits linear mixed models per module and exports `{PROJECT}_TMT-Module-Results.xlsx`.
5. **Visualization** – `src/5_Plotting` scripts build protein, module, PCA, and variance figures; `src/7_Shiny` holds the `module-profiles` Shiny application.
6. **Post-processing** – `postprocess/` notebooks aggregate replicates, highlight WT vs MUT trajectories, and prepare pathway summaries.

## Environment setup

### R (renv)

The R toolchain is pinned with `renv`.

1. Launch R in the repository root.
2. Set the working directory if necessary: `setwd('path/to/SpatialProteomics')`.
3. Install renv (first run only): `install.packages("renv")`.
4. Restore the environment tracked in `src/renv.lock`: `renv::restore(lockfile = "src/renv.lock")`.
5. Re-run `renv::restore()` whenever dependencies change. Missing packages can be installed directly or via `BiocManager`, as noted in script headers.

`src/install_packages.R` lists HPC-specific dependencies that occasionally need manual installation. Container recipes (`src/Dockerfile`, `src/rocker-r.sif`) capture a fully reproducible setup.

### Python

Python notebooks rely on pandas, numpy, scipy, seaborn, matplotlib, and gseapy. Create an environment from `src/environment_spatial.yml` (conda or mamba) or manage packages manually if preferred.

## Running the pipeline

### 1. Protein detection and QC (`src/1_WASH-iBioID`)

`1_lmTest-BioID-analysis.R` cleans BioID runs, enforces species filters, and exports QC tables. Inputs default to `data/`; adjust the paths at the top of the script for custom datasets. Outputs are written to `tables/`.

### 2. TMT normalization and statistical testing (`src/2_SWIP-TMT`)

`1_MSstatsTMT-analysis.R` (with accompanying notebooks) performs log2 transforms, normalization, and MSstatsTMT contrasts; results populate `transformeddata/` and `tables/`.

### 3. Network and clustering (`src/3_Clustering`)

`1_generate-network.R` constructs fraction-wise adjacency matrices. `2_leidenalg-clustering.py` detects communities, optionally leveraging the `neten` enhancement described in the references. Supporting notebooks such as `3.5_mapaccesions.ipynb` map IDs and merge metadata.

### 4. Module analysis (`src/4_Module-Analysis`)

`1_module-lmerTest-analysis.R` fits linear mixed models that contrast WT and MUT within every module. The exported `{PROJECT}_TMT-Module-Results.xlsx` files serve as inputs for downstream enrichment and plotting. `2_module-GSEA.R` runs GSEA-based enrichment on the same modules.

### 5. Visualization (`src/5_Plotting` and `src/7_Shiny`)

Scripts `00_generate-colors.R` through `10_plot-variance.R` produce publication-quality figures. They expect the module assignments generated in step 4 and emit PNG/PDF outputs into `tables/` (or a path you configure). The `7_Shiny` directory houses the `module-profiles` application for interactive exploration.

## Post-processing and reporting (`postprocess/`)

`postprocess/` aggregates utilities needed after the core pipeline:

- `WT_Mut_fractionPlot.ipynb` overlays WT and MUT replicates for each protein/module with both individual traces and averaged profiles.
- `LRRK2_MixtureM1/` and `LRRK2_MixturesAveraged/` demonstrate how replicate summaries feed into figures.
- `environment_plot.yml` defines a minimal conda environment for the plotting notebooks.
- `get_KEGGnums.ipynb` and `get_pathways.ipynb` fetch KEGG ortholog IDs and execute custom hypergeometric enrichment when a focused pathway view is required.

## Gene ontology resources

Choose the enrichment strategy that best matches your downstream analysis:

1. **Scripted GSEA** – `src/4_Module-Analysis/2_module-GSEA.R`.
2. **Automated ShinyGO scraping** – `ShinyGo/ShinyGO_analysis/module_pathways_1.ipynb` plus `cluster_analysis_2.py`.
3. **Custom KEGG workflows** – `geneontologies/get_KEGGnums.ipynb` followed by `geneontologies/get_pathways.ipynb`, which emit CSVs to `enrichments_DIY/` and `enrichments_GSEA/`.

Recommended usage:

1. Start from `{GENE}_TMT-Module-Results.xlsx` in `tables/`.
2. Supply module-to-gene mappings to the enrichment notebook or script of choice.
3. Interpret modules with adjusted p-values ≤ 0.05 and feed significant results to the plotting notebooks in `postprocess/`.

## Plotting WT vs MUT fractions

`postprocess/WT_Mut_fractionPlot.ipynb` generates the WT vs MUT fraction plots used across presentations. Update the configuration cell (input table, module of interest, color palette) and run the notebook sequentially. Helper functions at the end can be uncommented for targeted debugging.

## References

1. Courtland, J.L., Bradshaw, T.W.A., Waitt, G., Soderblom, E., Ho, T., Rajab, A., Vancini, R., Kim, I.H., & Soderling, S.H. (2021). Genetic disruption of WASHC4 drives endo-lysosomal dysfunction and cognitive-movement impairments in mice and humans. *eLife*, 10:e61590. [doi:10.7554/eLife.61590](https://elifesciences.org/articles/61590)
2. Traag, V.A., Waltman, L., & Van Eck, N.-J. (2019). From Louvain to Leiden: guaranteeing well-connected communities. *Scientific Reports*, 9, 5233. [doi:10.1038/s41598-019-41695-z](https://www.nature.com/articles/s41598-019-41695-z)
3. Wang, B., Pourshafeie, A., Zitnik, M., Zhu, J., Bustamante, C.D., Batzoglou, S., & Leskovec, J. (2018). Network enhancement as a general method to denoise weighted biological networks. *Nature Communications*, 9, 3108. [doi:10.1038/s41467-018-05469-x](https://www.nature.com/articles/s41467-018-05469-x)
4. Ge, S.X., Jung, D., & Yao, R. (2020). ShinyGO: a graphical gene-set enrichment tool for animals and plants. *Bioinformatics*, 36(8), 2628–2629. [doi:10.1093/bioinformatics/btz931](https://doi.org/10.1093/bioinformatics/btz931)
5. Fang, Z., Liu, X., & Peltz, G. (2022). GSEApy: a comprehensive package for performing gene set enrichment analysis in Python. *Bioinformatics*, btac757. [doi:10.1093/bioinformatics/btac757](https://doi.org/10.1093/bioinformatics/btac757)
