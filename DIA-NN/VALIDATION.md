# Validation of the public workflow

Validated on 16 November 2026 using the archived 45-run CoPPRA dataset, R 4.4.2, and the package versions in `reference_package_versions.csv`. Inputs and archived results were read only; regenerated outputs were written to a separate directory.

- Python, R, R Markdown, and shell syntax checks passed. The runner's help command worked and an unknown stage was rejected.
- The complete shell runner finished all 14 analysis stages with `COPPRA_PERMUTATIONS=100`. The published default remains 10,000; the full 10,000-permutation run was not repeated for this code-publication check.
- All six input-file SHA-256 hashes matched the archived manifest.
- All 67 QC, PCA, heatmap, and limma CSV tables were byte-identical to the archived outputs. All 18 PTM enrichment CSV tables and 21 other GO/enrichment CSV tables were also byte-identical.
- Twelve QC handoff components matched numerically at tolerance 1e-12, including sample assignments, precursor/biological matrices, primary peptide eligibility, annotations, and method settings. Input path strings change with the new directory layout.
- The first 100 stored joint permutation assignments and responsive-peptide counts matched the archived run exactly.
- Built-in validation passed for filter order, normalization, preservation of observed values, peptide-median imputation, limma response calls, and the RU-only protein definition. The primary matrix contained 42,028 peptides across 12 biological samples, with 2,276 imputed values.
- Plot generation completed, the 3D PCA renderer reported no sample-label overlaps, and the renamed heatmap title was visually inspected.

The three KEGG GSEA tables differed because the live query supplied a different tested pathway set: `hsa04900` was present and `hsa04918` was absent relative to the archived results. Statistics for shared pathways were unchanged except for some RU adjusted P/q values following the altered multiple-testing inputs. BH-significant pathway sets were unchanged in all three contrasts. The analysis sources and thresholds were preserved; the workflow records the gene sets used by each run.

Existing enrichment warnings remain visible in per-stage warning logs: some PTM categories have undefined enrichment statistics, a small fraction of gene identifiers do not map, and rank ties affect approximation. These were not reinterpreted as significant findings. Additional package-build and plotting deprecation warnings were nonfatal.

This validation covers the DIA-NN workflow.
