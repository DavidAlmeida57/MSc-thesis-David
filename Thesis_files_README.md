---
editor_options: 
  markdown: 
    wrap: 72
---

# Thesis_files — README

## Purpose

This repository contains the code and supporting data used for the
proteomics analyses in the thesis. It includes scripts to preprocess
Proteome Discoverer results, filter proteins, compute fold-changes and
abundance ratios, assign proteins to biological groups (Independent,
Phosphorylated, Dephosphorylated), compare detected interactors with
external PPI resources, plot enriched pathways, and perform PCA on
filtered proteins.

------------------------------------------------------------------------

## Repository structure (as provided)

```         
Thesis_files/
├── Data/                   # All input tables used by the scripts (Proteome Discoverer exports, BioGRID, APID, G:Profiler sheets, selected tables)
└── Scripts/
    ├── thesis_analysis_script.R
    └── pathway_plots_script.R
```

------------------------------------------------------------------------

## Contents & brief description

### `thesis_analysis_script.R`

Main analysis pipeline that: - Reads the Proteome Discoverer export
(expected in `Data/`). - Computes group-level mean abundances per
condition/day. - Applies quality filters (Master protein, contaminants,
FDR, unique peptides, coverage). - Computes adjusted fold-changes (FC)
vs `N1` and filters proteins with `FC > 1.5`. - Computes SE/SA and SA/SE
ratios and uses these to split intersections when relevant. - Builds
protein groups for Day 3 and Day 7: Independent (I), Phosphorylated (P),
Dephosphorylated (D). - Compares detected interactors with BioGRID and
APID lists. - Produces summary tables, Venn diagrams and barplots. -
Writes intermediate and final Excel outputs.

**Main outputs written (default names used in script)**: -
`Dados Proteinas antes da filtragem.xlsx` — intermediate table after
computing group means - `Lista de Proteinas(I,P,D)_FC_1,5.xlsx` — final
lists of Accessions (I, P, D) by day - Optionally: additional
`write_xlsx()` outputs if enabled in script

------------------------------------------------------------------------

### `pathway_plots_script.R`

Pathway plotting utilities that: - Load G:Profiler pathway enrichment
results exported to Excel (KEGG/WikiPathways/Reactome sheets). - Convert
numeric columns and optionally compute unified color/size limits across
Day 3 / Day 7 pairs. - Produce publication-ready horizontal segment+dot
plots where: - x = `-log10(adjusted p-value)` - point size =
intersection size (\# genes) - point color = adjusted p-value (FDR) -
Creates plots for exclusive and common pathways for SA, SE and Total,
and for curated lists (top 20). - Optionally you can enable `ggsave()`
lines to export PNGs.

------------------------------------------------------------------------

### PCA analysis (integrated in `thesis_analysis_script.R`)

The script includes a PCA module that: - Selects proteins filtered by FC
threshold (FC \> 1.5) for Day 3 and Day 7 and also a combined Day3+Day7
analysis. - Builds abundance matrices from the raw `F` columns
corresponding to each sample, then: - applies `log2(x + 1)` transform, -
transposes to have samples as rows and proteins as columns, - removes
proteins with zero variance, - runs `prcomp()` (center = TRUE, scale. =
TRUE). - Produces PCA plots (via `ggfortify::autoplot`) for Day 3, Day 7
and combined; groups and shapes indicate condition and day. - Extracts
top contributing proteins to PC1 and PC2 (by absolute loading), maps
Accessions → Gene Symbol, computes relative abundances (WT_rel, SE_rel,
SA_rel), and writes this table to: - `Top_proteins_PCA_abun_rel.xlsx`

------------------------------------------------------------------------

## Inputs (expected files in `Data/`)

-   Proteome Discoverer export: `Results_amostras.xlsx` (must contain
    abundance columns `Abundances (Grouped): F1` … `F18`, `Accession`,
    `Master`, `Contaminant`, `Protein FDR Confidence: Combined`,
    `# Unique Peptides`, `Coverage [%]`, `MW [kDa]`, `Gene Symbol`,
    etc.)
-   BioGRID tab file (used to check APP interactors)
-   APID Excel file `Interactome_APP_3_3_2_3_APID_single_prot.xlsx`
-   G:Profiler Excel sheets used by `pathway_plots_script.R` (as present
    in your `Data/` folder)
-   Any curated lists (e.g., `20_pathways_escolhidos.xlsx`, common
    pathway files) used by the plotting script

> **Tip:** place all input Excel / tab files under `Data/` and run
> scripts from the root `Thesis_files/` folder so relative paths match
> (or edit the top-of-script paths to point to
> `file.path(getwd(), "Data", "...")`).

------------------------------------------------------------------------

## How to run (examples)

1.  Install R (recommended ≥ 4.0) and required packages:

``` r
install.packages(c("readxl", "dplyr", "VennDiagram", "ggplot2", "writexl",
                   "scales", "ggfortify", "factoextra", "tibble"))
```

2.  Ensure working directory is the repo root (where `Data/` and
    `Scripts/` live). In R:

``` r
setwd("path/to/Thesis_files")
getwd()  # verify
```

3.  Run main analysis (from terminal or R):

``` bash
Rscript Scripts/thesis_analysis_script.R
```

4.  Run pathway plotting script:

``` bash
Rscript Scripts/pathway_plots_script.R
```

Alternatively open each script in RStudio and run interactively to view
plots and tweak parameters.

------------------------------------------------------------------------

## Key thresholds, rules and conventions used

-   Filter: `Protein FDR Confidence: Combined == "High"`.
-   Remove contaminant entries (`Contaminant == TRUE`) and a small
    hard-coded list of known contaminant Accessions.
-   Exclude proteins if (`MW [kDa] > 20` AND `# Unique Peptides == 1`)
    OR (`Coverage [%] < 20` AND `# Unique Peptides == 1`).
-   Fold-change: computed vs `N1`; proteins are kept for downstream
    per-condition sets when `FC > 1.5`.
-   Zero-handling for FC: denominator = 0 and numerator \> 0 → FC =
    `100` (sentinel for large/infinite).
-   Intersection enrichment: proteins in SA∩SE with `ratio_SA_SE >= 2`
    are moved to SA-only; those with `ratio_SE_SA >= 2` moved to
    SE-only. Remaining intersection proteins are considered "non-DE
    intersection" for Venn sizing.

------------------------------------------------------------------------

## Outputs produced (examples)

-   `Dados Proteinas antes da filtragem.xlsx` (intermediate)
-   `Lista de Proteinas(I,P,D)_FC_1,5.xlsx` (group Accessions per day)
-   `Top_proteins_PCA_abun_rel.xlsx` (top PCA contributors with relative
    abundances)
-   Multiple plot objects rendered to the R graphics device (barplots,
    Venn diagrams, pathway plots, PCA plots). Enable `ggsave()` in
    scripts to export PNGs.

------------------------------------------------------------------------

## Notes, reproducibility & tips

-   Scripts use relative paths by default but include a
    `base_dir <- getwd()` / `file.path()` pattern — verify your R
    working directory before running, or edit the path variables at the
    top of each script.
-   If plots overlap (e.g., `VennDiagram` objects), use
    `grid::grid.newpage()` before drawing each Venn or generate Venns to
    file (`venn.diagram(..., filename = "venn.png")`).
-   Scripts are modular — it is straightforward to run parts
    individually (preprocessing, FC calculation, pathway plotting, PCA).
-   If you plan to publish or share data publicly, verify any data usage
    licenses for BioGRID/APID resources and remove sensitive or
    restricted data from `Data/` before publishing.

------------------------------------------------------------------------

## Contact / author

For questions about how the scripts were implemented or for help
reproducing results, contact the thesis author (details in the thesis
document).
