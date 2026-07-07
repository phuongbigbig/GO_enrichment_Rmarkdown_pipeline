# GO & Pathway Enrichment Pipeline for Non-Model Organisms

A complete R Markdown pipeline for **Gene Ontology (GO)** and **pathway**
over-representation analysis (ORA), designed for any organism in NCBI Gene — no
OrgDb package required for the GO half. GO enrichment runs with
`clusterProfiler::enricher()` and DAG-aware `topGO`; pathway enrichment covers
**KEGG, Reactome, WikiPathways and PANTHER**, with automatic ortholog mapping to
a model-organism reference for the databases that need one.

Pipeline file: **`GO_Enrichment.concise.Rmd`**

---

## Table of Contents

1. [Overview](#overview)
2. [What the pipeline produces](#what-the-pipeline-produces)
3. [Directory Structure](#directory-structure)
4. [Prerequisites](#prerequisites)
5. [Step 1 — Find Your Organism's TaxID](#step-1--find-your-organisms-taxid)
6. [Step 2 — Build the Annotation Files](#step-2--build-the-annotation-files)
7. [Step 3 — Prepare Your DE Gene List](#step-3--prepare-your-de-gene-list)
8. [Step 4 — Configure Parameters](#step-4--configure-parameters)
9. [Step 5 — Install R Packages](#step-5--install-r-packages)
10. [Step 6 — Run the Pipeline](#step-6--run-the-pipeline)
11. [Output Files](#output-files)
12. [Gene-set eligibility filters](#gene-set-eligibility-filters)
13. [Pathway enrichment & ortholog mapping](#pathway-enrichment--ortholog-mapping)
14. [Adapting to Another Species](#adapting-to-another-species)
15. [Troubleshooting](#troubleshooting)
16. [Citation](#citation)

---

## Overview

The pipeline takes a list of differentially expressed (DE) genes (as NCBI Gene
IDs or gene symbols) and tests which GO terms and pathways are over-represented
relative to a background set. It is organised into six clear stages:

1. **Setup** — every user-configurable parameter in a single chunk at the top,
   plus packages and helper functions.
2. **Input data** — load pre-built annotation files and the DE / background
   gene lists.
3. **GO enrichment** — `enricher()` for BP / MF / CC, and DAG-aware `topGO`.
4. **Pathway enrichment** — KEGG (direct), and Reactome / WikiPathways / PANTHER
   via orthologs to a model-organism reference.
5. **Validation & benchmarking** — a manual hypergeometric test that reproduces
   `enricher()` for all three ontologies, plus an `enricher()`-vs-`topGO`
   overlap and significance concordance.
6. **Summary** — the parameters used and a single at-a-glance count of every
   enrichment result across all methods.

Two complementary GO methods are run:

- **`clusterProfiler::enricher()`** — one-sided hypergeometric (Fisher) test with
  multiple-testing correction, run separately for Biological Process (BP),
  Molecular Function (MF) and Cellular Component (CC).
- **`topGO`** — DAG-aware enrichment using the `weight01` algorithm, which
  accounts for the hierarchical structure of the Gene Ontology and reduces
  redundancy between parent and child terms.

---

## What the pipeline produces

- A **self-contained HTML report** (tabbed tables and a curated set of
  visualisations: dot plots, bar plots, enrichment maps, gene-concept networks,
  heatmaps, an enrichment volcano, topGO comparison plots, and combined
  cross-database figures).
- All results as **TSV files** and all figures as **300-dpi PNGs** under
  `GO_results/`.

Only two things are written to disk: the HTML report and the `GO_results/`
directory. Caching is **off**, so every run is reproducible from scratch.

---

## Directory Structure

```
project/
├── GO_Enrichment.concise.Rmd            # The R Markdown pipeline
├── Phypophthalmus_gene2go.tsv           # Organism GO annotations (required)
├── Phypophthalmus_gene_info.tsv         # Symbol-to-ID lookup (optional)
├── Phypophthalmus_to_drerio_orthologs.tsv  # Ortholog map (optional; auto-built & cached)
├── de_gene_list.txt                     # Your DE gene list (you provide this)
├── background_genes.txt                 # Optional custom background
└── GO_results/                          # Output directory (created by pipeline)
    ├── GO_enrichment_BP.tsv
    ├── GO_enrichment_MF.tsv
    ├── GO_enrichment_CC.tsv
    ├── GO_enrichment_combined.tsv
    ├── topGO_BP.tsv / topGO_MF.tsv / topGO_CC.tsv
    ├── pathway_results/
    │   ├── KEGG_enrichment.tsv
    │   ├── Reactome_enrichment.tsv
    │   ├── WikiPathways_enrichment.tsv
    │   ├── PANTHER_enrichment.tsv
    │   ├── pathway_enrichment_combined.tsv
    │   └── ortholog_mapping_used.tsv
    └── plots/
        ├── dotplot_BP.png, barplot_BP.png, emap_BP.png, cnetplot_BP.png, ...
        ├── integrated_summary.png
        └── pathways/
            ├── dotplot_KEGG.png, barplot_KEGG.png, dotplot_Reactome.png, ...
```

---

## Prerequisites

**System tools** (for the one-time file preparation):

- `bash`, `awk`, `gunzip`/`zcat` (standard on macOS and Linux)
- `wget` or `curl` (to download NCBI files)

**R packages** (installed in Step 5):

- CRAN: `tidyverse`, `DT`, `scales`, `ggrepel`, `rbioapi`
- Bioconductor (GO): `clusterProfiler`, `enrichplot`, `GO.db`, `AnnotationDbi`, `topGO`
- Bioconductor (pathways): `ReactomePA`, `DOSE`, `biomaRt`, `KEGGREST`, and a
  **reference-organism** OrgDb such as `org.Dr.eg.db` (zebrafish),
  `org.Gg.eg.db` (chicken) or `org.Hs.eg.db` (human) — match this to your
  `ORTHO_ORGDB` parameter.

> The GO half runs with only the "GO" packages. The pathway half is optional and
> can be turned off with `RUN_PATHWAY_ENRICHMENT <- FALSE`.

---

## Step 1 — Find Your Organism's TaxID

Every organism in NCBI has a unique Taxonomy ID (TaxID), needed to extract the
right rows from the NCBI reference files.

- **Search NCBI Taxonomy:** <https://www.ncbi.nlm.nih.gov/taxonomy/> — the TaxID
  is the number shown next to the species name.

| TaxID  | Species                          | Common name              |
|--------|----------------------------------|--------------------------|
| 310915 | *Pangasianodon hypophthalmus*    | Striped catfish          |
| 7998   | *Ictalurus punctatus*            | Channel catfish          |
| 8128   | *Oreochromis niloticus*          | Nile tilapia             |
| 7955   | *Danio rerio*                    | Zebrafish (common ref.)  |
| 9031   | *Gallus gallus*                  | Chicken (common ref.)    |
| 9606   | *Homo sapiens*                   | Human (common ref.)      |
| 10090  | *Mus musculus*                   | Mouse                    |

---

## Step 2 — Build the Annotation Files

The script **loads** pre-built, organism-filtered files; it does not download
anything itself. Prepare them once with the commands below (these mirror the
"Obtaining the annotation files" section inside the Rmd). Replace `310915` with
your TaxID.

**1. GO annotations — `*_gene2go.tsv` (required):**

```bash
wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz
gunzip gene2go.gz
# Keep the header (NR==1) plus rows for your TaxID (column 1)
awk -F'\t' 'NR==1 || $1==310915' gene2go > Phypophthalmus_gene2go.tsv
# Columns: TaxID GeneID GO_ID Evidence Qualifier GO_term PubMed Category
```

**2. Symbol ↔ ID lookup — `*_gene_info.tsv` (only if your DE list uses symbols):**

```bash
wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz
gunzip gene_info.gz
awk -F'\t' 'NR==1 || $1==310915' gene_info > gene_info_sp.txt
# Keep GeneID, Symbol, Synonyms, description, type_of_gene → Phypophthalmus_gene_info.tsv
```

**3. Ortholog map — `*_to_<ref>_orthologs.tsv` (only for pathway enrichment):**
Needed for Reactome / WikiPathways / PANTHER (and KEGG-via-orthologs). If absent,
the pipeline builds it automatically (Ensembl BioMart, falling back to symbol
matching) and caches it. To build it offline, the file needs just two columns,
`sp_geneid` (your species' NCBI GeneID) and `ref_entrez` (reference GeneID):

```bash
wget https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_orthologs.gz
gunzip gene_orthologs.gz
# Your species (310915) ↔ reference (zebrafish 7955), both directions:
{ echo -e "sp_geneid\tref_entrez"
  awk -F'\t' '$1==310915 && $4==7955 {print $2"\t"$5}' gene_orthologs
  awk -F'\t' '$1==7955 && $4==310915 {print $5"\t"$2}' gene_orthologs
} | sort -u > Phypophthalmus_to_drerio_orthologs.tsv
```

> After building the small filtered files, you can delete the large genome-wide
> downloads — only the per-organism TSVs are used.

---

## Step 3 — Prepare Your DE Gene List

A plain-text file with one gene per line. The pipeline auto-detects whether the
input is numeric NCBI Gene IDs or gene symbols.

```
# de_gene_list.txt — NCBI Gene IDs (no lookup table needed)
113547368
113547513
113547290
```

```
# de_gene_list_symbols.txt — gene symbols (requires the gene_info file)
tlr3
irf9
nfkb1
```

**Rules:** one gene per line; `#` starts a comment; inline comments after a `#`
are stripped; blank lines ignored; do not mix IDs and symbols.

**From DESeq2:**

```r
sig_genes <- res %>% as.data.frame() %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>% rownames()
writeLines(sig_genes, "de_gene_list.txt")
```

**Optional custom background** — the set of genes actually expressed/testable in
your experiment (recommended over the default of all annotated genes):

```r
all_tested <- rownames(res)[!is.na(res$padj)]
writeLines(all_tested, "background_genes.txt")
```

---

## Step 4 — Configure Parameters

**All user settings live in a single `params` chunk at the top of the Rmd** — no
need to hunt through the code. The main ones:

```r
## Organism
TAXID          <- 310915
SPECIES_NAME   <- "Pangasianodon hypophthalmus"

## Input files
GENE2GO_FILE    <- "Phypophthalmus_gene2go.tsv"
GENE_INFO_FILE  <- "Phypophthalmus_gene_info.tsv"
DE_GENE_FILE    <- "de_gene_list_symbols.txt"
BACKGROUND_FILE <- NULL                        # or "background_genes.txt"
ORTHO_FILE      <- "Phypophthalmus_to_drerio_orthologs.tsv"

## Statistical thresholds
PVALUE_CUTOFF <- 0.05
QVALUE_CUTOFF <- 0.2
PADJ_METHOD   <- "BH"

## Gene-set eligibility filters (see below)
MIN_DE_IN_TERM <- 2      # ≥ this many DE genes per term
MIN_GENE_SET   <- 5      # ≥ this many background genes per term
MAX_GENE_SET   <- 500    # ≤ this many background genes per term

## Which analyses to run
RUN_TOPGO              <- TRUE
RUN_PATHWAY_ENRICHMENT <- TRUE

## Pathway / ortholog settings (see "Adapting to Another Species")
KEGG_ORG           <- "phyp"
KEGG_VIA_ORTHOLOGS <- FALSE
REACTOME_ORG       <- "zebrafish"
WP_ORGANISM        <- "Danio rerio"
PANTHER_TAXID      <- 7955
ORTHO_ORGDB           <- "org.Dr.eg.db"
ORTHO_BIOMART_PREFIX  <- "drerio"
ORTHO_ENSEMBL_DATASET <- "drerio_gene_ensembl"
SPECIES_BIOMART_REGEX <- "hypophthalmus|pangasianodon"
```

---

## Step 5 — Install R Packages

Run once in R / RStudio (swap the reference OrgDb to match `ORTHO_ORGDB`):

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "enrichplot", "GO.db", "AnnotationDbi",
                       "topGO", "ReactomePA", "DOSE",
                       "org.Dr.eg.db",           # reference OrgDb — e.g. org.Gg.eg.db, org.Hs.eg.db
                       "biomaRt", "KEGGREST"))

install.packages(c("tidyverse", "DT", "scales", "ggrepel", "rbioapi"))
```

---

## Step 6 — Run the Pipeline

**RStudio:** open the `.Rmd` and click **Knit** (Ctrl/Cmd+Shift+K).

**R console:**

```r
rmarkdown::render("GO_Enrichment.concise.Rmd")
```

**Terminal:**

```bash
Rscript -e 'rmarkdown::render("GO_Enrichment.concise.Rmd")'
```

Produces a self-contained HTML report and writes all results/plots to
`GO_results/`.

---

## Output Files

| File | Description |
|------|-------------|
| `GO_Enrichment.concise.html` | Full interactive HTML report |
| `GO_results/GO_enrichment_{BP,MF,CC}.tsv` | `enricher()` results per ontology |
| `GO_results/GO_enrichment_combined.tsv` | All ontologies + fold enrichment |
| `GO_results/topGO_{BP,MF,CC}.tsv` | topGO results (classic + weight01 + elim) |
| `GO_results/pathway_results/KEGG_enrichment.tsv` | KEGG pathways |
| `GO_results/pathway_results/Reactome_enrichment.tsv` | Reactome pathways |
| `GO_results/pathway_results/WikiPathways_enrichment.tsv` | WikiPathways |
| `GO_results/pathway_results/PANTHER_enrichment.tsv` | PANTHER pathways |
| `GO_results/pathway_results/pathway_enrichment_combined.tsv` | All databases combined |
| `GO_results/pathway_results/ortholog_mapping_used.tsv` | The ortholog map applied |
| `GO_results/plots/*.png` | GO plots (300-dpi) + `integrated_summary.png` |
| `GO_results/plots/pathways/*.png` | Pathway plots (300-dpi) |

---

## Gene-set eligibility filters

Two thresholds control which terms are testable/reported (both user-editable at
the top):

- **`MIN_DE_IN_TERM`** (default 2) — a term must be hit by at least this many of
  your DE genes to be reported. `enricher()` / `enrichKEGG()` / etc. have no such
  argument (they can report a term hit by a single gene), so the pipeline
  enforces it as a post-hoc filter on the `Count` column across **every** method
  (GO, KEGG, Reactome, WikiPathways, PANTHER, topGO, and the manual test).
- **`MIN_GENE_SET`** / **`MAX_GENE_SET`** (default 5 / 500) — a term's
  *background* size must fall in this range to be tested. Passed to `enricher()`
  et al. as `minGSSize` / `maxGSSize` and to `topGO` as `nodeSize`.

---

## Pathway enrichment & ortholog mapping

KEGG can run **directly** when your species is a registered KEGG organism
(`KEGG_VIA_ORTHOLOGS = FALSE`). Reactome, WikiPathways and PANTHER cover only
model organisms, so the pipeline maps your genes to a **reference organism** via
orthologs. The ortholog map is built once (Ensembl BioMart → symbol-matching
fallback) and cached to `ORTHO_FILE`; delete that file to rebuild it.

Ordering note: ortholog mapping runs **before** KEGG, so KEGG can also be run on
orthologs (`KEGG_VIA_ORTHOLOGS = TRUE`, with `KEGG_ORG` set to the reference
code, e.g. `gga`) when your species is absent from KEGG.

> **Different references, different backgrounds.** Each database's p-values are
> relative to its own universe. PANTHER in particular uses its own whole-genome
> reference, so cross-database comparisons in the Summary are qualitative, not
> like-for-like.

---

## Adapting to Another Species

The pipeline is not catfish-specific — change **parameters, not code**. The key
decision is the **reference organism**, which must be supported by each pathway
database you want to use.

**Example: a non-model bird → chicken (*Gallus gallus*).**

**Step 1.** Rebuild `*_gene2go.tsv` (and `*_gene_info.tsv`) for your bird's
TaxID as in [Step 2](#step-2--build-the-annotation-files).

**Step 2.** Check reference-organism support per database:

| Database | Check supported organisms | Chicken? |
|---|---|---|
| KEGG | <https://www.genome.jp/kegg/catalog/org_list.html> (code `gga`) | ✅ `gga` |
| Reactome (`ReactomePA`) | fixed set: human, mouse, rat, zebrafish, fly, celegans, yeast | ❌ use **human** |
| WikiPathways | `clusterProfiler::get_wp_organisms()` | ✅ `Gallus gallus` |
| PANTHER | `rbioapi::rba_panther_info()` | ✅ TaxID `9031` |

Because chicken is **not** in `ReactomePA`, map Reactome to **human** (build a
second ortholog file, or set `REACTOME_ORG` and skip if unavailable), while
KEGG / WikiPathways / PANTHER use chicken.

**Step 3.** Set the parameters:

```r
TAXID                 <- <your_bird_taxid>
SPECIES_NAME          <- "<Your bird species>"
GENE2GO_FILE          <- "<yourbird>_gene2go.tsv"
GENE_INFO_FILE        <- "<yourbird>_gene_info.tsv"
ORTHO_FILE            <- "<yourbird>_to_ggallus_orthologs.tsv"

KEGG_ORG              <- "gga"        # chicken KEGG code
KEGG_VIA_ORTHOLOGS    <- TRUE         # if your bird is not in KEGG
REACTOME_ORG          <- "human"      # chicken unsupported by ReactomePA
WP_ORGANISM           <- "Gallus gallus"
PANTHER_TAXID         <- 9031

ORTHO_ORGDB           <- "org.Gg.eg.db"          # BiocManager::install("org.Gg.eg.db")
ORTHO_BIOMART_PREFIX  <- "ggallus"
ORTHO_ENSEMBL_DATASET <- "ggallus_gene_ensembl"
SPECIES_BIOMART_REGEX <- "<yourbird_ensembl_pattern>"   # or "" to skip BioMart
```

**Step 4.** Build the ortholog file — let the pipeline build it (needs `biomaRt`
+ your `ORTHO_ORGDB`), or build it offline with the `gene_orthologs` recipe using
your TaxID and `9031`. Find your species' Ensembl dataset name with:

```r
library(biomaRt); ens <- useMart("ensembl")
subset(listDatasets(ens), grepl("finch|taeniopygia", dataset, ignore.case = TRUE))
```

If your species is **not in Ensembl**, set `SPECIES_BIOMART_REGEX <- ""` to use
symbol matching only. For truly non-model species, generate orthologs with
**eggNOG-mapper** or **OrthoFinder** and save the two-column
`sp_geneid`/`ref_entrez` TSV.

**Step 5.** Knit — the ortholog-mapping and pathway chunks read all of the above
from the parameters; no chunk edits required.

---

## Troubleshooting

### "File not found" errors
Ensure the annotation files are in the same directory as the `.Rmd` (or use
absolute paths). Check with `list.files(pattern = "\\.tsv$")`.

### `fisher.test`/hypergeometric error with a custom background
Fixed: the manual validation restricts each term's gene set to the universe
(matching `enricher()`) and uses `phyper()`, so a smaller custom background no
longer produces an invalid contingency table.

### `AnnotationDbi::select` masks `dplyr::select`
Handled automatically (`select <- dplyr::select` in the setup chunk).

### `enrichplot::barplot` not found
The pipeline uses a custom `enrich_barplot()`; no action needed.

### KEGG / Reactome / PANTHER return nothing or error
These call remote APIs — transient failures are caught and the section is
skipped rather than halting the knit. For Reactome, confirm your `REACTOME_ORG`
is in the supported set. If your species has few orthologs, the
`ENOUGH_ORTHOLOGS` guard skips the ortholog-based databases.

### topGO / no significant terms
Small or sparsely annotated gene lists may yield nothing, especially under the
conservative `weight01`. Consider relaxing `PVALUE_CUTOFF`, lowering
`MIN_DE_IN_TERM`, supplementing annotations (eggNOG-mapper / InterProScan), or
using the default all-annotated background if your custom one is very small.

### Caching
Caching is **off** (`cache = FALSE`) so results never go stale between runs. A
full run may take a few minutes; the API-backed pathway steps dominate the time.

---

## Citation

If you use this pipeline, please cite the underlying tools:

- **clusterProfiler:** Wu T, et al. (2021). clusterProfiler 4.0: A universal
  enrichment tool for interpreting omics data. *The Innovation*, 2(3), 100141.
- **topGO:** Alexa A, Rahnenführer J (2009). Gene set enrichment analysis with
  topGO. *Bioconductor Vignettes*.
- **ReactomePA:** Yu G, He QY (2016). ReactomePA: an R/Bioconductor package for
  reactome pathway analysis and visualization. *Mol. BioSyst.*, 12(2), 477–479.
- **PANTHER / rbioapi:** Rezwani M, et al. (2022). rbioapi: user-friendly R
  interface to biologic web services' API. *Bioinformatics*, 38(10), 2952–2953.
- **Gene Ontology:** Gene Ontology Consortium (2023). The Gene Ontology
  knowledgebase in 2023. *Genetics*, 224(1), iyad031.
```
