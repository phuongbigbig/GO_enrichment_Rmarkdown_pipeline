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
16. [Changelog — robustness and edge-case fixes](#changelog--robustness-and-edge-case-fixes)
17. [Citation](#citation)

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
- Bioconductor (pathways): `ReactomePA`, `DOSE`, `biomaRt`, and a
  **reference-organism** OrgDb such as `org.Dr.eg.db` (zebrafish),
  `org.Gg.eg.db` (chicken) or `org.Hs.eg.db` (human) — match this to your
  `ORTHO_ORGDB` parameter. `KEGGREST` is optional: `enrichKEGG()` downloads
  KEGG data itself.
- `biomaRt` and the OrgDb are only needed to **build** the ortholog map. If you
  already have `ORTHO_FILE` (e.g. built offline from NCBI `gene_orthologs`),
  the ortholog-based databases run without them.

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

A map built by the pipeline also gets a `ref_orgdb` column (e.g.
`org.Dr.eg.db`) that records which reference it was built for. If you later
change `ORTHO_ORGDB`, the pipeline stops with a clear message instead of reusing
the wrong map. Hand-built files without this column are accepted as they are,
so make sure they match `ORTHO_ORGDB` yourself.

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
are stripped; blank lines ignored; duplicate lines are removed. A list that
mixes numeric IDs and symbols is accepted: numeric IDs are kept as they are and
symbols are converted (this needs the gene_info file). The report prints how
many of each were found.

DE genes that are not in the background are dropped (the report prints how
many). The knit stops with a clear message if **no** DE gene is left, or if
**every** background gene is a DE gene, because over-representation cannot be
tested in either case.

**From DESeq2:**

```r
sig_genes <- res %>% as.data.frame() %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1) %>% rownames()
writeLines(sig_genes, "de_gene_list.txt")
```

**Optional custom background** — the set of genes actually expressed/testable in
your experiment (recommended over the default of all annotated genes). If
`BACKGROUND_FILE` is set but the file does not exist (e.g. a typo), the pipeline
prints a warning, falls back to all annotated genes, and the Summary shows the
background that was **actually** used:

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

**One reference per run.** The pipeline builds a single ortholog map, to
`ORTHO_ORGDB`. Every ortholog-based database therefore receives that organism's
Entrez IDs, so `REACTOME_ORG`, `WP_ORGANISM`, `PANTHER_TAXID` (and `KEGG_ORG`
when `KEGG_VIA_ORTHOLOGS = TRUE`) must name the **same** organism. For the
known reference packages (`org.Dr.eg.db`, `org.Gg.eg.db`, `org.Hs.eg.db`,
`org.Mm.eg.db`, `org.Rn.eg.db`) the parameters chunk checks this. A database
that points at a different organism is **skipped** with a printed warning,
instead of silently returning nothing or nonsense. For other OrgDb packages
no check is made.

**PANTHER is over-representation only.** PANTHER's Fisher test reports both
over- and under-represented pathways. The pipeline keeps only over-represented
ones (`plus_minus == "+"`, or `fold_enrichment > 1` if that column is
missing), so the results match the other ORA databases. Only PANTHER pathways
with FDR < `PVALUE_CUTOFF` enter the combined pathway plots and table.

**Paralogs.** One-to-many orthologs are kept, so a catfish gene that maps to
both zebrafish paralogs (e.g. *gene_a* and *gene_b*) can count twice in a
pathway. This is the usual practice, but check pathway hits that depend on
paralog pairs.

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

Because chicken is **not** in `ReactomePA`, a chicken run cannot include
Reactome. With a chicken ortholog map, `REACTOME_ORG <- "human"` would receive
chicken IDs, so the pipeline detects the mismatch and **skips Reactome** with a
warning. To get Reactome for a bird, knit the document a second time with
`ORTHO_ORGDB <- "org.Hs.eg.db"`, `REACTOME_ORG <- "human"`,
`ORTHO_BIOMART_PREFIX <- "hsapiens"`, `ORTHO_ENSEMBL_DATASET <-
"hsapiens_gene_ensembl"` and a **separate** `ORTHO_FILE` (e.g.
`<yourbird>_to_hsapiens_orthologs.tsv`).

**Step 3.** Set the parameters:

```r
TAXID                 <- <your_bird_taxid>
SPECIES_NAME          <- "<Your bird species>"
GENE2GO_FILE          <- "<yourbird>_gene2go.tsv"
GENE_INFO_FILE        <- "<yourbird>_gene_info.tsv"
ORTHO_FILE            <- "<yourbird>_to_ggallus_orthologs.tsv"

KEGG_ORG              <- "gga"        # chicken KEGG code
KEGG_VIA_ORTHOLOGS    <- TRUE         # if your bird is not in KEGG
REACTOME_ORG          <- "human"      # chicken unsupported by ReactomePA → auto-skipped in this run
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
skipped rather than halting the knit. **The reason is printed in the report**
as a line starting with `Note:` (or `WARNING:` in the parameters chunk) —
look there first. Common reasons:

- the database's organism does not match `ORTHO_ORGDB` (see
  [One reference per run](#pathway-enrichment--ortholog-mapping));
- fewer than `MIN_GENE_SET` DE genes mapped to orthologs (the
  `ortholog-mapping` chunk prints the counts);
- `rbioapi` / `ReactomePA` not installed;
- the API was down or changed its output format (re-run later).

### `object 'fdr' not found` in the `panther-plot` chunk
Fixed in the current version. It happened when PANTHER was skipped or returned
nothing: the empty placeholder table had no columns, and the plot chunk
filtered on `fdr` before checking whether the table had rows.

### `Ortholog file ... was built for org.X but ORTHO_ORGDB is org.Y`
You changed the reference organism but kept the old cached ortholog map.
Delete `ORTHO_FILE` (it will be rebuilt) or point `ORTHO_FILE` at a new name.

### `Every background gene is a DE gene`
Your DE list and background are the same set, so there is nothing to test.
Check that `DE_GENE_FILE` and `BACKGROUND_FILE` point at the right files.

### Pandoc warning `Deprecated: --mathjax`
Harmless and now silenced: the YAML header sets `mathjax: null` because the
report has no LaTeX math in the text (plot axis labels are drawn by R).

### topGO / no significant terms
Small or sparsely annotated gene lists may yield nothing, especially under the
conservative `weight01`. Consider relaxing `PVALUE_CUTOFF`, lowering
`MIN_DE_IN_TERM`, supplementing annotations (eggNOG-mapper / InterProScan), or
using the default all-annotated background if your custom one is very small.

### Caching
Caching is **off** (`cache = FALSE`) so results never go stale between runs. A
full run may take a few minutes; the API-backed pathway steps dominate the time.

---

## Changelog — robustness and edge-case fixes

The pipeline ran correctly on the striped catfish data set before this review.
The fixes below cover cases where it could **stop the knit**, **skip an
analysis without saying why**, or — most importantly — **give a wrong result
without any error**. Every change is marked with a `[FIX]` comment in
`GO_Enrichment.concise.Rmd`; search for `[FIX]` to find them. All R chunks
were checked to parse; the full document should be re-knitted on real data
to confirm.

### A. Fixes that change results (wrong output without any error)

| # | Chunk | Problem | When it happened | Fix |
|---|---|---|---|---|
| A1 | `load-gene2go` | NCBI `gene2go` rows with a **NOT** qualifier (e.g. `NOT enables`, `NOT involved_in`) were counted as positive annotations. These mean the gene is explicitly **not** linked to the term. | Any species whose `gene2go` contains NOT rows. | Rows whose `Qualifier` starts with `NOT` are removed before analysis; the number removed is printed. Rows with missing `GeneID` / `GO_ID` are also dropped. |
| A2 | `panther-run` | PANTHER's Fisher test returns both **over-** and **under-represented** pathways. Strongly *depleted* pathways were reported as significant hits. | Any run where PANTHER found a depleted pathway. | Keep only `plus_minus == "+"` (or `fold_enrichment > 1` if that column is absent). |
| A3 | `topgo-run` | `GenTable(topNodes = min(50, ...))` kept only the **top 50** terms per ontology. topGO counts were capped at 50, and the enricher-vs-topGO overlap was biased toward "enricher only". | Any ontology with > 50 scored terms. | `topNodes` = all scored nodes. |
| A4 | `topgo-run` | `GenTable` cuts term names to **40 characters** by default, so long names were truncated in tables and plots. | Any GO term name > 40 characters. | `numChar = 1000`. |
| A5 | `topgo-run` | topGO results had no upper size limit, unlike `enricher()` (`maxGSSize`), so the two methods tested different sets of terms. | Terms with > `MAX_GENE_SET` annotated genes. | Filter `Annotated <= MAX_GENE_SET`. P-values reported as `"< 1e-30"` are now parsed robustly (`as.character()` before removing `<`). |
| A6 | `topgo-plots` | The weight01-vs-classic scatter said "below diagonal = **more** significant after DAG correction". With weight01 on the y-axis, points below the diagonal are **less** significant (redundant parent terms are down-weighted). | Always (interpretation only). | Comment and subtitle corrected. |
| A7 | `go-combined`, `integrated-summary` | When the **same term name** appeared in two sources (e.g. "Ribosome" in GO-CC and KEGG, or a name shared by two ontologies), `geom_col()` **stacked** them into one bar whose length was the **sum** of both −log10 values. | Any duplicated description among the plotted terms. | Duplicated names get a suffix, e.g. `Ribosome [KEGG]`, `... [BP]`. The dot plot / bar plot helpers and the PANTHER plot also use `make.unique()` so factor levels are never duplicated. |
| A8 | `helpers` and all plots | A p-value that underflows to **0** gives `-log10(0) = Inf`; ggplot silently **drops** non-finite values, so the *strongest* term disappeared from the plot. A log10 colour scale also fails on 0. | Very strong enrichments (p < ~1e-308) or APIs that return 0. | New helper `nlog10()` floors p at `.Machine$double.xmin` before taking −log10; used in every −log10 calculation. `enrich_barplot()` floors `p.adjust` the same way. |
| A9 | `pathway-combined-build` | `panther_clean` holds **every** tested PANTHER pathway (`cutoff = 1`), while KEGG / Reactome / WikiPathways contribute only significant terms. PANTHER flooded the combined dot plot, volcano and integrated overview with non-significant pathways. | Any run with PANTHER results. | Only PANTHER rows with `fdr < PVALUE_CUTOFF` enter `all_pathway_df`. `ID` is coerced to character so `bind_rows()` never fails on type differences. |
| A10 | all top-N selections | `slice_min(..., n = k)` keeps **ties**, so "top 10" could return many more rows when several terms share a p-value, distorting plot sizes. | Tied p-values (common with small gene counts). | `with_ties = FALSE` everywhere. |
| A11 | `input-gene-list` | In a list that **mixed** numeric IDs and symbols, the numeric IDs were sent through symbol matching, never matched, and were **silently lost**. Duplicate lines also inflated the reported counts. | Mixed or duplicated input lists. | Input is de-duplicated; numeric IDs are kept as they are and only the symbols are converted; the split is printed. Background is also de-duplicated. |
| A12 | `params` | If an ortholog-based database pointed at a **different organism** than `ORTHO_ORGDB` (e.g. chicken map but `REACTOME_ORG = "human"`), it received the wrong organism's IDs and returned nothing or nonsense. | Mixed reference set-ups, as the old README suggested for birds. | New reference-consistency check (`REF_LOOKUP` table) for zebrafish, chicken, human, mouse and rat OrgDbs. Mismatched databases are **skipped** with a printed `WARNING:`. Unknown OrgDbs are not checked. |

### B. Fixes that stop the knit from failing or skipping work

| # | Chunk | Problem | When it happened | Fix |
|---|---|---|---|---|
| B1 | `panther-run`, `panther-plot` | `panther_clean <- tibble()` had **no columns**. `panther-plot` ran `filter(fdr < PVALUE_CUTOFF)` before its `nrow() > 0` guard → `object 'fdr' not found`. | PANTHER skipped, failed, or returned nothing. | Typed empty placeholder `panther_empty` (all expected columns) and the guard now runs **before** filtering. |
| B2 | `panther-run` | Cleaning of the PANTHER response was outside `tryCatch`; a renamed column in the API output stopped the whole knit. | PANTHER API format change. | Cleaning is wrapped in `tryCatch`; on failure a note is printed and the empty table is used. |
| B3 | `manual-bp/mf/cc` | `manual_results[["CC"]] %>% filter(...)` stopped the knit when an ontology had **no annotations** for the species (`NULL`). | Sparse annotations (e.g. no CC terms). | New NULL-safe `manual_table(ont)` helper used by all three chunks. |
| B4 | `pathway-packages` | `RUN_ORTHO_PATHWAYS` required **both** `biomaRt` and the OrgDb, even when a ready `ORTHO_FILE` existed. All ortholog-based databases were skipped for no reason. | Offline-built ortholog file without biomaRt installed. | Runs if `ORTHO_FILE` exists, **or** biomaRt is available, **or** OrgDb + gene_info allow symbol matching. The BioMart step checks for `biomaRt` itself. |
| B5 | `pathway-packages`, `kegg-run` | KEGG was gated on `KEGGREST` being installed, although `enrichKEGG()` does not need it. | KEGGREST not installed. | KEGG no longer depends on KEGGREST. |
| B6 | `ortholog-mapping` | A `NULL` ortholog map (all strategies failed) crashed the following `filter()`. Hand-built files with blank / non-numeric IDs gave type mismatches. | No map could be built, or a hand-made TSV. | `NULL` → typed empty table; `sp_geneid` / `ref_entrez` coerced to integer and NA rows dropped. |
| B7 | `topgo-run` | The gene factor was built without fixed levels; topGO needs a 2-level factor. An ontology with **no GO nodes** after pruning caused `GenTable` errors. | Edge cases with tiny gene sets. | `factor(..., levels = c(0, 1))`; ontologies with 0 nodes are skipped with a note. |
| B8 | `load-gene2go`, `load-gene-info`, `ortholog-mapping` | `read.delim()` treats `"` as a quote character. A stray quote in a GO term name or gene description **merged lines silently**. | Any text field containing `"` or `'`. | `quote = ""` (and `comment.char = ""`) for gene2go, gene_info and the ortholog file. |
| B9 | `load-gene2go`, `load-gene-info` | Missing columns in the gene_info file only failed much later with an unclear error; an empty gene2go after filtering was not caught. | Wrong or truncated input files. | Early, explicit `stop()` messages naming the missing columns / empty file. |
| B10 | `input-gene-list` | No check for a degenerate test where every background gene is DE. | DE file = background file by mistake. | Clear `stop()` (checked after the "no DE genes left" check). |
| B11 | `params` | Logical switches were not validated. | Typo such as `RUN_TOPGO <- "yes"`. | `stopifnot(is.logical(...))` for `KEGG_VIA_ORTHOLOGS`, `RUN_TOPGO`, `RUN_PATHWAY_ENRICHMENT`. |

### C. Fixes to reporting (the report now says what really happened)

| # | Chunk | Problem | Fix |
|---|---|---|---|
| C1 | `helpers` (used everywhere) | The setup chunk sets `message = FALSE`, so **every** `message("... skipped / failed ...")` was invisible in the HTML. It was impossible to see *why* PANTHER or Reactome was empty. | New helper `say()` prints `Note: ...` into the chunk output; all status/failure messages use it (only "Saved: file" messages stay hidden). |
| C2 | `input-gene-list`, `summary-params` | A misspelled `BACKGROUND_FILE` fell back silently to all annotated genes, while the Summary still showed the file name as the background used. | Prints a `WARNING:`; new `BACKGROUND_USED` variable is what the Summary reports. |
| C3 | `ortholog-mapping` | Changing `ORTHO_ORGDB` silently reused the old cached map (e.g. zebrafish IDs sent to chicken databases). | Newly built maps get a `ref_orgdb` column; a mismatch stops with a clear message. Maps built before this change have no tag and cannot be checked. |
| C4 | `reactome-run`, `wp-run`, `panther-run`, `pathway-summary-table`, `pathway-combined-plots` | Captions, subtitles and the Strategy column were hard-coded as "zebrafish orthologs" / "Direct (phyp)". | All derived from the parameters (`REF_LABEL`, `KEGG_ORG`, `KEGG_VIA_ORTHOLOGS`). |
| C5 | `go-export` | `imap()` returned a list of `NULL`s that was auto-printed into the report as `$BP NULL $CC NULL ...`. | `iwalk()`. |
| C6 | YAML header | Pandoc 3 printed `[WARNING] Deprecated: --mathjax`. | `mathjax: null` (the report has no LaTeX math in the text). |
| C7 | Rmd text, "Adapting to another species" | The text suggested mixing a chicken map with human Reactome in one run, which A12 now (correctly) skips. | Rewritten: one reference per run; run the document again for human Reactome. |

### Known limitations (not changed)

- **Paralogs:** one-to-many orthologs are kept, so a gene can count twice in a
  pathway through paralog pairs (see
  [Pathway enrichment & ortholog mapping](#pathway-enrichment--ortholog-mapping)).
- **Reference check coverage:** only the five OrgDb packages listed in A12 are
  checked; other references are trusted as configured.
- **Untagged ortholog caches:** files built before this version (or by hand)
  have no `ref_orgdb` column, so C3 cannot detect a mismatch for them.
- **PANTHER background:** PANTHER still uses its own whole-genome reference
  rather than your custom background, so its p-values are not directly
  comparable with the other databases.

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
