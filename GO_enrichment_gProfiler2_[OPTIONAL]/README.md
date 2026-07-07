# Gene Ontology & Pathway Enrichment Analysis Pipeline

## Overview

This R Markdown pipeline performs **functional over-representation analysis (ORA)** on a user-supplied gene list using the [g:Profiler](https://biit.cs.ut.ee/gprofiler/gost) web service via its R client package `gprofiler2`. It tests whether specific biological functions, pathways, or regulatory elements are statistically over-represented among your genes compared to a genomic background.

The pipeline queries multiple annotation databases simultaneously:

| Source   | Database                  | What it covers                              |
|----------|---------------------------|---------------------------------------------|
| GO:BP    | Gene Ontology             | Biological processes (e.g. apoptosis)       |
| GO:CC    | Gene Ontology             | Cellular components (e.g. nucleus)          |
| GO:MF    | Gene Ontology             | Molecular functions (e.g. kinase activity)  |
| KEGG     | KEGG Pathways             | Metabolic and signalling pathways           |
| REAC     | Reactome                  | Curated reaction pathways                   |
| WP       | WikiPathways              | Community-curated pathways                  |
| TF       | TRANSFAC                  | Transcription factor binding motifs         |
| MIRNA    | miRTarBase                | microRNA target genes                       |
| CORUM    | CORUM                     | Protein complexes                           |
| HPA      | Human Protein Atlas       | Tissue/cell-type expression                 |
| HP       | Human Phenotype Ontology  | Disease-associated phenotypes               |

The availability of each source depends on the selected organism. Model organisms (human, mouse, zebrafish) have the most comprehensive annotations.


## Requirements

### Software

- **R** >= 4.1.0
- **RStudio** (recommended, for interactive knitting)
- **Internet connection** (g:Profiler queries its server remotely)

### R Packages

The following packages are required. The script includes an installation chunk (set `eval=TRUE` to run it once):

| Package        | Purpose                                            |
|----------------|----------------------------------------------------|
| `gprofiler2`   | R client for the g:Profiler enrichment API         |
| `ggplot2`      | Publication-quality plots                          |
| `dplyr`        | Data manipulation                                  |
| `tidyr`        | Data reshaping                                     |
| `DT`           | Interactive HTML tables                            |
| `scales`       | Axis formatting helpers                            |
| `forcats`      | Factor reordering for plot aesthetics              |
| `stringr`      | String wrapping and truncation                     |
| `plotly`       | Interactive plotly-based Manhattan plot             |
| `RColorBrewer` | Colour palettes                                    |
| `ggtext`       | Rich-text facet strip labels                       |
| `tidytext`     | Per-facet axis reordering (`reorder_within`)       |

Install all at once:

```r
install.packages(c(
  "gprofiler2", "ggplot2", "dplyr", "tidyr",
  "DT", "scales", "forcats", "stringr",
  "plotly", "RColorBrewer", "ggtext", "tidytext"
))
```


## Input Preparation

### 1. Gene List File (Required)

Prepare a text file containing your gene identifiers — one per line or as the first column of a tabular file. The file format is auto-detected by extension.

#### Option A: Plain text file (`.txt`)

One gene identifier per line. Blank lines and lines starting with `#` are ignored.

```
# Example: my_genes.txt
ENSG00000141510
ENSG00000012048
ENSG00000139618
TP53
BRCA1
```

#### Option B: Tab-separated file (`.tsv`)

The pipeline reads the **first column** as gene IDs. Additional columns (log2FC, p-values, etc.) are ignored. A header row is expected by default.

```
GeneID    log2FC    padj
TP53      2.31      0.0004
BRCA1     -1.87     0.012
EGFR      3.45      0.0001
```

If your file has **no header row**, set `gene_file_header <- FALSE` in the script.

#### Option C: Comma-separated file (`.csv`)

Same rules as TSV — first column is used, header expected by default.

```
GeneID,log2FC,padj
TP53,2.31,0.0004
BRCA1,-1.87,0.012
```

#### Gene identifier types

g:Profiler auto-detects the identifier namespace. You can use:

- **Ensembl gene IDs**: `ENSG00000141510`, `ENSGALG00000006385`
- **Gene symbols**: `TP53`, `BRCA1`, `MYC`
- **Mixed types**: Ensembl IDs and symbols in the same file

Ensembl IDs are generally preferred as they are unambiguous. Gene symbols can be ambiguous across species or have aliases that may not resolve correctly.

#### Typical sources of gene lists

- Differentially expressed genes (DEGs) from DESeq2 or edgeR
- Genes from a specific cluster in single-cell analysis
- Genes in a QTL region or GWAS hit list
- Genes with specific variant types (e.g. high-impact variants from VEP/SnpEff)


### 2. Background Gene List (Optional)

By default, g:Profiler uses all annotated genes in the selected organism's genome as the statistical background. This is appropriate for most RNA-seq analyses.

However, you should provide a **custom background** when your gene list was drawn from a specific subset of the genome — for example:

- **Only genes detected in your RNA-seq experiment** (i.e. genes with non-zero counts). Using the whole genome as background inflates significance because many undetected genes dilute the expected overlap.
- **Only genes on a microarray chip**
- **Only genes in a specific genomic region** (e.g. a QTL interval)

The background file follows the same format rules as the gene list file (`.txt`, `.tsv`, or `.csv`; first column used).

**Important:** When a custom background is supplied, the script automatically sets `domain_scope = "custom"` in the g:Profiler query. Without this, g:Profiler would silently ignore the background file and use its default annotated genome.


### 3. Organism Selection

Set the `organism` parameter to the g:Profiler organism code for your species. Common codes:

| Code           | Species                                    |
|----------------|--------------------------------------------|
| `"hsapiens"`   | *Homo sapiens* (human)                     |
| `"mmusculus"`  | *Mus musculus* (mouse)                     |
| `"rnorvegicus"`| *Rattus norvegicus* (rat)                  |
| `"drerio"`     | *Danio rerio* (zebrafish)                  |
| `"ggallus"`    | *Gallus gallus* (chicken)                  |
| `"sscrofa"`    | *Sus scrofa* (pig)                         |
| `"btaurus"`    | *Bos taurus* (cattle)                      |
| `"olatipes"`   | *Oryzias latipes* (medaka)                 |
| `"icpunctatus"`| *Ictalurus punctatus* (channel catfish)    |
| `"oniloticus"` | *Oreochromis niloticus* (tilapia)          |
| `"dmelanogaster"` | *Drosophila melanogaster* (fruit fly)   |
| `"celegans"`   | *Caenorhabditis elegans* (nematode)        |
| `"athaliana"`  | *Arabidopsis thaliana* (thale cress)       |

Full list: <https://biit.cs.ut.ee/gprofiler/page/organism-list>


## How to Run

### Step 1: Edit Parameters

Open `GO_Enrichment_gprofiler2.Rmd` in RStudio. Edit these three lines in the **"Gene List & Organism"** section:

```r
gene_file        <- "my_gene_list.txt"    # path to your gene list file
gene_file_header <- FALSE                 # TRUE if TSV/CSV has a header row
organism         <- "hsapiens"            # your organism code
```

If using a custom background:

```r
bg_file <- "my_background_genes.txt"      # uncomment and set the path
```

### Step 2: Place Files

Put your gene list file (and background file, if used) in the **same directory** as the `.Rmd` file, or provide a full/relative path in `gene_file`.

### Step 3: Knit to HTML

**In RStudio:**
- Click the **Knit** button (or press `Ctrl+Shift+K` / `Cmd+Shift+K`)
- Select "Knit to HTML"

**From the R console:**

```r
rmarkdown::render("GO_Enrichment_gprofiler2.Rmd")
```

**From the terminal:**

```bash
Rscript -e 'rmarkdown::render("GO_Enrichment_gprofiler2.Rmd")'
```

The output is an HTML report file in the same directory.

### Step 4: Collect Outputs

After knitting, the working directory will contain:

| File | Description |
|------|-------------|
| `GO_Enrichment_gprofiler2.html` | Interactive HTML report with all tables and figures |
| `gprofiler2_enrichment_<organism>_full_results.tsv` | Full results table (tab-separated) |
| `gprofiler2_enrichment_<organism>_full_results.csv` | Full results table (comma-separated) |
| `gprofiler2_enrichment_<organism>_GO_BP.tsv` | Per-source table: GO Biological Process |
| `gprofiler2_enrichment_<organism>_GO_CC.tsv` | Per-source table: GO Cellular Component |
| `gprofiler2_enrichment_<organism>_GO_MF.tsv` | Per-source table: GO Molecular Function |
| `gprofiler2_enrichment_<organism>_KEGG.tsv` | Per-source table: KEGG Pathways |
| `gprofiler2_enrichment_<organism>_REAC.tsv` | Per-source table: Reactome |
| ... | (one TSV per source with results) |

To additionally save publication-quality PNG/PDF figures, set `eval=TRUE` in the **"Save Publication Figures"** chunk and re-knit.


## Analysis Steps Explained

### Step 1 — Gene List Reading and Validation

The script reads the gene file using a helper function `read_gene_list()` that auto-detects format by file extension. It strips blank lines, comment lines (starting with `#`), trims whitespace, removes duplicates, and reports how many unique identifiers were loaded. The same function handles the optional background file.

### Step 2 — Enrichment Query via g:Profiler

The gene list is sent to the g:Profiler server using `gprofiler2::gost()`. The server performs a **hypergeometric test** (equivalent to a one-tailed Fisher's exact test) for each annotation term in each database. This test asks: given *N* genes in the background, *K* annotated to this term, and *n* genes in my query, is the observed overlap *k* larger than expected by chance?

Key parameters:

- **`correction_method = "g_SCS"`**: g:Profiler's own multiple-testing correction method (Set Counts and Sizes), which accounts for the hierarchical and overlapping structure of GO terms. Alternatives are `"fdr"` (Benjamini-Hochberg) and `"bonferroni"`.
- **`user_threshold = 0.05`**: Adjusted p-value cutoff for significance.
- **`domain_scope`**: Set to `"annotated"` when using the full genome background, or `"custom"` when a user-supplied background list is provided. This switch is handled automatically by the script.
- **`evcodes = TRUE`**: Requests the list of overlapping genes for each significant term, which is useful for downstream interpretation but makes the query slower.

### Step 3 — Result Processing and Metric Computation

The raw g:Profiler results are processed to compute additional metrics not returned by the server:

- **Gene ratio** = `intersection_size / query_size` — the fraction of your query genes that belong to each term.
- **Fold enrichment** = `(intersection_size / query_size) / (term_size / effective_domain_size)` — how many times more often your genes appear in this term than expected by chance. A fold enrichment of 5 means your genes are 5× more likely to be in this term than a random gene set of the same size.
- **-log10(p-value)** — a transformation that makes small (highly significant) p-values appear as large positive numbers, useful for visualisation.
- **log10(fold enrichment)** — log-transformed fold enrichment for bar chart x-axes; fold enrichment values below 1 are floored to 1 (log10 = 0) since the analysis tests for over-representation only.

**Note on `effective_domain_size`:** This value differs between annotation sources. For example, GO:BP might have 18,000 annotated human genes while KEGG has ~8,400. The fold enrichment formula uses the per-row `effective_domain_size` so each term is compared against the correct source-specific background. Using a single global value would inflate or deflate fold enrichments for terms from different databases.

### Step 4 — Table Export

Full results are exported as both TSV (tab-separated, preferred for bioinformatics tools) and CSV (comma-separated, for spreadsheet software). Per-source TSV files are also generated for convenience. The exported tables include all computed metrics (gene ratio, fold enrichment) and, when `evcodes = TRUE`, a comma-separated list of overlapping gene IDs for each term.

### Step 5 — g:Profiler Built-in Visualisations

Three visualisations from the gprofiler2 package:

1. **Interactive Manhattan plot** (`gostplot`, interactive = TRUE) — a plotly-based zoomable, hoverable plot where each dot is an enriched term, arranged by source on the x-axis and -log10(p-value) on the y-axis. Best for exploratory browsing within the HTML report.

2. **Static Manhattan plot** (`gostplot`, interactive = FALSE) — a ggplot2 version suitable for PDF export or journal figures. The top 3 terms per source are labelled using `publish_gostplot()`.

3. **Enrichment table** (`publish_gosttable`) — a formatted table rendered as a plot, showing the top 5 terms per source with colour-coded source indicators.

### Step 6 — Publication-Quality ggplot2 Visualisations

Custom ggplot2 figures designed for journal submission, each addressing different aspects of the enrichment results:

1. **GO bar chart (BP + CC + MF combined)** — horizontal bars coloured by p-value, with log10(fold enrichment) on the x-axis, faceted by GO sub-ontology. Two versions are generated: one using `tidytext::reorder_within()` for correct per-facet ordering, and an alternative that avoids the tidytext dependency.

2. **Pathway-specific bar charts** — separate bar charts for KEGG (green palette), Reactome (purple palette), and WikiPathways (orange palette), each showing the top enriched pathways.

3. **Dot plot (all sources)** — encodes three variables simultaneously: gene ratio on x-axis, dot size proportional to gene count, and colour representing p-value. Horizontal segments connect each dot to the baseline for readability.

4. **Lollipop plot (all sources)** — shows -log10(p-value) on the x-axis with points connected to the y-axis by line segments. Dot size encodes gene count.

5. **Heatmap** — a tile plot showing -log10(p-value) as fill colour for the top 5 terms per source, with gene counts overlaid as text. Terms are grouped by source on the y-axis.

6. **Intersection size bar chart** — the top 25 terms ranked by the number of overlapping genes, coloured by source. Useful for identifying which terms capture the most genes from your query.

All faceted plots use `facet_grid(... , space = "free_y")` to ensure bars have uniform thickness across panels regardless of how many terms each source contributes.


## Advanced Configuration

### Restricting Annotation Sources

To query only specific databases, set the `sources` parameter:

```r
sources <- c("GO:BP", "GO:CC", "GO:MF")           # GO only
sources <- c("GO:BP", "KEGG", "REAC")              # GO:BP + pathways
sources <- c("GO:BP", "GO:MF", "KEGG", "REAC", "WP")  # common combination
```

### Changing the Number of Top Terms

Adjust `top_n` to show more or fewer terms per source in the plots:

```r
top_n <- 20    # show top 20 terms per source (default is 10)
```

### Switching Multiple-Testing Correction

```r
correction_method <- "fdr"          # Benjamini-Hochberg (widely used)
correction_method <- "bonferroni"   # most conservative
correction_method <- "g_SCS"        # g:Profiler default (recommended)
```

The g:SCS method is specifically designed for enrichment analysis where GO terms overlap hierarchically. Standard FDR and Bonferroni treat all tests as independent, which is not true for GO — parent and child terms share genes. g:SCS accounts for this dependency structure and is generally the recommended choice for g:Profiler analyses.

### Saving High-Resolution Figures

Set `eval=TRUE` in the **"Save Publication Figures"** chunk. Figures are saved at 300 DPI as both PNG and PDF. Adjust dimensions in the `ggsave()` calls to meet specific journal requirements.


## Interpreting Results

### Key Columns in the Results Table

| Column              | Meaning                                                        |
|---------------------|----------------------------------------------------------------|
| `source`            | Annotation database (GO:BP, KEGG, REAC, etc.)                 |
| `term_id`           | Database accession (e.g. GO:0006915, KEGG:04210)               |
| `term_name`         | Human-readable description                                     |
| `p_value`           | Adjusted p-value (after multiple-testing correction)           |
| `term_size`         | Total genes annotated to this term in the database             |
| `query_size`        | Number of your query genes that have any annotation            |
| `intersection_size` | Number of your query genes that overlap with this term         |
| `precision`         | intersection_size / query_size                                 |
| `recall`            | intersection_size / term_size                                  |
| `fold_enrich`       | How many times more enriched than expected by chance           |
| `gene_ratio`        | Same as precision; fraction of query genes in this term        |
| `intersection`      | Comma-separated list of overlapping gene IDs (if evcodes=TRUE) |

### What Makes a Result Biologically Meaningful?

A statistically significant result (low p-value) is not always biologically meaningful. Consider:

- **Fold enrichment**: A term with p = 0.01 but fold enrichment = 1.2 is barely enriched. Terms with fold enrichment > 2-3 are more likely to reflect genuine biological signal.
- **Intersection size**: A term with only 1-2 overlapping genes may be a statistical artifact. Terms with 3+ overlapping genes are more robust.
- **Term size**: Very broad terms (term_size > 1000, e.g. "protein binding") are less informative than specific terms (term_size 10-200).
- **Redundancy**: GO terms are hierarchical — a significant parent term (e.g. "apoptotic process") and its significant child term (e.g. "intrinsic apoptotic signaling pathway") are not independent findings. The g:SCS correction partially addresses this, but biological interpretation should still consider term relationships.


## Troubleshooting

| Problem | Solution |
|---------|----------|
| "No significant enrichment terms found" | Relax `user_threshold` (e.g. 0.1), check gene ID format, verify organism code, try removing the background file |
| "Gene file not found" | Check that the file path is correct relative to the `.Rmd` file location |
| Very few terms enriched | Your gene list may be too small (< 10 genes), or the organism has limited annotations |
| g:Profiler server timeout | Check your internet connection; the g:Profiler server may be temporarily down — try again later |
| `ggtext` or `tidytext` not found | Run `install.packages(c("ggtext", "tidytext"))` |
| Fold enrichment values seem implausibly high | This can happen for very small terms (term_size = 2-3); these are statistically fragile and should be interpreted cautiously |


## Citation

If you use this pipeline in a publication, cite the g:Profiler tool:

> Kolberg L, Raudvere U, Kuzmin I, Adler P, Vilo J, Peterson H (2023). "g:Profiler — interoperable web service for functional enrichment analysis and gene identifier mapping (2023 update)." *Nucleic Acids Research*, 51(W1), W199-W204. doi: [10.1093/nar/gkad347](https://doi.org/10.1093/nar/gkad347)

For the R package:

> Kolberg L, Raudvere U, Kuzmin I, Vilo J, Peterson H (2020). "gprofiler2 — an R package for gene list functional enrichment analysis and namespace conversion toolset g:Profiler." *F1000Research*, 9:709. doi: [10.12688/f1000research.24956.2](https://doi.org/10.12688/f1000research.24956.2)


## License

This pipeline script is provided as-is for academic research use.
