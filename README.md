# Diff_Promoter_Analysis
This pipeline performs differential binding analysis on ChIP-seq/ChIP-like seq (e.g. ATACseq, Cut&amp;Run seq) data across multiple histone marks, identifying differentially regulated promoters and visualizing results using genome browser tracks.

## Scripts

### 1. `diff_bind_3.R` - Differential Binding Analysis
Performs LIMMA-based differential analysis on promoter regions for multiple histone marks.

**Key Functions:**
- `run_promoter_limma()` - Runs differential analysis for a specific histone mark
- `annotate_all_contrasts()` - Annotates significant results with gene symbols

**Analysis Steps:**
1. Loads promoter count data (3kb regions)
2. Filters low-expression promoters (mean > 2)
3. Log10 transformation
4. LIMMA differential analysis (treatment vs Control)
5. Exports UP/DOWN regulated genes for each contrast

**Histone Marks Analyzed:**
- A1, A2, P, R1, R2

**Treatment Groups:**
- GAL1, G3BP, POSTN (all compared to Ctrl)

### 2. `diff_bind_4.R` - Gene Visualization
Identifies top differentially regulated genes and creates genome browser visualizations.

**Key Functions:**
- `get_top_genes_for_mark()` - Extracts genes with highest/lowest logFC
- `plot_gene_gviz()` - Creates Gviz tracks showing ChIP-seq signal

**Visualization Features:**
- Genome tracks for control and treatment samples
- Gene models with exon structure
- BigWig signal overlays (control in red, treatment in purple)
- ±3kb promoter context

## Input Files

### Required Data
- `hg38_promoters_chr1-22_X_3kb.tsv` - Promoter count matrix
- `sample_name_replacement.csv` - Sample name mapping
- `hg38_promoters_chr1-22_X_3kb.bed` - Promoter coordinates
- `geneModels.csv` - Gene annotation
- BigWig files in `./data/bw/`

### Sample Naming Convention
```
{GROUP}_{HISTONE}_{REP}
Example: Ctrl_A1_r1, GAL1_A2_r2, G3BP_P_r3
```

## Output Structure

```
./
├── data/
│   ├── promoter/
│   │   ├── A1/
│   │   │   ├── GAL1_Ctrl_UP_genes.csv
│   │   │   ├── GAL1_Ctrl_DOWN_genes.csv
│   │   │   └── bw_limma_*.csv
│   │   ├── A2/
│   │   ├── P/
│   │   ├── R1/
│   │   └── R2/
│   └── pr_metadata.csv
└── analysis/
    └── promoter/
        ├── A1/
        │   ├── hist_raw_values.pdf
        │   ├── hist_log10_values.pdf
        │   ├── top_genes_summary_A1.csv
        │   └── Gviz_*.pdf
        └── [A2, P, R1, R2]/
```

## Key Parameters

### Differential Analysis
- **Low expression filter:** mean count > 2
- **Normalization:** log10(x + 1)
- **Significance cutoff:** P-value < 0.05
- **Log fold-change cutoff:** |logFC| > 0

### Visualization
- **Promoter window:** ±3000 bp from gene boundaries
- **Y-axis:** Auto-scaled to maximum signal
- **Control color:** Red
- **Treatment color:** Purple

## Dependencies

```r
DESeq2
limma
stringr
ChIPseeker
TxDb.Hsapiens.UCSC.hg38.knownGene
org.Hs.eg.db
rtracklayer
Gviz
AnnotationHub
```

## Usage

### Run Complete Pipeline
```r
# 1. Differential analysis
source("diff_bind_3.R")

# 2. Generate visualizations
source("diff_bind_4.R")
```

### Analyze Single Histone Mark
```r
# Run LIMMA
res <- run_promoter_limma(df, histoneMark = "A1")

# Annotate results
summary <- annotate_all_contrasts(res, histoneMark = "A1")
```

### Plot Specific Gene
```r
plot_gene_gviz(
  histoneMark = "A1",
  gene = "MYC",
  sample_type = "G3BP",
  replacement_file = "./data/sample_name_replacement.csv",
  geneModels = geneModels,
  bw_dir = "./data/bw/"
)
```

## Output Interpretation

### UP/DOWN Gene Lists
- **logFC:** Log2 fold-change (treatment/control)
- **P.Value:** Unadjusted p-value
- **SYMBOL:** Gene symbol
- **chr, start, end:** Genomic coordinates

### Gviz Plots
- **Top track:** Chromosome ideogram
- **Second track:** Genomic coordinates
- **Third track:** Gene model with exons
- **Bottom tracks:** ChIP-seq signal (one per replicate)

## Notes

- Script uses custom library path for package installation
- All contrasts use Control as reference group
- Gene models exclude exon 1 in visualizations
- Missing genes in annotation are skipped automatically