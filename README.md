# RNA-seq-analysis
A streamlined RNA-Seq pipeline for DESeq2 differential expression, automated gene annotation (biotypes), and interactive visualization.
# RNA-Seq Analysis with DESeq2

A pipeline for identifying differentially expressed genes, adding gene annotations/biotypes, and generating interactive volcano plots.

## Requirements
Requires R with the following packages: `DESeq2`, `biomaRt`, `tidyverse`, and `plotly`.

## Quick Start
1. Load your counts matrix and sample metadata.
2. Run `deseq2_analysis.R` to perform the differential expression.
3. Open `plots/interactive_volcano_plot.html` to explore the significant genes and their biotypes.

## Outputs
* `significant_genes_annotated.csv` - Filtered DEGs with gene symbols and biotypes.
* `interactive_volcano_plot.html` - Interactive visualization.
