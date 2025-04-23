# Immune Suppression DEG Analysis in SCLC

This repository contains the full R-based workflow for identifying and analyzing differentially expressed genes (DEGs) related to immune suppression in Small Cell Lung Cancer (SCLC) subtypes using single-cell RNA-seq data.

## Structure

- `immune_suppression_analysis.Rmd`: Main reproducible R Markdown script.
- `results/`: Output files including filtered genes, enrichment results, and plots.

## Data Sources & References

- MSigDB C7 Gene Sets (Immunologic Signatures)
  Gene sets retrieved using the msigdbr R package (v1.10.0) from the C7: Immunologic Signatures collection of MSigDB.
  [https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp#C7](https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp?collection=IMMUNESIGDB)

- Software Package
  Dudoit, S., & Tan, Y. (2021). msigdbr: MSigDB Gene Sets for Multiple Organisms in a Tidy Data Format.
  Bioconductor. [https://bioconductor.org/packages/msigdbr](https://bioconductor.org/packages/release/data/experiment/html/msigdb.html)
