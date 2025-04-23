# Subtype Classification Module

This folder contains code and documentation related to the classification of SCLC subtypes based on bulk DEG signatures.

## Description

The goal of this module is to assign subtype labels (SCLC-A, SCLC-N, SCLC-P) to single cells from the SCLC dataset, using gene signatures derived from the bulk RNA-seq data published in George et al.

Although the original method was written in **Python**, this module implements the same approach in **R**, using a Seurat `.rds` object for input.

Steps include:
- Scoring each cell using `AddModuleScore()` or `score_genes()` based on bulk DEG lists
- Selecting training cells for each subtype
- Running classification using a distance-based approach in diffusion map space (e.g., via Phenograph or custom methods)
- Saving subtype predictions in the `meta.data` slot of the Seurat object

## Note

This classification method is **adapted from published work** and only used for subtype labeling prior to independent downstream analysis. All subsequent analyses (e.g., clustering, DEGs, GSEA) are performed independently.
