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

## Data and Method Citation

- The subtype-specific gene signatures (SCLC-A, N, P) used in this project are based on the study by George et al., 2015:
  > George, J., et al. *Comprehensive genomic profiles of small cell lung cancer*. Nature 524, 47–53 (2015). [https://doi.org/10.1038/nature14664](https://www.nature.com/articles/nature14664)

- The SCLC subtype classification workflow used in this project is adapted from:

> Chan, J.M., Quintanal-Villalonga, Á., Gao, V.R., Xie, Y., Allaj, V., Chaudhary, O., ... Rudin, C.M.  
> *Signatures of plasticity, metastasis, and immunosuppression in an atlas of human small cell lung cancer*.  
> Cancer Cell, 2021. [https://doi.org/10.1016/j.ccell.2021.04.005](https://www.sciencedirect.com/science/article/pii/S1535610821004979?via%3Dihub#app2)

The authors provided a Python-based classification pipeline (`Classify_SCLC_subtype.py`), which was reimplemented here using R and integrated into a Seurat-based analysis workflow.
