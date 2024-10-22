# SCALT: automatic identification of cell types from single-cell RNA sequencing data
SCALT (Single Cell Annotation Likelihood Tool) is an innovative method which introduces a paradigm-shift for the analysis of scRNAseq data. In this approach, cells are annotated to a specific type at individual level, by using a simple but elegant method based on maximum likelihood, without the need for clustering, dimensionality reduction or manual annotation. SCALT leverages a collection of 471 lists of cell-type specific genes, constructed by extensive re-analysis of comprehensive and expert curated catalogues (HPA and DISCO).

For further details, we advise to read the paper available at DOI

As mentioned in the paper, SCALT is composed of other parallel utilities that allow to perform more complex analysis such as:

1. build cell type specific lists of genes in a deterministic fashion starting from a counts matrix and the correspoding annotation for each cell;
2. build cell type specific lists of genes starting from a count matrix and a gathering of user-defined cell type specific lists of genes making the use of an hypergeometric test.


## SCALT manual
The manual of SCALT is available at [SCALT documentation](https://scalt.readthedocs.io/en/latest/)
