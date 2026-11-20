Differential Gene Expression Analysis of Breast Cancer

Project Overview

This project explores differences in gene expression between normal and tumour breast tissue using publicly available microarray data.

The analysis uses the GSE15852 dataset from the NCBI Gene Expression Omnibus (GEO).

Objectives

- Compare gene expression between normal and breast tumour samples.
- Identify significantly differentially expressed genes.
- Annotate microarray probes with gene symbols.
- Explore biological functions and pathways associated with selected genes.

Dataset

- Dataset: GSE15852
- Platform: Affymetrix Human Genome U133A (GPL96)
- Samples: 86
- Normal samples: 43
- Tumour samples: 43
- Study design: 43 matched normal-tumour pairs.

Tools and Packages

- R and RStudio
- GEOquery
- limma
- Bioconductor
- hgu133a.db
- AnnotationDbi
- DAVID Functional Annotation Tool

Methodology

1. Retrieved public breast cancer microarray data from GEO.
2. Organized sample information and identified matched normal-tumour pairs.
3. Prepared the expression matrix and applied log2 transformation.
4. Performed paired differential gene expression analysis using the limma package.
5. Annotated microarray probe IDs with gene symbols.
6. Identified upregulated and downregulated genes.
7. Explored functional enrichment using DAVID.

Key Results

- 416 significantly differentially expressed probes were identified.
- These corresponded to 347 unique gene symbols after annotation and deduplication.
- Significance criteria: adjusted p-value below 0.05 and absolute log2 fold change of at least 1.
- Upregulated and downregulated gene lists were generated for further interpretation.

Functional Enrichment

DAVID functional enrichment analyses were explored using the gene lists.

Note: The DAVID reports currently included in this repository are preliminary and were generated using earlier versions of the gene lists. They have not been updated using the final regenerated lists.

Repository Contents

- "breast_cancer_analysis.R" – R analysis script.
- "GSE15852_Differential_Expression.csv" – Differential expression results.
- "GSE15852_Significant_Genes.csv" – Significant results.
- "GSE15852_Unique_Significant_Genes.csv" – Unique annotated gene list.
- Upregulated and downregulated gene lists.
- Volcano plot and other visualization files.
- Preliminary DAVID enrichment reports.

Limitations

- This is an exploratory analysis of a public microarray dataset.
- Findings require further biological validation.
- Functional enrichment results are preliminary and should not be treated as final validated findings.

Learning Outcomes

This project provided practical exposure to public biological databases, microarray gene expression data, R-based statistical analysis, gene annotation and functional enrichment.

Acknowledgements

Dataset obtained from the NCBI Gene Expression Omnibus (GEO): GSE15852.

This is an undergraduate personal bioinformatics project.
