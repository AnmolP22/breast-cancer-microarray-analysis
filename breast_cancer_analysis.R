# Breast Cancer Gene Expression Analysis
# Store the dataset name
dataset_id <-"GSE15852"
#Display the dataset name
print(dataset_id)
# Install BiocManager if it is not already installed
if(!requireNamespace("BiocManager",quietly=TRUE)){
  install.packages("BiocManager")
}
# Install tools for bioinformatics analysis
BiocManager::install(c("GEOquery","limma"))
requireNamespace("GEOquery", quietly = TRUE)
requireNamespace("limma",quietly=TRUE)
print(requireNamespace("limma", quietly=TRUE))
print(requireNamespace("GEOquery", quietly=TRUE))
gse <- getGEO("GSE15852", GSEMatrix = TRUE)
sample_info <- data.frame(
  titles = titles,
  patient_id = patient_id,
  tissue = ifelse(grepl("Normal", titles),"Normal", "Tumour")
)
expr_matrix <- exprs(gse[[1]])
sample_info$patient_id <-
  factor(sample_info$patient_id)

sample_info$tissue <-
  factor(sample_info$tissue,
                             levels =
           c("Normal", "Tumour"))
sample_info$patient_id <-
  factor(sample_info$patient_id)
design <- model.matrix(~ patient_id + tissue, data = sample_info)
fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit)

library(limma)
fit <- lmFit(expr_matrix, design)
fit <- eBayes(fit)
results <- topTable(fit, coef = "tissueTumour", number = Inf, adjust.method = "BH")
head(results)
png(
  filename = "Top_20_Differentially_Expressed_Genes.png",
  width = 1600,
  height = 1000,
  res = 150
)

par(mar = c(10, 5, 4, 2))

barplot(
  top_20_genes$logFC,
  names.arg = top_20_genes$GeneSymbol,
  las = 2,
  col = ifelse(top_20_genes$logFC > 0, "red", "blue"),
  main = "Top 20 Differentially Expressed Genes",
  ylab = "Log2 Fold Change",
  cex.names = 0.7
)

dev.off()
upregulated <- unique_genes[
  unique_genes$logFC > 0,
]

downregulated <- unique_genes[
  unique_genes$logFC < 0,
]
write.table(
  upregulated$GeneSymbol,
  "C:/users/anmol/OneDrive/Documents/Upregulated_Genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  downregulated$GeneSymbol,
  "C:/users/anmol/OneDrive/Documents/Downregulated_Genes.txt",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
up_DAVID <- read.csv(
  "C:/users/anmol/OneDrive/Documents/DAVIDChartReport_Upregulated_Genes_2026-11-15.csv"
)

down_DAVID <- read.csv(
  "C:/users/anmol/OneDrive/Documents/DAVIDChartReport_Downregulated_Genes_2026-11-15.csv"
)
up_kegg <- subset(
  up_DAVID,
  Category == "KEGG_PATHWAY" & Benjamini < 0.05
)

down_kegg <- subset(
  down_DAVID,
  Category == "KEGG_PATHWAY" & Benjamini < 0.05
)
png(
  "C:/users/anmol/OneDrive/Documents/KEGG_Pathway_Comparison.png",
  width = 1600,
  height = 1000,
  res = 150
)

par(mar = c(5, 5, 4, 2))

barplot(
  c(nrow(up_kegg), nrow(down_kegg)),
  names.arg = c("Upregulated", "Downregulated"),
  col = c("red", "blue"),
  main = "Significant KEGG Pathways by Expression Group",
  ylab = "Number of Significant Pathways",
  ylim = c(0, 6)
)

dev.off()
fit <- lmFit(expr_matrix, design)

fit <- eBayes(fit)
results <- topTable(
  fit,
  coef = "tissueTumour",
  number = Inf,
  adjust.method = "BH"
)
significant <- results[
  results$adj.P.Val < 0.05 &
    abs(results$logFC) >= 1,
]

gene_symbols <- AnnotationDbi::mapIds(
  hgu133a.db,
  keys = rownames(significant),
  column = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"
)

significant$GeneSymbol <- unname(gene_symbols)

significant_annotated <- significant[
  !is.na(significant$GeneSymbol) &
    significant$GeneSymbol != "",
]

upregulated <- significant_annotated[
  significant_annotated$logFC > 0,
]

downregulated <- significant_annotated[
  significant_annotated$logFC < 0,
]

write.csv(upregulated,
          "Upregulated_Genes.csv",
          row.names = TRUE)

write.csv(downregulated,
          "Downregulated_Genes.csv",
          row.names = TRUE)
up_unique <- upregulated[
  order(upregulated$adj.P.Val),
]

up_unique <- up_unique[
  !duplicated(up_unique$GeneSymbol),
]

down_unique <- downregulated[
  order(downregulated$adj.P.Val),
]

down_unique <- down_unique[
  !duplicated(down_unique$GeneSymbol),
]

write.table(
  up_unique$GeneSymbol,
  "Upregulated_Unique_Genes.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)

write.table(
  down_unique$GeneSymbol,
  "Downregulated_Unique_Genes.txt",
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)