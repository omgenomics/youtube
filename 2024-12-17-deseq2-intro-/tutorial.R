# ==============================================================================
# RNA-Seq differential expression tutorial with DESeq2
# GXA Project: https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5244/Results
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup
# ------------------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

# ------------------------------------------------------------------------------
# Download data
# ------------------------------------------------------------------------------

# Download raw counts
counts = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/DifferentialSecondaryDataFiles.RnaSeq/raw-counts")
head(counts)

# Download metadata
metadata = read.delim("https://www.ebi.ac.uk/gxa/experiments-content/E-MTAB-5244/resources/ExperimentDesignFile.RnaSeq/experiment-design")
head(metadata)


# ------------------------------------------------------------------------------
# Wrangle data for DESeq2
# ------------------------------------------------------------------------------

# DESeq expects the counts to have gene IDs as row names
head(counts)
rownames(counts) = counts$Gene.ID
head(counts)

# Remove unused columns (gene ID and gene name)
genes = counts[, c("Gene.ID", "Gene.Name")]
counts = counts[, -c(1, 2)]
head(counts)

# DESeq expects the metadata matrix to have sample IDs in the rownames
head(metadata)
rownames(metadata) = metadata$Run
head(metadata)

# Only keep columns of interest
metadata = metadata[, c("Sample.Characteristic.genotype."), drop=FALSE]
# Look at metadata to see how the variables change with respect to each other
metadata

# Rename column
colnames(metadata) = c("genotype")
metadata

# Remove spaces in names to avoid DESeq warnings
metadata$genotype[metadata$genotype == 'wild type genotype'] = 'wildtype'
metadata$genotype[metadata$genotype == 'Snai1 knockout'] = 'knockout'
metadata

# Turn genotype into a factor
metadata$genotype = factor(metadata$genotype, levels=c("wildtype", "knockout"))
metadata$genotype


# ------------------------------------------------------------------------------
# Spot check expression for knockout gene SNAI1
# ------------------------------------------------------------------------------

gene_id = genes$Gene.ID[ genes$Gene.Name == 'SNAI1' ]
gene_counts = counts[gene_id, ]
gene_counts

gene_data = cbind(metadata, counts=as.numeric(gene_counts))
gene_data

library(ggplot2)
ggplot(gene_data, aes(x = genotype, y = counts, fill = genotype)) + geom_boxplot()


# ------------------------------------------------------------------------------
# Run DESeq
# ------------------------------------------------------------------------------

dds <- DESeqDataSetFromMatrix(countData=counts, colData=metadata, design=~genotype)

# Ignore genes with low counts
dds <- dds[rowSums(counts(dds)) > 10, ]

# Run DESeq
dds <- DESeq(dds)

# Compare expression
res = results(dds, contrast=c("genotype", "knockout", "wildtype"), alpha=1e-5)
res

# Other examples of design formula:
# - https://www.ebi.ac.uk/gxa/experiments/E-MTAB-8031/Results
#     * RNA-seq on lung tumors
#     * Variables:
#         1. Location: tumor cells or normal tissue nearby
#         2. Sex of the patient
#     * Formula = ~sex + location
#
# - https://www.ebi.ac.uk/gxa/experiments/E-MTAB-10041/Results
#     * RNA-seq on mouse colon cells
#     * Variables:
#         1. Compound: cancer drug or not
#         2. Genotype: wild type mouse, and 2 other genotypes of TP53 mutations
#     * Formula = ~genotype + compound  doesn't make sense since want to compare each genotype separately
#     * Formula = ~group
#         metadata$group = factor(paste(metadata$genotype, metadata$compound))


# Sidenote: "~" is not a DESeq specific operator
head(iris)
model = lm(Petal.Width ~ Petal.Length, iris)
plot(iris$Petal.Length, iris$Petal.Width)
abline(model)


# ------------------------------------------------------------------------------
# Spot checks
# ------------------------------------------------------------------------------

# Compare results to GXA
# https://www.ebi.ac.uk/gxa/experiments/E-MTAB-5244/Results?specific=true&geneQuery=%255B%255D&filterFactors=%257B%257D&cutoff=%257B%2522foldChange%2522%253A11%252C%2522pValue%2522%253A0.01%257D&regulation=%2522UP_DOWN%2522

# Merge gene name into data frame so can compare to GXA UI using gene names
res_df = as.data.frame(res)
head(res_df)
head(genes)
res_df = merge(res_df, genes, by='row.names')
head(res_df)

genes_to_check = c("THY1", "SFMBT2", "PASD1", "SNAI1")
res_df[res_df$Gene.Name %in% genes_to_check, ]


# ------------------------------------------------------------------------------
# Visualization
# ------------------------------------------------------------------------------

# MA plot
plotMA(res)

# Volcano plot
BiocManager::install('EnhancedVolcano')
library(EnhancedVolcano)
EnhancedVolcano(res, lab=rownames(res), x='log2FoldChange', y='pvalue')

# Circos plot
BiocManager::install("biomaRt")
library(biomaRt)

# Find dataset name in Ensembl
ensembl <- useEnsembl(biomart="genes")
datasets = listDatasets(ensembl)
head(datasets)

dataset_nb = grep("human", datasets$description, ignore.case=TRUE)
dataset_nb

dataset = datasets$dataset[dataset_nb]
dataset

ensembl <- useDataset(dataset=dataset, mart=ensembl)

# Get coordinates of all genes
attributes <- c("ensembl_gene_id", "chromosome_name", "start_position", "end_position")
values <- list(ensembl_gene_id=c())
all.genes <- getBM(attributes=attributes, values=values, mart=ensembl)
head(all.genes)

# Rename column so it matches res_df
head(res_df)
colnames(all.genes)[1] = "Gene.ID"
head(all.genes)

# Merge the DESeq output with the Ensembl gene coordinates
merged_data <- merge(all.genes, res_df, by="Gene.ID")
head(merged_data)

# Add chr prefix to chromosome names
merged_data$chromosome_name <- paste("chr", merged_data$chromosome_name, sep = "")
head(merged_data)
write.csv(merged_data, "~/Documents/rna-seq/deseq.csv", row.names = FALSE)

# Same for subset
merged_data_subset = merged_data[merged_data$Gene.Name %in% genes_to_check, ]
head(merged_data_subset)
write.csv(merged_data_subset, "~/Documents/rna-seq/deseq_subset.csv", row.names = FALSE)
