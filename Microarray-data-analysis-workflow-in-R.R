#---
#title: "GLDS-3"
#author: "Abdus Salam"
#date: "2/25/2022"
#---

# Loading libraries
library(fdrtool)
library(lattice)
library(Matrix)
library(rpart)
library(affy)
library(affydata)
library(ArrayExpress)
library(limma)
library(Biobase)
library(Biostrings)
library(genefilter)
library(affy)
library(affyPLM)
library(ggplot2)
library(stats)
library(annotate)
library(biomaRt)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(pathview)
library(ReactomePA)
library(clusterProfiler)
library(RColorBrewer)
library(pheatmap)
library(EnhancedVolcano)
library(org.Dm.eg.db)

#W e will start by defining path to the data directory.
celpath = "E:/Gene_Lab/GLDS3/"

# Reading all *.CEL (*.cel) files in your current working directory and storing them into the AffyBatch object expdata Where expdata is the new AffyBatch object that will be created at the end of the process operated by the function ReadAffy.
data = ReadAffy(celfile.path=celpath)

# Displaying the summary of dataset.
data

# Importing the phenotype data:
ph = data@phenoData
ph
ph$sample
pData(data)

# Retrieving probe annotation using affy
feat = featureNames(data)
feat[1:10]
length(featureNames(data))
length(probeNames(data))

# Lets add an additional column that will represent the characteristics of individual file. .
ph@data[ ,1] = c("Flight", "Flight", "Flight", "Ground", "Ground", "Ground")

# Visualizing the phenotype data:
ph@data

# Setting a path to custom directory.
setwd("E:/Gene_lab/GLDS1/")

# Displaying the distribution of intensities in each sample before normalization.
name = "boxplot_before_normalization.jpg"
jpeg(name)
boxplot(data,which='pm',col='red',names=ph@data$sample) 
dev.off()

# Creating MA plot for each array.
for (i in 1:6)
{
name = paste("MAplot",i,".jpg",sep="")
jpeg(name)
MAplot(data,which=i)
dev.off()
}

# Normalization of the datset using RMA approach.
data.rma <- rma(data)
data.matrix = exprs(data.rma)

# Displaying the distribution of intensities in each sample before normalization.
name = "boxplot_after_normalization.jpg"
jpeg(name)
boxplot(data.rma,which='pm',col='blue',names=ph@data$sample) 
dev.off()

# Creating MA plots after normalization.
for (i in 1:6)
{
name = paste("MAplotnorm",i,".jpg",sep="")
jpeg(name)
MAplot(data.rma,which=i)
dev.off()
}

# Displaying the summary of normalized dataset.
data.rma

# Saving expression values of probes in working directory.
write.exprs(data.rma, file="E:/Gene_lab/GLDS3/Expression_values.xls")

# PCA of nomralized dataset based on sample characteristics.
PCA <- prcomp(t(exprs(data.rma)))
color=c('red', 'red', 'red', 'blue', 'blue', 'blue')
plot(PCA$x[,1:2], col = color, main = "PCA plot for Flight vs. Ground")
legend("right", legend = c("Flight", "Ground"), col = c("red", "blue"), pch =1, cex=0.8)

# Creating groups based on sample characteristics (Flight vs. Ground)
ph@data[ ,2] = c("Flight", "Flight", "Flight", "Ground", "Ground", "Ground")
colnames(ph@data)[2]="source"
ph@data
groups = ph@data$source
f = factor(groups,levels=c("Flight", "Ground"))

# Designing diffrential matrix
design = model.matrix(~ 0 + f)
colnames(design) = c("Flight","Ground")

# Obtaining summary
summary(exprs(data.rma))

# Calculate median expression level
cutoff <- median(exprs(data.rma))

# TRUE or FALSE for whether each gene is "expressed" in each sample
is_expressed <- exprs(data.rma) > cutoff

#Identify genes expressed in more than 2 samples
keep <- rowSums(is_expressed) > 2

# Check how many genes are removed / retained.
table(keep)

# Subset to just those expressed genes
gse <- data.rma[keep,]

# Fiting a linear model for each gene in the expression dataset by using designing matrix.
fit <- lmFit(exprs(gse), design)
head(fit$coefficients)

# Designing contrast matrix for diffrential expression of Cancer probes against controls.
contrast.matrix = makeContrasts(Flight-Ground, levels=design)
fit2 = contrasts.fit(fit, contrast.matrix)

# Now we are going to calculate the differential expression by empirical Bayes moderated t-statistics:
fit2 <- eBayes(fit2)
topTable(fit2)

# How many genes are differentially-expressed overall?
decideTests(fit2)
table(decideTests(fit2))

# Assigning row names to ID.
full_results <- topTable(fit2, number=Inf)
full_results <- tibble::rownames_to_column(full_results,"ID")

# Plotting a volcanoplot for visualization.
EnhancedVolcano(full_results, lab = rownames(full_results), x = 'logFC', y = 'P.Value', 
                FCcutoff = 0.5, pCutoff = 0.05, legendPosition = "right", 
                col=c('#C0C0C0', '#1E90FF', '#FFD700', '#FF6347'), 
                legendLabels=c('Not sig.','log2FC','adj.P',
                               'adj.P & log2FC'), border = 'full', borderWidth = 0.5, 
                labCol = '#FF6347', selectLab = "NA", 
                legendLabSize = 10, labSize = 0.00, xlim = c(-4,4), ylim = c(0,10), 
                title = "Flight vs Ground")

# Annotating probes with an IDs 
ensembl = useMart("ensembl")
ensembl = useDataset("dmelanogaster_gene_ensembl",mart=ensembl)
IDs <- getBM(attributes=c('affy_drosophila_2','entrezgene_id'), 
             filters = 'affy_drosophila_2', values = full_results$ID, mart = ensembl)
IDs <- IDs[!duplicated(IDs$affy_drosophila_2), ]
full_results <- full_results[!duplicated(full_results$ID), ]
IDs <- IDs %>% rename(ID = "affy_drosophila_2")
full_results <- dplyr::inner_join(full_results, IDs, by="ID")
full_results <- full_results %>% filter(full_results$P.Value < 0.05)

# Saving the DEGs reuslts in a working directory.
write.csv(full_results, file="E:/Gene_Lab/GLDS3/full_results.csv", row.names = F, quote = F)

# Making the list ready for gene set enrciment analysis adn filtering for significantly expressed genes.
genelist <- cbind(full_results$entrezgene_id.x, full_results$logFC)
colnames(genelist) <- c("ID", "logFC")
genelist <- as.data.frame(genelist)
geneList <- as.numeric(genelist[,2])
names(geneList) <- as.character(genelist[,1])
geneList <- sort(geneList, decreasing = TRUE)
head(geneList)
gene <- names(geneList)[abs(geneList) > 1]

# GSEA using msigdbr
msigdbr_species()

# Get all fruitfly gene sets:
m_df <- msigdbr(species = "Drosophila melanogaster")
head(m_df, 2) %>% as.data.frame

# Gene set names (gs_name) and gene identifiers (entrez_gene). 
m_df_H_t2g <- msigdbr(species = "Drosophila melanogaster", category = "H") %>% dplyr::select(gs_name, entrez_gene)
head(m_df_H_t2g)

# View what collections there are
collections <- msigdbr_collections()

# Select (subset) the Hallmark collection
m_df_H <- m_df[m_df$gs_cat=="H",]
head(m_df_H)

# For enrichemnt analysis we need a dataframe with two columns
msigdbr_t2g = m_df_H %>% dplyr::select(gs_name, entrez_gene) %>% as.data.frame()

# Enrichment by hypergeometric test implemented in the enricher function
res_enricher <- enricher(gene = gene, TERM2GENE = msigdbr_t2g)

# Gene set enrichment analysis is implemented in the GSEA funcion
res_GSEA <- GSEA(geneList, TERM2GENE = msigdbr_t2g)

# Convert it into a dataFrame:
res_GSEA_df <- as.data.frame(res_GSEA)

#Visualization
dotplot(res_GSEA, showCategory=20, font=5)

# Convert gene ID to Symbol
edox <- setReadable(res_GSEA, 'org.Dm.eg.db', 'ENTREZID')

# Visualization of network and a heatmap
cnetplot(edox, foldChange=geneList)

# KEGG analysis
kegg_organism = "dme"
kk2 <- gseKEGG(geneList     = geneList,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

# To look at the columns:
head(kk2)
kk2[c(1:10),c(2,3:4, 6)]

# Visualisation of the top KEGG pathway
hsa04510 <- pathview(gene.data  = geneList,
                     pathway.id = "dme00190",
                     species    = "dme",
                     limit      = list(gene=max(abs(geneList)), cpd=1))

# Dot plot visualising the first KEGG result from above
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)

