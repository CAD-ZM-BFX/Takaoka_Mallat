#!/usr/local/bin/Rscript
# SLX-22500 (GTC280), NRP1 KO/WT with cWD and iWD
# Mouse 
# each type has 6 replicates.
# CAD_mt709_0001; Ensemble GRCm39
# Information about your samples related to bioinformatics:
# 1. These were Illumina TruSeq Stranded mRNA libraries.
# 2. PE50 was sequenced on the Illumina NovaSeq6000 platform.
# 3. The average library size (including adapters) was 377 bp.
# 4. IDT for Illumina TruSeq RNA UD indexes - 8,8 indexes. 
#    The incorrect index type was IDT for Illumina DNA/RNA UD Indexes - 10,10 indexes
# Marcella mentioned that sample57,58,64 and 66 Index can not find in the pool.
# Sherin demultiplexed the lost reads pooled data and the above 4 samples have less number of reads
# 57-4.8Mb (cWD_WT_R1), 58-1.55Gb (cWD_WT_R2), 64-351.7Mb(cWD_KO_R2), 66-880.4Mb(cWD_KO_R4)
# Link to publication
# TO ADD ONCE AVAILABLE
#
#
# Script available from:
# https://github.com/CAD-BFX/CAD_mt709_0001
#
#
# Analysis Performed by Xiaohui Zhao
# Department of medicine, University of Cambridge, UK
# Copyright Xiaohui Zhao (xz289@cam.ac.uk)
#
#------------------------------------------------------------------------------
# License Information
# This program is free software: you can redistribute it and/or modify  
# it under the terms of the GNU General Public License as published by  
# the Free Software Foundation, version 3.
#
# This program is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#------------------------------------------------------------------------------

message("+-------------------------------------------------------------------------------+")
message("+  General settings and libraries calling for DESeq2 analysis                   +")
message("+-------------------------------------------------------------------------------+")

suppressPackageStartupMessages({
  library('DESeq2')
  library('ggplot2')
  library('RColorBrewer')
  library("cowplot")
  library("reshape")
  library("pheatmap")
  library("ggrepel")
  library("reshape2")
  library("biomaRt")
  library("matrixStats")
  library("plyr")
  library("BiocParallel")
  library("dplyr")
  library("ggalt")
  library("limma")
  library("apeglm")
  library("gdata") 
  library("ComplexHeatmap")
  library("seriation")
  library("methods")
  library("utils")
  library("Matrix")
  library("useful")
  library("edgeR")
  library("GetoptLong")
  library("UpSetR")
  library("circlize")
  library("VennDiagram")
  library("eulerr")
  library("readr")
  library("sva")
  library("gridExtra")
  library("grid")
  library("hexbin")
  ## GO package
  library("clusterProfiler")
  library("DOSE")
  library("GSEABase")
  library("AnnotationHub")
  library("org.Mm.eg.db")
  library("gage")
  library("gageData")
  library("enrichplot")
  library("ggraph")
  library("ggforce")
  ## secreted and tissueEnrich
  library("TissueEnrich")
  library("UniProt.ws")
  library("readxl")
  library("corrplot") 
  library("clusterSim") 
  ## scRNASeq
  library("scran")
  library("scater")
  library("SingleCellExperiment")
  library("bigmemory")
  library("mltools")
  library("rhdf5")
  library("recommenderlab")
  library("Seurat")
})


register(MulticoreParam(2))
basedir     <- "/Users/xz289/Documents/Minoru/Project_mt709_0001"
setwd(basedir)

Project         <- "SLX-22500"
Tissues         <- c("Nrp1")
significance    <- 0.05
l2fc            <- 1 
elementTextSize <- 10
TOPNUM          <- 2000

source("./Scripts/Utility_Function.R")


message("+-----             Load the data and ensemble                     -----------+")

load("./Data/References_Data/GRCm39/Ensembl_mmusculus_ID_Name_Des_Chr_GRCm39.RData")
load("./Data/SLX-22500/deseq2.dds.RData")

cts.batch1                   <- assay(dds)
cts.batch1                   <- cts.batch1[,order(colnames(cts.batch1))]
samT.batch1                  <- read.csv("./Data/SLX-22500/CAD_mt709_0001_NRP1-SampleTable.csv", header=T)
samT.batch1$Treatment        <- as.factor(samT.batch1$Treatment)
samT.batch1$Treatment        <- relevel(samT.batch1$Treatment, "WT")

samT.batch1                  <- samT.batch1[order(samT.batch1$sample),]
cmat <- cts.batch1

message("+--- Generate correlation heatmap for all samples ----------------+")
rownames(samT.batch1) <-samT.batch1$sample 
dds   <- DESeqDataSetFromMatrix(countData=cmat, colData=samT.batch1, design=~Condition+Treatment)
dds   <- DESeq(dds, parallel=TRUE)
dds   <- estimateSizeFactors(dds)

vsd <- vst(dds,     blind=F)
colData(vsd)

sampleDists      <- log(dist(t(assay(vsd))))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd$sample
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(10, "RdBu")) )(255)

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
## hclust
rv <- rowVars(assay(vsd))
o <- order(rv,decreasing=TRUE)
dists <- dist(t(assay(vsd)[head(o,2000),]))
hc <- hclust(dists)
allhc <- plot(hc, labels=vsd$sample)

## pca

pcaData    <- plotPCA(vsd, ntop=TOPNUM, intgroup=c("Group2"), returnData=TRUE)
pcaData$samplename <- colData(vsd)$sample
pcaData$condition  <- colData(vsd)$Condition
pcaData$treatment  <- colData(vsd)$Treatment

percentVar <- round(100 * attr(pcaData, "percentVar"))

pdf(paste0("./Figures_Tables/PCA_Figures/", Project, "-NRP1_all_Fig.PCA.T2000.pdf"),width=10,height=10)
par(bg=NA)
plt.all <-  ggplot(pcaData, aes(PC1, PC2, color=Group2)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(Group2), group=paste0(Group2), 
                        label=Group2), alpha=0.1, label.fontsize =10, label.buffer = unit(1, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("PCA (Top ", TOPNUM, " Most Variable Genes)"),
       caption = "CAD_mt709_0001: Minoru Takaoka" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1) 
plt_cWD <-  ggplot(pcaData[pcaData$condition=="cWD",], aes(PC1, PC2, color=treatment, label=samplename)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(treatment), group=paste0(treatment), 
                        label=treatment), alpha=0.1, label.fontsize =10, label.buffer = unit(1, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  geom_text_repel(aes(label=samplename), show.legend = FALSE, size=2)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("cWD PCA"),
       caption = "" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1) 

plt_iWD <-  ggplot(pcaData[pcaData$condition=="iWD",], aes(PC1, PC2, color=treatment, label=samplename)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(treatment), group=paste0(treatment), 
                        label=treatment), alpha=0.1, label.fontsize =10, label.buffer = unit(1, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  geom_text_repel(aes(label=samplename), show.legend = FALSE, size=2)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("iWD: PCA"),
       caption = "" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1)  
bottom_row <- plot_grid(plt_cWD, plt_iWD, labels = c('B', 'C'), align = 'h', rel_widths = c(1, 1))

plot_grid(plt.all, bottom_row, labels = c('A', ''), ncol = 1, align = 'vh', rel_heights = c(1, 1))
dev.off()

message("+--- remove 2 samples with lowest number of reads, cWD_WT_REP1 and cWD_KO_REP4---------------------+")

dds0   <- DESeqDataSetFromMatrix(countData=cmat[,-c(4,7)], colData=samT.batch1[-c(4,7),], design=~Condition+Treatment)
dds0   <- DESeq(dds0, parallel=TRUE)
dds0   <- estimateSizeFactors(dds0)

vsd0 <- vst(dds0,     blind=F)
colData(vsd0)

## heatmap clustering
sampleDists      <- log(dist(t(assay(vsd0))))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- vsd0$sample
colnames(sampleDistMatrix) <- NULL

pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
## hclust
rv0 <- rowVars(assay(vsd0))
o   <- order(rv0,decreasing=TRUE)
dists0 <- dist(t(assay(vsd0)[head(o,2000),]))
hc0    <- hclust(dists0)
rm2hc  <- plot(hc0, labels=vsd0$sample)

## pca plot
pcaData0    <- plotPCA(vsd0, ntop=TOPNUM, intgroup=c("Group2"), returnData=TRUE)
pcaData0$samplename <- colData(vsd0)$sample
pcaData0$condition  <- colData(vsd0)$Condition
pcaData0$treatment  <- colData(vsd0)$Treatment

percentVar0 <- round(100 * attr(pcaData0, "percentVar"))

pdf(paste0("./Figures_Tables/PCA_Figures/", Project, "-NRP1_rm2_Fig.PCA.T2000.pdf"),width=10,height=10)
par(bg=NA)
plt.all0 <-  ggplot(pcaData0, aes(PC1, PC2, color=Group2)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(Group2), group=paste0(Group2), 
                        label=Group2), alpha=0.1, label.fontsize =10, label.buffer = unit(1, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  xlab(paste0("PC1: ",percentVar0[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar0[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("PCA (Top ", TOPNUM, " Most Variable Genes)"),
       caption = "CAD_mt709_0001: Minoru Takaoka" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1) 
plt_cWD0 <-  ggplot(pcaData[pcaData0$condition=="cWD",], aes(PC1, PC2, color=treatment, label=samplename)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(treatment), group=paste0(treatment), 
                        label=treatment), alpha=0.1, label.fontsize =10, label.buffer = unit(1, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  geom_text_repel(aes(label=samplename), show.legend = FALSE, size=2)+
  xlab(paste0("PC1: ",percentVar0[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar0[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("cWD rm2 PCA"),
       caption = "" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1) 

plt_iWD0 <-  ggplot(pcaData[pcaData0$condition=="iWD",], aes(PC1, PC2, color=treatment, label=samplename)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(treatment), group=paste0(treatment), 
                        label=treatment), alpha=0.1, label.fontsize =10, label.buffer = unit(1, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  geom_text_repel(aes(label=samplename), show.legend = FALSE, size=2)+
  xlab(paste0("PC1: ",percentVar0[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar0[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("iWD rm2 PCA"),
       caption = "" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1)  
bottom_row0 <- plot_grid(plt_cWD0, plt_iWD0, labels = c('B', 'C'), align = 'h', rel_widths = c(1, 1))

plot_grid(plt.all0, bottom_row0, labels = c('A', ''), ncol = 1, align = 'v', rel_heights = c(1, 1))
dev.off()



message("+--- remove 4 samples with low number of reads---------------------+")

## Remove the six samples ("cWD_WT_REP1", "cWD_WT_REP2", "cWD_KO_REP2", "cWD_KO_REP4", "iWD_KO_REP5", "iWD_KO_REP6)
cmat1 <- cts.batch1[, -c(2,4,7,8)]

## Run DESeq2
rownames(samT.batch1) <-samT.batch1$sample 
dds1   <- DESeqDataSetFromMatrix(countData=cmat1, colData=samT.batch1[-c(2,4,7,8),], design=~Condition+Treatment)
dds1   <- DESeq(dds1, parallel=TRUE)
dds1   <- estimateSizeFactors(dds1)

vsd1 <- vst(dds1,     blind=F)
colData(vsd1)


## hclust
rv1 <- rowVars(assay(vsd1))
o   <- order(rv1,decreasing=TRUE)
dists1 <- dist(t(assay(vsd1)[head(o,2000),]))
hc1    <- hclust(dists1)
rm4hc  <- plot(hc1, labels=vsd1$sample)

## pca

pcaData1    <- plotPCA(vsd1, ntop=TOPNUM, intgroup=c("Group2"), returnData=TRUE)
pcaData1$samplename <- colData(vsd1)$sample
pcaData1$condition  <- colData(vsd1)$Condition
pcaData1$treatment  <- colData(vsd1)$Treatment

percentVar1 <- round(100 * attr(pcaData1, "percentVar"))

pdf(paste0("./Figures_Tables/PCA_Figures/", Project, "-NRP1_rm4_Fig.PCA.T2000.pdf"),width=10,height=10)
par(bg=NA)
plt.all1 <-  ggplot(pcaData1, aes(PC1, PC2, color=Group2)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(Group2), group=paste0(Group2), 
                        label=Group2), alpha=0.1, label.fontsize =10, label.buffer = unit(1, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  xlab(paste0("PC1: ",percentVar1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar1[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("PCA (Top ", TOPNUM, " Most Variable Genes)"),
       caption = "CAD_mt709_0001: Minoru Takaoka" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1) 
plt_rm4_cWD <-  ggplot(pcaData1[pcaData1$condition=="cWD",], aes(PC1, PC2, color=treatment, label=samplename)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(treatment), group=paste0(treatment), 
                        label=treatment), alpha=0.1, label.fontsize =10, label.buffer = unit(1, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  geom_text_repel(aes(label=samplename), show.legend = FALSE, size=2)+
  xlab(paste0("PC1: ",percentVar1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar1[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("Rm4: PCA"),
       caption = "" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1) 

plt_rm4_iWD <-  ggplot(pcaData1[pcaData1$condition=="iWD",], aes(PC1, PC2, color=treatment, label=samplename)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(treatment), group=paste0(treatment), 
                        label=treatment), alpha=0.1, label.fontsize =10, label.buffer = unit(1, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  geom_text_repel(aes(label=samplename), show.legend = FALSE, size=2)+
  xlab(paste0("PC1: ",percentVar1[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar1[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("Rm4: PCA"),
       caption = "" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1)  
bottom_row <- plot_grid(plt_rm4_cWD, plt_rm4_iWD, labels = c('B', 'C'), align = 'h', rel_widths = c(1, 1))

plot_grid(plt.all1, bottom_row, labels = c('A', ''), ncol = 1, align = 'v', rel_heights = c(1, 1))
dev.off()
message("+-------------------------------------------------------------------------------+")
message("+---       Pairwise DESeq Analysis, KO vs WT                               -----+")
message("+-------------------------------------------------------------------------------+")

samT.rm4       <- as.data.frame(samT.batch1[-c(2,4,7,8),], stringsAsFactors = TRUE)
control.names  <- c("cWD_WT", "iWD_WT")
compare.names  <- c("cWD_KO", "iWD_KO")
selGroups      <- c("cWD", "iWD")

for(i in 1:2){
  print(i)
  samT.sub     <- subset(samT.rm4, Condition==selGroups[i])
  samT.sub$Treatment  <- relevel(samT.sub$Treatment, "WT")
  mat.sub      <- cmat1[, colnames(cmat1)%in%samT.sub$sample]
  colnames(mat.sub)  <- as.factor(colnames(mat.sub))
  dds.sub      <- DESeqDataSetFromMatrix(countData = mat.sub,
                                         colData = samT.sub,
                                         design= ~ Treatment)
  dds.sub      <- estimateSizeFactors(dds.sub)
  sizeFactors(dds.sub)
  
  dds.sub      <- DESeq(dds.sub, parallel=TRUE)
  vst.sub      <- vst(dds.sub,     blind=F)
  resultsNames(dds.sub)
  colData(vst.sub)
  save(mat.sub, samT.sub, dds.sub, vst.sub, 
       file = paste0("./Data/SLX-22500/Rm4_",selGroups[i], "_DESeq2_dds_vstN.RData"))
  gc()
}

customPCA     <- function(sampleTBL, RLD, TOPNUM, ensEMBL2id) {
  RLD                 <- as.data.frame(RLD) 
  RLD.ori             <- RLD
  RLD.ori$ensembl_gene_id <- rownames(RLD.ori)
  RLD.mer             <- merge(RLD.ori, ensEMBL2id, by="ensembl_gene_id") 
  RLD.mer             <- RLD.mer[-which(RLD.mer$external_gene_name==""),]## remove no gene names
  RLD.mer             <- RLD.mer[-(which(duplicated(RLD.mer$external_gene_name)==T)),] # remove duplicated genes
  RLD.new             <- RLD.mer[,colnames(RLD)]
  rownames(RLD.new)   <- RLD.mer$external_gene_name
  colnames(RLD.new)   <- sampleTBL$sample
  rv     <- rowVars(as.matrix(RLD.new))
  select <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca    <- prcomp(t(RLD.new[select, ]))
  
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  
  scores    <- data.frame(sample=sampleTBL$sample, pca$x, 
                          condition=sampleTBL$Treatment)
  
  plt.pca <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, label=sample) ) +
    geom_point(size = 3, alpha=0.75) + 
    geom_text_repel(aes(label=sample), show.legend = FALSE, size=2) +
    xlab(pc1lab) + ylab(pc2lab) + 
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +             
    theme(text = element_text(size=elementTextSize)) + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  plt.pca.nl <- ggplot(scores, aes(x = PC1, y = PC2, colour=condition, fill=condition, label=sample) ) +
    geom_point(size = 3, alpha=0.75 ) + 
    geom_encircle(alpha = 0.1, show.legend = FALSE, aes(colour=condition, fill=condition, group = condition)) +
    xlab(pc1lab) + ylab(pc2lab) + 
    ggtitle(paste0("PCA Top ", TOPNUM, " MV")) +
    theme(text = element_text(size=elementTextSize))+ 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  
  
  loadings                 <- as.data.frame(pca$rotation)
  
  
  pca.1         <- loadings[order(loadings$PC1,decreasing=TRUE), ]
  pca.1.25      <- pca.1[c(1:25),]
  pca.1.25.plot <- ggplot(data=pca.1.25, aes(x=factor(rownames(pca.1.25),levels=unique(rownames(pca.1.25))), y=PC1)) + 
    geom_point(size = 3 ) + xlab("") + 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  
  pca.2         <- loadings[order(loadings$PC2,decreasing=TRUE), ]
  pca.2.25      <- pca.2[c(1:25),]
  pca.2.25.plot <- ggplot(data=pca.2.25, aes(x=factor(rownames(pca.2.25),levels=unique(rownames(pca.2.25))), y=PC2)) + 
    geom_point(size = 3 ) + xlab("") + 
    theme()+ 
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          axis.text.x = element_text(angle = 90, hjust = 1, size=8))
  
  return(list(plt.pca,plt.pca.nl,pca.1.25.plot, pca.2.25.plot))
  
}
for(i in 1:2){
  load(paste0("./Data/SLX-22500/Rm4_",selGroups[i], "_DESeq2_dds_vstN.RData"))
  pdf(paste0("./Figures_tables/PCA_Figures/DESeq2_pairwisePCA_Rm4_", compare.names[i], "_", control.names[i], "_N.pdf"), width=8, height= 6)
  pca <- customPCA(samT.sub, assay(vst.sub), TOPNUM, ensEMBL2id) 
  plot_grid(pca[[1]], pca[[2]], pca[[3]], pca[[4]], nrow=2, byrow=T)
  dev.off()
}


DES.files <- paste0("./Data/SLX-22500/Rm4_",selGroups, "_DESeq2_dds_vstN.RData")
RES.files  <- paste0("./Data/SLX-22500/Rm4_",selGroups,  "_res_all_merE_summaryN.csv")
RES.files1 <- paste0("./Data/SLX-22500/Rm4_",selGroups,  "_res_all_lfcshrink_merE_summaryN.csv")
RESsig.files  <- paste0("./Data/SLX-22500/Rm4_",selGroups,  "_res_Sig_merE_summaryN.csv")
RESsig.files1 <- paste0("./Data/SLX-22500/Rm4_",selGroups,  "_res_Sig_lfcshrink_merE_summaryN.csv")

merECRes_function <- function(resfile, ensEMBL2id, ddsfile,outfile){
  resdat   <- as.data.frame(resfile)
  resdat$ensembl_gene_id <- rownames(resdat)
  ### remove the padj missing
  resdat   <- subset(resdat, !is.na(padj)) 
  resmerE  <- merge(resdat, ensEMBL2id, by = "ensembl_gene_id")
  
  normCnt  <- as.data.frame(counts(ddsfile, normalized=TRUE))
  normCnt$ensembl_gene_id <- rownames(normCnt)
  resmerEC <- merge(resmerE, normCnt, by = "ensembl_gene_id")
  write.csv(resmerEC, file=outfile, row.names=F)
}


for(i in 1:2){
  load(DES.files[i])
  res       <- results(dds.sub, contrast = c("Treatment", "KO", "WT"))
  res       <- subset(res, !is.na(padj))
  res.sig   <- subset(res, padj <= significance & abs(log2FoldChange) >= l2fc)
  print(nrow(res.sig)); print(table(sign(res.sig$log2FoldChange)))
  write.csv(res.sig, file = RESsig.files[i])
  res.merEC <- merECRes_function(res, ensEMBL2id, dds.sub,RES.files[i])
  res1      <- lfcShrink(dds.sub, coef=2)
  res1.sig  <- subset(res1, padj <= significance & abs(log2FoldChange) >= l2fc)
  print(nrow(res1.sig)); print(table(sign(res1.sig$log2FoldChange)))
  write.csv(res1.sig, file = RESsig.files1[i])
  res.merEC1 <- merECRes_function(res1, ensEMBL2id, dds.sub, RES.files1[i])
}


## cWD, 80(60-, 20+) res; lfcshrink 57 (-1, 43; 1 14);
## iWD, 131(66-, 65+) res; lfcshrink 104 (-1, 48; 1 56);
message("+--- remove 6 samples with low number of reads---------------------+")

cmat2 <- cts.batch1[, -c(2,4,7,8,17,18)]

dds2   <- DESeqDataSetFromMatrix(countData=cmat2, colData=samT.batch1[-c(2,4,7,8,17:18),], design=~Condition+Treatment)
dds2   <- DESeq(dds2, parallel=TRUE)
dds2   <- estimateSizeFactors(dds2)

vsd2 <- vst(dds2,     blind=F)
colData(vsd2)

## hclust
rv2 <- rowVars(assay(vsd2))
o   <- order(rv2,decreasing=TRUE)
dists2 <- dist(t(assay(vsd2)[head(o,2000),]))
hc2    <- hclust(dists2)
rm6hc  <- plot(hc2, labels=vsd2$sample)


## pca
pcaData2    <- plotPCA(vsd2, ntop=TOPNUM, intgroup=c("Group2"), returnData=TRUE)
pcaData2$samplename <- colData(vsd2)$sample
pcaData2$condition  <- colData(vsd2)$Condition
pcaData2$treatment  <- colData(vsd2)$Treatment

percentVar2 <- round(100 * attr(pcaData2, "percentVar"))

pdf(paste0("./Figures_Tables/PCA_Figures/", Project, "-NRP1_rm6_Fig.PCA.T2000.pdf"),width=10,height=10)
par(bg=NA)
plt.all2 <-  ggplot(pcaData2, aes(PC1, PC2, color=Group2)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(Group2), group=paste0(Group2), 
                        label=Group2), alpha=0.1, label.fontsize =10, label.buffer = unit(1, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  xlab(paste0("PC1: ",percentVar2[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar2[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("PCA (Top ", TOPNUM, " Most Variable Genes)"),
       caption = "CAD_mt709_0001: Minoru Takaoka" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1) 
plt_rm6_cWD <-  ggplot(pcaData2[pcaData2$condition=="cWD",], aes(PC1, PC2, color=treatment, label=samplename)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(treatment), group=paste0(treatment), 
                        label=treatment), alpha=0.1, label.fontsize =10, label.buffer = unit(1, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  geom_text_repel(aes(label=samplename), show.legend = FALSE, size=2)+
  xlab(paste0("PC1: ",percentVar2[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar2[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("Rm6 cWD: PCA "),
       caption = "" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1) 

plt_rm6_iWD <-  ggplot(pcaData2[pcaData2$condition=="iWD",], aes(PC1, PC2, color=treatment, label=samplename)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(treatment), group=paste0(treatment), 
                        label=treatment), alpha=0.1, label.fontsize =10, label.buffer = unit(1, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  geom_text_repel(aes(label=samplename), show.legend = FALSE, size=2)+
  xlab(paste0("PC1: ",percentVar2[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar2[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("Rm6 iWD: PCA"),
       caption = "" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1)  


bottom_row <- plot_grid(plt_rm6_cWD, plt_rm6_iWD, labels = c('B', 'C'), align = 'h', rel_widths = c(1, 1))

plot_grid(plt.all2, bottom_row, labels = c('A', ''), ncol = 1, align = 'vh', rel_heights = c(1, 1))
dev.off()

message("+-------------------------------------------------------------------------------+")
message("+---       Pairwise DESeq Analysis, KO vs WT                               -----+")
message("+-------------------------------------------------------------------------------+")

samT.rm6       <- as.data.frame(samT.batch1[-c(2,4,7,8,17:18),], stringsAsFactors = TRUE)
control.names  <- c("cWD_WT", "iWD_WT")
compare.names  <- c("cWD_KO", "iWD_KO")
selGroups      <- c("cWD", "iWD")

for(i in 1:2){
  print(i)
  samT.sub     <- subset(samT.rm6, Condition==selGroups[i])
  samT.sub$Treatment  <- relevel(samT.sub$Treatment, "WT")
  mat.sub      <- cmat2[, colnames(cmat2)%in%samT.sub$sample]
  colnames(mat.sub)  <- as.factor(colnames(mat.sub))
  dds.sub      <- DESeqDataSetFromMatrix(countData = mat.sub,
                                         colData = samT.sub,
                                         design= ~ Treatment)
  dds.sub      <- estimateSizeFactors(dds.sub)
  sizeFactors(dds.sub)
  
  dds.sub      <- DESeq(dds.sub, parallel=TRUE)
  vst.sub      <- vst(dds.sub,     blind=F)
  resultsNames(dds.sub)
  colData(vst.sub)
  save(mat.sub, samT.sub, dds.sub, vst.sub, 
       file = paste0("./Data/SLX-22500/",selGroups[i], "_DESeq2_dds_vstN.RData"))
  gc()
}

for(i in 1:2){
  load(paste0("./Data/SLX-22500/",selGroups[i], "_DESeq2_dds_vstN.RData"))
  pdf(paste0("./Figures_tables/PCA_Figures/DESeq2_pairwisePCA_", compare.names[i], "_", control.names[i], "_N.pdf"), width=8, height= 6)
  pca <- customPCA(samT.sub, assay(vst.sub), TOPNUM, ensEMBL2id) 
  plot_grid(pca[[1]], pca[[2]], pca[[3]], pca[[4]], nrow=2, byrow=T)
  dev.off()
}


DES.files <- paste0("./Data/SLX-22500/Rm6_",selGroups, "_DESeq2_dds_vstN.RData")
RES.files  <- paste0("./Data/SLX-22500/Rm6_",selGroups,  "_res_all_merE_summaryN.csv")
RES.files1 <- paste0("./Data/SLX-22500/Rm6_",selGroups,  "_res_all_lfcshrink_merE_summaryN.csv")
RESsig.files  <- paste0("./Data/SLX-22500/Rm6_",selGroups,  "_res_Sig_merE_summaryN.csv")
RESsig.files1 <- paste0("./Data/SLX-22500/Rm6_",selGroups,  "_res_Sig_lfcshrink_merE_summaryN.csv")


for(i in 1:2){
  load(DES.files[i])
  res       <- results(dds.sub, contrast = c("Treatment", "KO", "WT"))
  res       <- subset(res, !is.na(padj))
  res.sig   <- subset(res, padj <= significance & abs(log2FoldChange) >= l2fc)
  print(nrow(res.sig)); print(table(sign(res.sig$log2FoldChange)))
  write.csv(res.sig, file = RESsig.files[i])
  res.merEC <- merECRes_function(res, ensEMBL2id, dds.sub,RES.files[i])
  res1      <- lfcShrink(dds.sub, coef=2)
  res1.sig  <- subset(res1, padj <= significance & abs(log2FoldChange) >= l2fc)
  print(nrow(res1.sig)); print(table(sign(res1.sig$log2FoldChange)))
  write.csv(res1.sig, file = RESsig.files1[i])
  res.merEC1 <- merECRes_function(res1, ensEMBL2id, dds.sub, RES.files1[i])
}

## cWD, 80(60-, 20+) res; lfcshrink 57 (-1, 43; 1 14);
## iWD, 108(58-, 50+) res; lfcshrink 75 (-1, 32; 1 43);

message("+------------Overlap sigDEGs checking ---------------------+")
## Rm4
cWDrm4.res <- read.csv(paste0("./Data/SLX-22500/Rm4_",selGroups[1],  "_res_all_merE_summaryN.csv"), header=T)
iWDrm4.res <- read.csv(paste0("./Data/SLX-22500/Rm4_",selGroups[2],  "_res_all_merE_summaryN.csv"), header=T)

cWDrm6.res <- read.csv(paste0("./Data/SLX-22500/Rm6_",selGroups[1],  "_res_all_merE_summaryN.csv"), header=T)
iWDrm6.res <- read.csv(paste0("./Data/SLX-22500/Rm6_",selGroups[2],  "_res_all_merE_summaryN.csv"), header=T)

## select significant ones
cWDrm4.sig <- subset(cWDrm4.res, padj <= 0.05 & abs(log2FoldChange)>=1)
iWDrm4.sig <- subset(iWDrm4.res, padj <= 0.05 & abs(log2FoldChange)>=1)

cWDrm6.sig <- subset(cWDrm6.res, padj <= 0.05 & abs(log2FoldChange)>=1) ## the same as rm4
iWDrm6.sig <- subset(iWDrm6.res, padj <= 0.05 & abs(log2FoldChange)>=1)

cWDrm4.sig <- cWDrm4.sig[order(-cWDrm4.sig$log2FoldChange),]
iWDrm4.sig <- iWDrm4.sig[order(-iWDrm4.sig$log2FoldChange),]
iWDrm6.sig <- iWDrm6.sig[order(-iWDrm6.sig$log2FoldChange),]

## 
write.xlsx(cWDrm4.sig, file = "./Data/SLX-22500/Significant_DEGs_rm4_rm6_07_Feb.xlsx", append=F,
           row.names= F, sheetName = "Rm4_cWD_KOvsWT")
write.xlsx(iWDrm4.sig, file = "./Data/SLX-22500/Significant_DEGs_rm4_rm6_07_Feb.xlsx", append=T,
           row.names= F, sheetName = "Rm4_iWD_KOvsWT")
write.xlsx(iWDrm6.sig, file = "./Data/SLX-22500/Significant_DEGs_rm4_rm6_07_Feb.xlsx", append=T,
           row.names= F, sheetName = "Rm6_iWD_KOvsWT")

cWD_iWD_commonRm4  <- merge(cWDrm4.sig[,2:9], iWDrm4.sig[,1:8], by = "external_gene_name")
cWD_iWD_common1 <- cWD_iWD_commonrm4[sign(cWD_iWD_common$log2FoldChange.x)==sign(cWD_iWD_common$log2FoldChange.y),]

colnames(cWD_iWD_common1) <- c("GeneSymbol", "baseMean_cWD", "l2fc_cWD", "lfcSE_cWD", "stat_cWD", "pvalue_cWD",
                               "padj_cWD", "description", "GeneID", "baseMean_iWD", "l2fc_iWD", "lfcSE_iWD", 
                               "stat_iWD", "pvalue_iWD", "padj_iWD")
write.xlsx(cWD_iWD_common1, file = "./Data/SLX-22500/Significant_DEGs_rm4_rm6_07_Feb.xlsx", append=T,
           row.names= F, sheetName = "Rm4_cWD_iWD_commonGenes")

cWD_iWD_commonRm6  <- merge(cWDrm6.sig[,2:9], iWDrm6.sig[,1:8], by = "external_gene_name")
cWD_iWD_common2    <- cWD_iWD_commonRm6[sign(cwD_iWD_commonRm6$log2FoldChange.x)==sign(cwD_iWD_commonRm6$log2FoldChange.y),]

colnames(cWD_iWD_common2) <- c("GeneSymbol", "baseMean_cWD", "l2fc_cWD", "lfcSE_cWD", "stat_cWD", "pvalue_cWD",
                               "padj_cWD", "description", "GeneID", "baseMean_iWD", "l2fc_iWD", "lfcSE_iWD", 
                               "stat_iWD", "pvalue_iWD", "padj_iWD")
write.xlsx(cWD_iWD_common2, file = "./Data/SLX-22500/Significant_DEGs_rm4_rm6_07_Feb.xlsx", append=T,
           row.names= F, sheetName = "Rm6_cWD_iWD_commonGenes")

## unique genes
cWD_unique_rm4 <- cWDrm4.sig[cWDrm4.sig$external_gene_name%in%cWD_iWD_common1$GeneSymbol==F,]
iWD_unique_rm4 <- iWDrm4.sig[iWDrm4.sig$external_gene_name%in%cWD_iWD_common1$GeneSymbol==F,]
iWD_unique_rm6 <- iWDrm6.sig[iWDrm6.sig$external_gene_name%in%cWD_iWD_common2$GeneSymbol==F,]
write.xlsx(cWD_unique_rm4, file = "./Data/SLX-22500/Significant_DEGs_rm4_rm6_07_Feb.xlsx", append=T,
           row.names= F, sheetName = "Rm4_cWD_unique_Genes")
write.xlsx(iWD_unique_rm4, file = "./Data/SLX-22500/Significant_DEGs_rm4_rm6_07_Feb.xlsx", append=T,
           row.names= F, sheetName = "Rm4_iWD_unique_Genes")
write.xlsx(iWD_unique_rm6, file = "./Data/SLX-22500/Significant_DEGs_rm4_rm6_07_Feb.xlsx", append=T,
           row.names= F, sheetName = "Rm6_iWD_unique_Genes")

## function
functionPlotDECorrelation <- function(clustA, clustB, clustA.lab, clustB.lab, sig_cut, logfc_cut, topN, outdir) {
  
  rownames(clustA)  <- clustA$ensembl_gene_id
  rownames(clustB)  <- clustB$ensembl_gene_id
  
  clustA <- clustA[,c("log2FoldChange",  "padj",  "external_gene_name")]
  clustB <- clustB[,c("log2FoldChange",  "padj",  "external_gene_name")]
  
  colnames(clustA)  <- c("log2FoldChange.A",  "padj.A",  "external_gene_name.A")
  colnames(clustB)  <- c("log2FoldChange.B",  "padj.B",  "external_gene_name.B")
  
  clustA.topN <- subset(clustA, clustA$log2FoldChange.A > 0      & clustA$padj.A<=sig_cut)[1:topN,] 
  clustA.botN <- tail(subset(clustA, clustA$log2FoldChange.A < 0 & clustA$padj.A<=sig_cut),topN)
  
  clustB.topN <- subset(clustB, clustB$log2FoldChange.B > 0      & clustB$padj.B<=sig_cut)[1:topN,] 
  clustB.botN <- tail(subset(clustB, clustB$log2FoldChange.B < 0 & clustB$padj.B<=sig_cut),topN)
  
  compare.A.B <- merge(clustA, clustB, by='row.names', all=TRUE)
  compare.A.B[ is.na(compare.A.B)] <- 0
  compare.A.B$colour[(compare.A.B$padj.A <= sig_cut & compare.A.B$padj.B <= sig_cut)] <- "purple"# & abs(compare.A.B$log2FoldChange.A) >= logfc_cut & abs(compare.A.B$log2FoldChange.B) >= logfc_cut)] <- "purple"
  
  compare.A.B$colour[(compare.A.B$padj.A > sig_cut  & compare.A.B$padj.B <= sig_cut)] <- "blue"# & abs(compare.A.B$log2FoldChange.A) < logfc_cut & abs(compare.A.B$log2FoldChange.B) >= logfc_cut)] <- "blue"
  compare.A.B$colour[(compare.A.B$padj.A == 0       & compare.A.B$padj.B <= sig_cut)] <- "blue"# & abs(compare.A.B$log2FoldChange.A) < logfc_cut & abs(compare.A.B$log2FoldChange.B) >= logfc_cut)] <- "blue"
  
  compare.A.B$colour[(compare.A.B$padj.A <= sig_cut & compare.A.B$padj.B > sig_cut)] <- "darkgreen"# & abs(compare.A.B$log2FoldChange.A) >= logfc_cut & abs(compare.A.B$log2FoldChange.B) < logfc_cut)] <- "darkgreen"
  compare.A.B$colour[(compare.A.B$padj.A <= sig_cut & compare.A.B$padj.B == 0      )] <- "darkgreen"#& abs(compare.A.B$log2FoldChange.A) >= logfc_cut & abs(compare.A.B$log2FoldChange.B) < logfc_cut)] <- "darkgreen"
  
  compare.A.B$colour[(compare.A.B$padj.A > sig_cut & compare.A.B$padj.B > sig_cut)] <- "grey"
  compare.A.B$colour[(compare.A.B$padj.A > sig_cut & compare.A.B$padj.B == 0 )]     <- "grey"
  compare.A.B$colour[(compare.A.B$padj.A == 0      & compare.A.B$padj.B > sig_cut)] <- "grey"
  
  
  compare.A.B$colour[( abs(compare.A.B$log2FoldChange.A) < logfc_cut & abs(compare.A.B$log2FoldChange.B) < logfc_cut)] <- "grey"
  
  clustB.blue      <- subset(compare.A.B, compare.A.B$colour == "blue" & compare.A.B$padj.B < sig_cut)
  clustB.blue.topN <- clustB.blue[order(clustB.blue$log2FoldChange.B,decreasing=TRUE),] 
  clustB.blue.topN <- clustB.blue.topN[c(1:5),]
  clustB.blue.botN <- clustB.blue[order(clustB.blue$log2FoldChange.B,decreasing=FALSE),] 
  clustB.blue.botN <- clustB.blue.botN[c(1:5),]
  
  write.csv(clustB.blue, file=paste0(outdir, "/", Project, "_Table.DEvsDE.cWD-cWD_KO_vs_WT_VS.iWD_KO_vs_WT", "_Blue", ".ann.csv"))
  
  
  clustA.green      <- subset(compare.A.B, compare.A.B$colour == "darkgreen" & compare.A.B$padj.A < sig_cut)
  clustA.green.topN <- clustA.green[order(clustA.green$log2FoldChange.A,decreasing=TRUE),] 
  clustA.green.topN <- clustA.green.topN[c(1:5),]
  clustA.green.botN <- clustA.green[order(clustA.green$log2FoldChange.A,decreasing=FALSE),] 
  clustA.green.botN <- clustA.green.botN[c(1:5),]
  
  write.csv(clustA.green, file=paste0(outdir, "/", Project, "_Table.DEvsDE.iWD-cWD_KO_vs_WT_VS.iWD_KO_vs_WT", "_Green", ".ann.csv"))
  
  
  compare.A.B$label <- 0
  compare.A.B[compare.A.B$Row.names %in% unlist(rownames(clustA.topN)), ]$label  <- 1
  compare.A.B[compare.A.B$Row.names %in% unlist(rownames(clustA.botN)), ]$label  <- 1
  compare.A.B[compare.A.B$Row.names %in% unlist(rownames(clustB.topN)), ]$label  <- 1
  compare.A.B[compare.A.B$Row.names %in% unlist(rownames(clustB.botN)), ]$label  <- 1
  
  compare.A.B[compare.A.B$Row.names %in% unlist(clustB.blue.topN$Row.names),  ]$label  <- 1
  compare.A.B[compare.A.B$Row.names %in% unlist(clustB.blue.botN$Row.names),  ]$label  <- 1
  compare.A.B[compare.A.B$Row.names %in% unlist(clustA.green.topN$Row.names), ]$label  <- 1
  compare.A.B[compare.A.B$Row.names %in% unlist(clustA.green.botN$Row.names), ]$label  <- 1
  
  minFC <- -8 #min( compare.A.B.wt.ko$avg_logFC.A, compare.A.B.wt.ko$avg_logFC.B  )
  maxFC <-  4 #max( compare.A.B.wt.ko$avg_logFC.A, compare.A.B.wt.ko$avg_logFC.B  )
  
  cor.plt <- ggplot(data=compare.A.B, aes(x=log2FoldChange.A, y=log2FoldChange.B, colour=colour, label=external_gene_name.B)) +
    geom_vline(xintercept = -(logfc_cut), linetype="dashed", colour="black", size=0.25, alpha=0.5) +
    geom_vline(xintercept = logfc_cut,    linetype="dashed", colour="black", size=0.25, alpha=0.5) +
    geom_hline(yintercept = -(logfc_cut), linetype="dashed", colour="black", size=0.25, alpha=0.5) +
    geom_hline(yintercept = logfc_cut,    linetype="dashed", colour="black", size=0.25, alpha=0.5) +
    geom_abline(intercept = 0, linetype="solid", colour="black", alpha=0.5) +
    geom_point(alpha=0.5, size=1) +  
    geom_label_repel( data=subset(compare.A.B, compare.A.B$label == 1 ), force=10, fill='white',
                      show.legend = FALSE, nudge_x=0.0, nudge_y=0.0, segment.size = 0.25, size=3) +
    coord_fixed() +
    scale_x_continuous(limits=c(minFC, maxFC), breaks=seq(minFC, maxFC,2)) + 
    scale_y_continuous(limits=c(minFC, maxFC), breaks=seq(minFC, maxFC,2)) +
    xlab(paste0("log2FC (", clustA.lab, " C Vs H) ")) + 
    ylab(paste0("log2FC (", clustB.lab, " C Vs H) ")) +
    scale_colour_manual(name="", values=c("purple"="purple", "blue"="blue", "darkgreen"="darkgreen", "grey"="grey"), 
                        labels=c("purple"=paste0(clustA.lab, " & ", clustB.lab), 
                                 "blue"=paste0(clustB.lab, " only"), "darkgreen"=paste0(clustA.lab, " only"))) +
    ggtitle(paste0(clustA.lab, " (KO/WT) Vs ", clustB.lab, " (KO/WT) [log2FC=", logfc_cut,"]")) +
    theme(text=element_text(family="sans"), legend.position="bottom", aspect.ratio=1 )
  
  return(cor.plt)
  
}


#test <- functionPlotDECorrelation(res.group.BP_C_vs_BP_H.ann,res.group.EPL_C_vs_EPL_H.ann, "BP", "EPL", significance, 1.5, topN)
#test
outdir <- paste0(basedir, "/Data/SLX-22500/")
test <- functionPlotDECorrelation(cWDrm4.res,iWDrm4.res, "cWD", "iWD", 0.05, 1, 20, outdir)

message("+------------Hcluster plot for all checking ---------------------+")

pdf("./Data/SLX-22500/Hcluster_all_rm2_rm4_rm6_plot_8_02_2023.pdf", width = 8, height=20)
par(mfrow=c(4,1))
plot(hc, labels=vsd$sample)
plot(hc0, labels=vsd0$sample)
plot(hc1, labels=vsd1$sample)
plot(hc2, labels=vsd2$sample)
dev.off()

message("+------------GeneOntology analysis ---------------------+")

load("./Data/References_Data/GRCm39/Ensembl_GRCm39_entrezID_ensembID_exterName.RData")

cWDrm4.resm   <- merge(cWDrm4.res[,c(1,3,7,8)], ensEMBL2id.GO[,c(1,3)], by = "ensembl_gene_id")
iWDrm4.resm   <- merge(iWDrm4.res[,c(1,3,7,8)], ensEMBL2id.GO[,c(1,3)], by = "ensembl_gene_id")
iWDrm6.resm   <- merge(iWDrm6.res[,c(1,3,7,8)], ensEMBL2id.GO[,c(1,3)], by = "ensembl_gene_id")

resdfbk.cWD     <- subset(cWDrm4.resm, padj>0.05 & abs(log2FoldChange)<1)
resdfbk.iWDrm4  <- subset(iWDrm4.resm, padj>0.05 & abs(log2FoldChange)<1)
resdfbk.iWDrm6  <- subset(iWDrm6.resm, padj>0.05 & abs(log2FoldChange)<1)

resdf.cWD     <- subset(cWDrm4.resm, padj<=0.05 & abs(log2FoldChange)>=1)
resdf.iWDrm4  <- subset(iWDrm4.resm, padj<=0.05 & abs(log2FoldChange)>=1)
resdf.iWDrm6  <- subset(iWDrm6.resm, padj<=0.05 & abs(log2FoldChange)>=1)    

ColNames <- c("ENTREZID", "ENSEMBL", "SYMBOL", "L2FC") 
resdfbk.cWD1 <- resdfbk.cWD[,c(5,1,4,2)] ; colnames(resdfbk.cWD1) <- ColNames; resdfbk.cWD1 <- resdfbk.cWD1[order(resdfbk.cWD1$L2FC),]
resdfbk.iWD1 <- resdfbk.iWDrm4[,c(5,1,4,2)] ; colnames(resdfbk.iWD1) <- ColNames; resdfbk.iWD1 <- resdfbk.iWD1[order(resdfbk.iWD1$L2FC),]
resdfbk.iWD2 <- resdfbk.iWDrm6[,c(5,1,4,2)] ; colnames(resdfbk.iWD2) <- ColNames; resdfbk.iWD2 <- resdfbk.iWD2[order(resdfbk.iWD2$L2FC),]

resdf.cWD1 <- resdf.cWD[,c(5,1,4,2)] ; colnames(resdf.cWD1) <- ColNames; resdf.cWD1 <- resdf.cWD1[order(resdf.cWD1$L2FC),]
resdf.iWD1 <- resdf.iWDrm4[,c(5,1,4,2)] ; colnames(resdf.iWD1) <- ColNames; resdf.iWD1 <- resdf.iWD1[order(resdf.iWD1$L2FC),]
resdf.iWD2 <- resdf.iWDrm6[,c(5,1,4,2)] ; colnames(resdf.iWD2) <- ColNames; resdf.iWD2 <- resdf.iWD2[order(resdf.iWD2$L2FC),]


geneList.cWD  <- resdf.cWD1$L2FC; names(geneList.cWD) <- resdf.cWD1$SYMBOL
geneList.iWD1 <- resdf.iWD1$L2FC; names(geneList.iWD1) <- resdf.iWD1$SYMBOL
geneList.iWD2 <- resdf.iWD2$L2FC; names(geneList.iWD2) <- resdf.iWD2$SYMBOL


geneListZ.cWD  <- resdf.cWD1$L2FC; names(geneList.cWD) <- resdf.cWD1$ENTREZID
geneListZ.iWD1 <- resdf.iWD1$L2FC; names(geneList.iWD1) <- resdf.iWD1$ENTREZID
geneListZ.iWD2 <- resdf.iWD2$L2FC; names(geneList.iWD2) <- resdf.iWD2$ENTREZID

## enrichGO analysis

GO.cWD1 <- enrichGO(gene = resdf.cWD1$SYMBOL, OrgDb = org.Mm.eg.db,
                    #universe = resdfbk.cWD1$SYMBOL,
                    keyType = "SYMBOL", ont = "ALL",
                    pAdjustMethod = "BH", pvalueCutoff  = 0.05 )
GO.iWD1 <- enrichGO(gene = resdf.iWD1$SYMBOL, OrgDb = org.Mm.eg.db,
                    #universe = resdfbk.iWD1$SYMBOL,
                    keyType = "SYMBOL", ont = "ALL",
                    pAdjustMethod = "BH", pvalueCutoff  = 0.05 )
GO.iWD2 <- enrichGO(gene = resdf.iWD2$SYMBOL, OrgDb = org.Mm.eg.db,
                    #universe = resdfbk.cWD1$SYMBOL,
                    keyType = "SYMBOL", ont = "ALL",
                    pAdjustMethod = "BH", pvalueCutoff  = 0.05 )

GO.KEGG.iWD2 <-  enrichKEGG(resdf.iWD2$ENTREZID, organism = "mmu")
GO.KEGG.iWD2 <- setReadable(GO.KEGG.iWD2, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")

save(GO.cWD1, GO.iWD1, GO.iWD2,GO.KEGG.iWD2, file = "./Data/SLX-22500/GeneOntology_Analysis_cWDrm4_iWDrm46_08_02_2023.RData")

write.xlsx(as.data.frame(GO.cWD1), file = "./Data/SLX-22500/GeneOntology_Analysis_cWDrm4_iWDrm46_08_02_2023.xlsx",
           append = F, row.names=F, sheetName = "cWD_KOvsWT_GO")
write.xlsx(as.data.frame(GO.iWD1), file = "./Data/SLX-22500/GeneOntology_Analysis_cWDrm4_iWDrm46_08_02_2023.xlsx",
           append = T, row.names=F, sheetName = "iWD_rm4_KOvsWT_GO")
write.xlsx(as.data.frame(GO.iWD2), file = "./Data/SLX-22500/GeneOntology_Analysis_cWDrm4_iWDrm46_08_02_2023.xlsx",
           append = T, row.names=F, sheetName = "iWD_rm6_KOvsWT_GO")
write.xlsx(as.data.frame(GO.KEGG.iWD2), file = "./Data/SLX-22500/GeneOntology_Analysis_cWDrm4_iWDrm46_08_02_2023.xlsx",
           append = T, row.names=F, sheetName = "iWD_rm6_KOvsWT_KEGG")


message("+--- After discussing with Ziad on 9th Feb, 2023, the sig foldchange level will be 0.6 instead of 1----+")

## select significant ones
cWDrm4.sig <- subset(cWDrm4.res, padj <= 0.05 & abs(log2FoldChange)>=0.6)
iWDrm4.sig <- subset(iWDrm4.res, padj <= 0.05 & abs(log2FoldChange)>=0.6)

cWDrm6.sig <- subset(cWDrm6.res, padj <= 0.05 & abs(log2FoldChange)>=0.6) ## the same as rm4
iWDrm6.sig <- subset(iWDrm6.res, padj <= 0.05 & abs(log2FoldChange)>=0.6)

cWDrm4.sig <- cWDrm4.sig[order(-cWDrm4.sig$log2FoldChange),]
iWDrm4.sig <- iWDrm4.sig[order(-iWDrm4.sig$log2FoldChange),]
iWDrm6.sig <- iWDrm6.sig[order(-iWDrm6.sig$log2FoldChange),]

## 
write.xlsx(cWDrm4.sig, file = "./Data/SLX-22500/Significant_DEGs_rm4_rm6_l2fc06_09_Feb.xlsx", append=F,
           row.names= F, sheetName = "Rm4_cWD_KOvsWT")
write.xlsx(iWDrm4.sig, file = "./Data/SLX-22500/Significant_DEGs_rm4_rm6_l2fc06_09_Feb.xlsx", append=T,
           row.names= F, sheetName = "Rm4_iWD_KOvsWT")
write.xlsx(iWDrm6.sig, file = "./Data/SLX-22500/Significant_DEGs_rm4_rm6_l2fc06_09_Feb.xlsx", append=T,
           row.names= F, sheetName = "Rm6_iWD_KOvsWT")

cWD_iWD_commonRm4  <- merge(cWDrm4.sig[,2:9], iWDrm4.sig[,1:8], by = "external_gene_name")
cWD_iWD_common1 <- cWD_iWD_commonRm4[sign(cWD_iWD_commonRm4$log2FoldChange.x)==sign(cWD_iWD_commonRm4$log2FoldChange.y),]

colnames(cWD_iWD_common1) <- c("GeneSymbol", "baseMean_cWD", "l2fc_cWD", "lfcSE_cWD", "stat_cWD", "pvalue_cWD",
                               "padj_cWD", "description", "GeneID", "baseMean_iWD", "l2fc_iWD", "lfcSE_iWD", 
                               "stat_iWD", "pvalue_iWD", "padj_iWD")
write.xlsx(cWD_iWD_common1, file = "./Data/SLX-22500/Significant_DEGs_rm4_rm6_l2fc06_09_Feb.xlsx", append=T,
           row.names= F, sheetName = "Rm4_cWD_iWD_commonGenes")

cWD_iWD_commonRm6  <- merge(cWDrm6.sig[,2:9], iWDrm6.sig[,1:8], by = "external_gene_name")
cWD_iWD_common2    <- cWD_iWD_commonRm6[sign(cwD_iWD_commonRm6$log2FoldChange.x)==sign(cwD_iWD_commonRm6$log2FoldChange.y),]

colnames(cWD_iWD_common2) <- c("GeneSymbol", "baseMean_cWD", "l2fc_cWD", "lfcSE_cWD", "stat_cWD", "pvalue_cWD",
                               "padj_cWD", "description", "GeneID", "baseMean_iWD", "l2fc_iWD", "lfcSE_iWD", 
                               "stat_iWD", "pvalue_iWD", "padj_iWD")
write.xlsx(cWD_iWD_common2, file = "./Data/SLX-22500/Significant_DEGs_rm4_rm6_l2fc06_09_Feb.xlsx", append=T,
           row.names= F, sheetName = "Rm6_cWD_iWD_commonGenes")

## unique genes
cWD_unique_rm4 <- cWDrm4.sig[cWDrm4.sig$external_gene_name%in%cWD_iWD_common1$GeneSymbol==F,]
iWD_unique_rm4 <- iWDrm4.sig[iWDrm4.sig$external_gene_name%in%cWD_iWD_common1$GeneSymbol==F,]
iWD_unique_rm6 <- iWDrm6.sig[iWDrm6.sig$external_gene_name%in%cWD_iWD_common2$GeneSymbol==F,]
write.xlsx(cWD_unique_rm4, file = "./Data/SLX-22500/Significant_DEGs_rm4_rm6_l2fc06_09_Feb.xlsx", append=T,
           row.names= F, sheetName = "Rm4_cWD_unique_Genes")
write.xlsx(iWD_unique_rm4, file = "./Data/SLX-22500/Significant_DEGs_rm4_rm6_l2fc06_09_Feb.xlsx", append=T,
           row.names= F, sheetName = "Rm4_iWD_unique_Genes")
write.xlsx(iWD_unique_rm6, file = "./Data/SLX-22500/Significant_DEGs_rm4_rm6_l2fc06_09_Feb.xlsx", append=T,
           row.names= F, sheetName = "Rm6_iWD_unique_Genes")


message("+---GeneOntology Analysis-------------------------------+")
cWDrm4.resm   <- merge(cWDrm4.res[,c(1,3,7,8)], ensEMBL2id.GO[,c(1,3)], by = "ensembl_gene_id")
iWDrm4.resm   <- merge(iWDrm4.res[,c(1,3,7,8)], ensEMBL2id.GO[,c(1,3)], by = "ensembl_gene_id")
iWDrm6.resm   <- merge(iWDrm6.res[,c(1,3,7,8)], ensEMBL2id.GO[,c(1,3)], by = "ensembl_gene_id")

resdfbk.cWD     <- subset(cWDrm4.resm, padj>0.05 & abs(log2FoldChange)<0.6)
resdfbk.iWDrm4  <- subset(iWDrm4.resm, padj>0.05 & abs(log2FoldChange)<0.6)
resdfbk.iWDrm6  <- subset(iWDrm6.resm, padj>0.05 & abs(log2FoldChange)<0.6)

resdf.cWD     <- subset(cWDrm4.resm, padj<=0.05 & abs(log2FoldChange)>=0.6)
resdf.iWDrm4  <- subset(iWDrm4.resm, padj<=0.05 & abs(log2FoldChange)>=0.6)
resdf.iWDrm6  <- subset(iWDrm6.resm, padj<=0.05 & abs(log2FoldChange)>=0.6)    

ColNames <- c("ENTREZID", "ENSEMBL", "SYMBOL", "L2FC") 
resdfbk.cWD1 <- resdfbk.cWD[,c(5,1,4,2)] ; colnames(resdfbk.cWD1) <- ColNames; resdfbk.cWD1 <- resdfbk.cWD1[order(resdfbk.cWD1$L2FC),]
resdfbk.iWD1 <- resdfbk.iWDrm4[,c(5,1,4,2)] ; colnames(resdfbk.iWD1) <- ColNames; resdfbk.iWD1 <- resdfbk.iWD1[order(resdfbk.iWD1$L2FC),]
resdfbk.iWD2 <- resdfbk.iWDrm6[,c(5,1,4,2)] ; colnames(resdfbk.iWD2) <- ColNames; resdfbk.iWD2 <- resdfbk.iWD2[order(resdfbk.iWD2$L2FC),]

resdf.cWD1 <- resdf.cWD[,c(5,1,4,2)] ; colnames(resdf.cWD1) <- ColNames; resdf.cWD1 <- resdf.cWD1[order(resdf.cWD1$L2FC),]
resdf.iWD1 <- resdf.iWDrm4[,c(5,1,4,2)] ; colnames(resdf.iWD1) <- ColNames; resdf.iWD1 <- resdf.iWD1[order(resdf.iWD1$L2FC),]
resdf.iWD2 <- resdf.iWDrm6[,c(5,1,4,2)] ; colnames(resdf.iWD2) <- ColNames; resdf.iWD2 <- resdf.iWD2[order(resdf.iWD2$L2FC),]


geneList.cWD  <- resdf.cWD1$L2FC; names(geneList.cWD) <- resdf.cWD1$SYMBOL
geneList.iWD1 <- resdf.iWD1$L2FC; names(geneList.iWD1) <- resdf.iWD1$SYMBOL
geneList.iWD2 <- resdf.iWD2$L2FC; names(geneList.iWD2) <- resdf.iWD2$SYMBOL


geneListZ.cWD  <- resdf.cWD1$L2FC; names(geneList.cWD) <- resdf.cWD1$ENTREZID
geneListZ.iWD1 <- resdf.iWD1$L2FC; names(geneList.iWD1) <- resdf.iWD1$ENTREZID
geneListZ.iWD2 <- resdf.iWD2$L2FC; names(geneList.iWD2) <- resdf.iWD2$ENTREZID

## enrichGO analysis

GO.cWD1 <- enrichGO(gene = resdf.cWD1$SYMBOL, OrgDb = org.Mm.eg.db,
                    #universe = resdfbk.cWD1$SYMBOL,
                    keyType = "SYMBOL", ont = "ALL",
                    pAdjustMethod = "BH", pvalueCutoff  = 0.05 )
GO.iWD1 <- enrichGO(gene = resdf.iWD1$SYMBOL, OrgDb = org.Mm.eg.db,
                    #universe = resdfbk.iWD1$SYMBOL,
                    keyType = "SYMBOL", ont = "ALL",
                    pAdjustMethod = "BH", pvalueCutoff  = 0.05 )
GO.iWD2 <- enrichGO(gene = resdf.iWD2$SYMBOL, OrgDb = org.Mm.eg.db,
                    #universe = resdfbk.cWD1$SYMBOL,
                    keyType = "SYMBOL", ont = "ALL",
                    pAdjustMethod = "BH", pvalueCutoff  = 0.05 )

GO.KEGG.iWD1 <-  enrichKEGG(resdf.iWD1$ENTREZID, organism = "mmu")
GO.KEGG.iWD1 <- setReadable(GO.KEGG.iWD1, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
GO.KEGG.cWD1 <-  enrichKEGG(resdf.cWD1$ENTREZID, organism = "mmu")
GO.KEGG.cWD1 <- setReadable(GO.KEGG.cWD1, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
save(GO.cWD1, GO.iWD1, GO.iWD2,GO.KEGG.cWD1,GO.KEGG.iWD1, file = "./Data/SLX-22500/GeneOntology_Analysis_cWDrm4_iWDrm46_l2fc06_09_02_2023.RData")

write.xlsx(as.data.frame(GO.cWD1), file = "./Data/SLX-22500/GeneOntology_Analysis_cWDrm4_iWDrm46_l2fc06_09_02_2023.xlsx",
           append = F, row.names=F, sheetName = "cWD_KOvsWT_GO")
write.xlsx(as.data.frame(GO.iWD1), file = "./Data/SLX-22500/GeneOntology_Analysis_cWDrm4_iWDrm46_l2fc06_09_02_2023.xlsx",
           append = T, row.names=F, sheetName = "iWD_rm4_KOvsWT_GO")
write.xlsx(as.data.frame(GO.iWD2), file = "./Data/SLX-22500/GeneOntology_Analysis_cWDrm4_iWDrm46_l2fc06_09_02_2023.xlsx",
           append = T, row.names=F, sheetName = "iWD_rm6_KOvsWT_GO")
write.xlsx(as.data.frame(GO.KEGG.iWD1), file = "./Data/SLX-22500/GeneOntology_Analysis_cWDrm4_iWDrm46_l2fc06_09_02_2023.xlsx",
           append = T, row.names=F, sheetName = "iWD_rm4_KOvsWT_KEGG")
write.xlsx(as.data.frame(GO.KEGG.cWD1), file = "./Data/SLX-22500/GeneOntology_Analysis_cWDrm4_iWDrm46_l2fc06_09_02_2023.xlsx",
           append = T, row.names=F, sheetName = "cWD_rm4_KOvsWT_KEGG")
message("+---- selected BP for cWD barplot------------------------------+")

selBPIDs <- c("GO:0006937", "Go:0045785", "GO:0051090", "GO:0032970", "GO:0046718",
            "GO:2001237", "GO:0030198", "GO:0051495", "GO:0018108", "GO:0007409",
            "GO:0018212", "GO:0070372", "GO:0060840")
selBPs   <- as.data.frame(GO.cWD1)[as.data.frame(GO.cWD1)$ID%in%selBPIDs, "Description"]
library(enrichplot)
pdf("/Users/xz289/Documents/Minoru/Project_mt709_0001/Data/SLX-22500/cWD_KOvsWT_Selcted_NRP1_highvslow_Barplot_N.pdf", width= 9, height = 7.5, onefile=T)
barplot(GO.cWD1, showCategory = selBPs) + 
  geom_bar(stat="identity", aes(alpha = -p.adjust), color="white")+
  scale_fill_gradient(low="darkred",high="tomato1")+
  ylab("") +
  xlab("Gene Counts") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=18, face= "bold"),
        axis.title.x = element_text(size=20, face= "bold"),
        axis.title.y = element_text(size=20, face= "bold"),
        axis.text.y = element_text(size=14, face="bold"),
        legend.title = element_text(size = 14, face="bold"), 
        legend.text = element_text(size = 14, face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = "white", colour = NA)) 
dev.off()
message("+----------FINISH 08/02/2023----------------------------------+")
