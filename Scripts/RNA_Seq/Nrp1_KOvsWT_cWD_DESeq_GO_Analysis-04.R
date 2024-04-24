#!/usr/local/bin/Rscript
# SLX-22500 (GTC280), NRP1 KO/WT with cWD diet
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
# https://github.com/CAD-BFX/Takaoka_Mallat
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

message("+-------------------------------------------------------------------+")
message("+  General settings and libraries calling for DESeq2 analysis       +")
message("+-------------------------------------------------------------------+")

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
  library("xlsx")
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
})


register(MulticoreParam(2))
basedir     <- "/Users/xz289/Documents/Minoru/Project_mt709_0001"
setwd(basedir)

Project         <- "SLX-22500"
Tissues         <- "Nrp1"
significance    <- 0.05
l2fc            <- 0.6 
elementTextSize <- 10
TOPNUM          <- 2000

##source("./Scripts/Utility_Function.R") ## check


message("+-----             Load the data and ensemble            -----------+")

load("./Data/References_Data/GRCm39/Ensembl_mmusculus_ID_Name_Des_Chr_GRCm39.RData")
load("./Data/SLX-22500/deseq2.dds.RData")

cts.batch    <- assay(dds)
cts.batch    <- cts.batch[,order(colnames(cts.batch))]
samT.batch   <- read.csv("./Data/SLX-22500/CAD_mt709_0001_NRP1-cWD_SampleTable.csv", header=T)
samT.batch$Treatment    <- as.factor(samT.batch$Treatment)
samT.batch$Treatment    <- relevel(samT.batch$Treatment, "WT")

samT.batch  <- samT.batch[order(samT.batch$sample),]
cmat        <- cts.batch

message("+--- Generate correlation heatmap for all samples   ----------------+")

rownames(samT.batch) <-samT.batch$sample 
dds   <- DESeqDataSetFromMatrix(countData=cmat, colData=samT.batch, 
                                design=~Treatment)
dds   <- DESeq(dds, parallel=TRUE)
dds   <- estimateSizeFactors(dds)

vsd <- vst(dds,     blind=F)
colData(vsd)

## customised pca plot to check the informal clustering with treatment


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

pdf("./Figures_tables/PCA_Figures/DESeq2-pairwisePCA_cWD_KOvsWT_TOP2000MV.pdf", width=8, height= 6)
pca <- customPCA(samT.batch, assay(vsd), TOPNUM=2000, ensEMBL2id) 
plot_grid(pca[[1]], pca[[2]], pca[[3]], pca[[4]], nrow=2, byrow=T)
dev.off()

message("+-------------------------------------------------------------------+")
message("+---       Pairwise DESeq Analysis, KO vs WT                   -----+")
message("+Remove PMC1,2,8,10 for the following DEGs analysis,low reads ------+")
message("+-------------------------------------------------------------------+")

rmIndex  <- c("PMC-1", "PMC-2", "PMC-8", "PMC-10")
samT.sub <- samT.batch[-which(samT.batch$Sample.name%in%rmIndex==T),]
cmat.sub <- cmat[,which(colnames(cmat)%in%samT.sub$sample==T)]
dds.sub  <- DESeqDataSetFromMatrix(countData=cmat.sub, colData=samT.sub, 
                                design=~Treatment)
dds.sub  <- DESeq(dds.sub, parallel=TRUE)
dds.sub  <- estimateSizeFactors(dds.sub)

vsd.sub <- vst(dds.sub,     blind=F)
colData(vsd.sub)

save(dds.sub, file = "./Data/SLX-22500/Rm4_cWD_DESeq2_dds_vstN.RData")


## Identified DEGs and merge with ensemble annotation


RES.files  <- "./Data/SLX-22500/Significant_DEGs_rm4_l2fc06_09_Feb.xlsx"

## merge ensemble function
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

res       <- results(dds.sub, contrast = c("Treatment", "KO", "WT"))
res       <- subset(res, !is.na(padj))
res.sig   <- subset(res, padj <= significance & abs(log2FoldChange) >= l2fc)
res.sig   <- res.sig[order(-res.sig$log2FoldChange),]
print(nrow(res.sig)); print(table(sign(res.sig$log2FoldChange)))
write.csv(res.sig, file = RESsig.file)
res.merEC <- merECRes_function(res, ensEMBL2id, dds.sub,RES.file)


message("+-------              GeneOntology analysis        -----------------+")

## load ensembl annotation with entrezID
load("./Data/References_Data/GRCm39/Ensembl_GRCm39_entrezID_ensembID_exterName.RData")

## prepare input data
cWDrm4.res    <- read.csv(RES.file, header = T)
cWDrm4.resm   <- merge(cWDrm4.res[,c(1,3,7,8)], ensEMBL2id.GO[,c(1,3)], by = "ensembl_gene_id")
resdfbk.cWD     <- subset(cWDrm4.resm, padj>0.05 & abs(log2FoldChange)<0.6)
resdf.cWD     <- subset(cWDrm4.resm, padj<=0.05 & abs(log2FoldChange)>=0.6)

ColNames <- c("ENTREZID", "ENSEMBL", "SYMBOL", "L2FC") 
resdfbk.cWD1 <- resdfbk.cWD[,c(5,1,4,2)] ; colnames(resdfbk.cWD1) <- ColNames; resdfbk.cWD1 <- resdfbk.cWD1[order(resdfbk.cWD1$L2FC),]
resdf.cWD1 <- resdf.cWD[,c(5,1,4,2)] ; colnames(resdf.cWD1) <- ColNames; resdf.cWD1 <- resdf.cWD1[order(resdf.cWD1$L2FC),]

geneList.cWD  <- resdf.cWD1$L2FC; names(geneList.cWD) <- resdf.cWD1$SYMBOL
geneListZ.cWD  <- resdf.cWD1$L2FC; names(geneList.cWD) <- resdf.cWD1$ENTREZID


## enrichGO analysis

GO.cWD1 <- enrichGO(gene = resdf.cWD1$SYMBOL, OrgDb = org.Mm.eg.db,
                    universe = resdfbk.cWD1$SYMBOL,
                    keyType = "SYMBOL", ont = "ALL",
                    pAdjustMethod = "BH", pvalueCutoff  = 0.05 )


GO.KEGG.cWD1 <-  enrichKEGG(resdf.cWD1$ENTREZID, organism = "mmu")
GO.KEGG.cWD1 <- setReadable(GO.KEGG.cWD1, OrgDb = org.Mm.eg.db, keyType = "ENTREZID")
save(GO.cWD1, GO.KEGG.cWD1, file = "./Data/SLX-22500/GeneOntology_Analysis_cWDrm4_l2fc06_09_02_2023.RData")

write.xlsx(as.data.frame(GO.cWD1), file = "./Data/SLX-22500/GeneOntology_Analysis_cWDrm4_l2fc06_09_02_2023.xlsx",
           append = F, row.names=F, sheetName = "cWD_KOvsWT_GO")
write.xlsx(as.data.frame(GO.KEGG.cWD1), file = "./Data/SLX-22500/GeneOntology_Analysis_cWDrm4_l2fc06_09_02_2023.xlsx",
           append = T, row.names=F, sheetName = "cWD_rm4_KOvsWT_KEGG")

message("+---Fig 4h: selected BP for cWD barplot NRP1 highvs low    ---------+")

selBPIDs <- c("GO:0006937", "GO:0045785",  "GO:0051090", "GO:0032970", "GO:2001237",
              "GO:0030198", "GO:0051495", "GO:0018108", "GO:0007409", "GO:0070372")
selBPs   <- as.data.frame(GO.cWD1)[as.data.frame(GO.cWD1)$ID%in%selBPIDs, "Description"]

selBPs.1 <- c("regulation of muscle contraction",
              "positive regulation of \ncell adhesion",
              "regulation of DNA-binding \ntranscription factor activity",
              "regulation of actin \nfilament-based process",
              "negative regulation of \nextrinsic apoptotic \nsignaling pathway", 
              "extracellular matrix organization", 
              "positive regulation of \ncytoskeleton organization",
              "peptidyl-tyrosine phosphorylation",
              "regulation of ERK1 and \nERK2 cascade",
              "axonogenesis")
MNRP1_GO   <- as.data.frame(GO.cWD1)
MBPplt_dat <- MNRP1_GO[MNRP1_GO$Description%in%selBPs,]
MBPplt_dat$Description <-selBPs.1
MBPplt_dat$Count <- as.numeric(unlist(lapply(MBPplt_dat$GeneRatio, 
                                             function(x) strsplit(x, split="/")[[1]][1])))



pdf("./Figures_Tables/Macrophage/Fig4h-cWD_KOvsWT_Selcted_NRP1_highvslow_Barplot_Apr_2024-rm3.pdf", 
    width= 13, height = 8, onefile=T)

MBPplt <- ggplot(MBPplt_dat, aes(x=reorder(Description, -p.adjust),
                                 y=Count, color=p.adjust, fill=p.adjust)) + 
  geom_bar(stat="identity", aes(alpha = -p.adjust), color="white")+
  coord_flip()+ 
  scale_fill_gradient(low="darkred",high="tomato1") +
  scale_y_continuous(limits = c(0,10)) +
  ylab("") +
  xlab("") +
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=18, face= "bold"),
        axis.title.x = element_text(size=20, face= "bold"),
        axis.title.y = element_text(size=20, face= "bold"),
        axis.text.y = element_text(size=20, face="bold"),
        legend.title = element_text(size = 20, face="bold"), 
        legend.text = element_text(size = 20, face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = "white", colour = NA)) +
  guides(fill=guide_legend("Padj"))
MBPplt
dev.off()


message("+----------FINISH 19/04/2023----------------------------------+")










