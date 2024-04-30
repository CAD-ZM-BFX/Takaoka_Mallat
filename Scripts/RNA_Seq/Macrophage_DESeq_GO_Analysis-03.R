#!/usr/local/bin/Rscript
# Macrophage-mt709_0001: RNASeq analysis for two Macrophage
# Mice deficient in the LDL receptor (LDLR KO mice) with high fat diet (HFD) feeding
# chow diet, intermittent(alternative) HFD, continuous HFD
# Macrophage: RNASeq for mouse with diet treatment, 2 groups intermittent HFD and continuous HFD.
# each has 3 replicates. Pool data, sample 1&2: 3 mice and sample 3: 4 mice
# Analysis Performed by Xiaohui Zhao
# School of Medicine, University of Cambridge, UK
# Copyright Xiaohui Zhao (xz289@cam.ac.uk)
#
# Script available from:
# https://github.com/CAD-BFX/Takaoka_Mallat
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
message("+General settings and libraries calling for Differential analysis   +")
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
  library("tidyverse")
  library("readxl")
  library("cividis")
  library('variancePartition')
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
  ## tissueEnrich
  library("TissueEnrich")
  library("UniProt.ws")
  library("corrplot") 
  library("clusterSim") 
 })


register(MulticoreParam(2))
basedir     <- "/Users/xz289/Documents/Minoru/Project_mt709_0001"
setwd(basedir)

Project         <- "CAD_mt709_0001"
Tissues         <- "Macrophage"
significance    <- 0.05
l2fc            <- 1 
elementTextSize <- 10
TOPNUM          <- 2000



message("+-----             Load the data and ensemble            -----------+")

load("./Data/References_Data/GRCm39/Ensembl_mmusculus_ID_Name_Des_Chr_GRCm39.RData")
load("./Data/Origin_Data/CAD_mt709_0001-RNASeq_Macrophage_BMonocyte_BMP_merge_Counts_SamT.RData")

cmat.mac  <- merge.cmatGN.new[rownames(merge.cmatGN.new)%in%ensEMBL2id$external_gene_name, 
                              grep("Mac", colnames(merge.cmatGN.new))]
## merge.cmat, samT.mer
rownames(samT.mer) <-samT.mer$SampleName 
samT.mac  <- samT.mer[grep("Mac", rownames(samT.mer)),]

pdf("./Figures_Tables/PCA_Figures/RNASeq_Macrophage_merge_Fig.PCA.T2000.pdf",width=10,height=10)
par(bg=NA)

plt.mac <-  ggplot(pcaData[1:12,], aes(PC1, PC2, color=condition, label=samplename)) +
  geom_mark_ellipse(aes(fill = NULL, color=paste0(batches), group=paste0(batches), 
                        label=batches), alpha=0.1, label.fontsize =10, label.buffer = unit(1, 'mm')) +
  geom_point(size=2, alpha=0.75) +
  geom_text_repel(aes(label=samplename), show.legend = FALSE, size=2)+
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
  coord_fixed() +
  theme_bw() +
  labs(title = paste0("Macrophage: PCA (Top ", TOPNUM, " Most Variable Genes)"),
       caption = "CAD_mt709_0001: Minoru Takaoka" ) +
  theme(text = element_text(size=12), 
        plot.subtitle=element_text(size=8,  face="italic", color="darkgrey"),
        legend.position='none', aspect.ratio=1) 

print(plt.mac)
dev.off()


message("+-------------------------------------------------------------------+")
message("+--- 0)  Pairwise DESeq Analysis, 2 batches                    -----+")
message("+-------------------------------------------------------------------+")

samT.mac       <- as.data.frame(samT.mac, stringsAsFactors = TRUE)
control.names  <- rep("cWD",2)
compare.names  <- rep("iWD",2)

control.Tnames <- c("Mac_N1", "Mac_N2")
compare.Tnames <- c("Mac_A1", "Mac_A2")

selGroups      <- c("Mac1", "Mac2")
rep.numR       <- c(3,3)
samT.mac$Group.Names <- rep(c("Mac1", "Mac2"), c(6,6))

for(i in c(1,2)){
  print(i)
  samT.sub     <- subset(samT.mac, Group.Names==selGroups[i])
  samT.sub$Condition        <- relevel(samT.sub$Condition , control.names[i])
  mat.sub      <- cmat.mac[, colnames(cmat.mac)%in%samT.sub$SampleName]
  colnames(mat.sub)  <- as.factor(colnames(mat.sub))
  dds.sub      <- DESeqDataSetFromMatrix(countData = mat.sub,
                                         colData = samT.sub,
                                         design= ~ Condition)
  dds.sub      <- estimateSizeFactors(dds.sub)
  sizeFactors(dds.sub)

  dds.sub      <- DESeq(dds.sub, parallel=TRUE)
  vst.sub      <- vst(dds.sub,     blind=F)
  resultsNames(dds.sub)
  colData(vst.sub)
  save(mat.sub, samT.sub, dds.sub, vst.sub, 
       file = paste0("./Data/Processed_Data/",selGroups[i], "_DESeq2_dds_vstN.RData"))
  gc()
}

message("+---           Pairwise PCAplot                                -----+")

## customised pca plot function
customPCA     <- function(sampleTBL, RLD, TOPNUM, model) {
  RLD    <- as.data.frame(RLD)
  rv     <- rowVars(as.matrix(RLD))
  select <- order(rv, decreasing = TRUE)[seq_len(min(TOPNUM, length(rv)))]
  pca    <- prcomp(t(RLD[select, ]))
  
  pc1var <- round(summary(pca)$importance[2,1]*100, digits=1)
  pc2var <- round(summary(pca)$importance[2,2]*100, digits=1)
  pc1lab <- paste0("PC1 (",as.character(pc1var),"%)")
  pc2lab <- paste0("PC2 (",as.character(pc2var),"%)")
  
  scores    <- data.frame(sample=sampleTBL$SampleName, pca$x, 
                          condition=sampleTBL$Condition)
  
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

for(i in c(1,2)){
  load(paste0("./Data/Processed_Data/",selGroups[i], "_DESeq2_dds_vstN.RData"))
  pdf(paste0("./Figures_Tables/PCA_Figures/DESeq2_pairwisePCA_", selGroups[i], "_",
             compare.names[i], "_", control.names[i], "_N.pdf"), 
      width=8, height= 6)
  pca <- customPCA(samT.sub, assay(vst.sub), TOPNUM, selGroups[i])
  plot_grid(pca[[1]], pca[[2]], pca[[3]], pca[[4]], nrow=2, byrow=T)
  dev.off()
}

message("+----        Differential analysis, batch1&2                    ----+")

## merge results files and ensembl annotation function. 
merECRes_function <- function(resfile, ensEMBL2id, ddsfile, colNames, outfile){
  resdat   <- as.data.frame(resfile)
  resdat$external_gene_name <- rownames(resdat)
  ### remove the padj missing
  resdat   <- subset(resdat, !is.na(padj)) 
  resmerE  <- merge(resdat, ensEMBL2id, by = "external_gene_name")
  colnames(ddsfile) <- colNames
  normCnt  <- as.data.frame(counts(ddsfile, normalized=TRUE))
  normCnt$ensembl_gene_id <- rownames(normCnt)
  resmerEC <- merge(resmerE, normCnt, by = "external_gene_name")
  write.csv(resmerEC, file=outfile, row.names=F)
}

## given the input and output files

## input
DES.files <- paste0("./Data/Processed_Data/",selGroups[1:2], "_DESeq2_dds_vstN.RData")

## output              
RES.files  <- paste0("./Figures_Tables/DEGs_summaryTables", selGroups[c(1,2)], "_res_all_merE_summaryN.csv")
RESsig.files  <- paste0("./Figures_Tables/DEGs_summaryTables", selGroups[c(1,2)], "_res_Sig_merE_summaryN.csv")

## columns selected
colNames <- list(paste0(rep(c("iWD","cWD"), each=3), "_", c(1,2,3,1,2,3)),
                 paste0(rep(c("iWD","cWD"), each=3), "_", c(1,2,3,1,2,3)))


for(i in 1:length(DES.files)){
  load(DES.files[i])
  res       <- results(dds.sub, contrast = c("Condition", compare.names[i], control.names[i]))
  res.sig   <- subset(res, padj <= significance & abs(log2FoldChange) >= l2fc)
  print(nrow(res.sig)); print(table(sign(res.sig$log2FoldChange)))
  write.csv(res.sig, file = RESsig.files[i])
  res.merEC <- merECRes_function_v2(res, ensEMBL2id, dds.sub, colNames[[i]], RES.files[i])
}

message("+-------------------------------------------------------------------+")
message("+---   1)     merge Mac1 & Mac2 DEGs analysis                  -----+")
message("+-------------------------------------------------------------------+")

## merge raw counts & generate merged sampleTable

mac.merMat <- merge.cmatGN.new[,grep("Mac", colnames(merge.cmatGN.new))]
mac.merSam <- samT.mer[samT.mer$Group=="Mac",]
mac.merSam$Group.Names <- rep(c("Mac1", "Mac2"), c(6,6))

## DEs analysis
dds.macmer <- DESeqDataSetFromMatrix(countData = mac.merMat,
                                     colData = mac.merSam,
                                     design= ~ Group.Names + Condition)
dds.macmer <- estimateSizeFactors(dds.macmer)
sizeFactors(dds.macmer)

dds.macmer <- DESeq(dds.macmer, parallel=TRUE)
vst.macmer <- vst(dds.macmer,     blind=F)
resultsNames(dds.macmer)
colData(vst.macmer)

save(mac.merMat, mac.merSam, dds.macmer, vst.macmer, 
     file ="./Data/Processed_Data/Macmer_DESeq2_dds_vst.RData")

## pca plot detection, with/without limma remove batch effects.

pca1 <- plotPCA(vst.macmer, intgroup="Condition")
mat  <- assay(vst.macmer)
mm   <- model.matrix(~Condition, colData(vst.macmer))
mat1  <- limma::removeBatchEffect(mat, batch=vst.macmer$Group.Names, design=mm)
assay(vst.macmer) <- mat1
pca2 <- plotPCA(vst.macmer, intgroup="Condition")

pdf("./Figures_Tables/PCA_Figures/CAD_mt709_0001-Mac_merge_PCA_rmBatch_plot_27_Jan_2022.pdf")

RLD1   <- mat; RLD2   <- mat1
model1 <- "mer.nrm"; model2 <- "mer.rm"
pca1.cust <- customPCA(mac.merSam, RLD1, TOPNUM, model1)
pca2.cust <- customPCA(mac.merSam, RLD2, TOPNUM, model2)
plot_grid(pca1.cust[[1]], pca1.cust[[2]], pca1.cust[[3]], pca1.cust[[4]], nrow=2, byrow=T)
plot(pca2.cust[[2]])
plot_grid(pca2.cust[[1]], pca2.cust[[2]], pca2.cust[[3]], pca2.cust[[4]], nrow=2, byrow=T)

dev.off()

## hierarchy clustering plot 

pdf("./Figures_Tables/PCA_Figures/Mac_merge_hclust_rmBatch_plot.pdf")
rv_all  <- rowVars(mat)
o_all   <- order(rv_all,decreasing=TRUE)
dists_all <- dist(t(mat[head(o_all,TOPNUM),]))
hc_all    <- hclust(dists_all)
plot(hc_all, labels=vst.macmer$SampleName)
dev.off()


message("+-----                    DEGs identification               --------+") 

colnames.comb <- colData(dds.macmer)$SampleName
REScomb.file <- "./Figures_Tables/Macrophage/Mac_Mer_Result_res_merE_summary.csv"
REScomb.sigfile <- "./Figures_Tables/Macrophage/Mac_Mer_Result_res_Sig_merE_summary.csv"

## use result function 
resmer       <- results(dds.macmer, contrast = c("Condition", "iWD", "cWD"))
resmer.sig   <- subset(resmer, padj <= significance & abs(log2FoldChange) >= l2fc)
resmer.sig   <- resmer.sig[order(-resmer.sig$log2FoldChange),]
print(nrow(resmer.sig)); print(table(sign(resmer.sig$log2FoldChange)))
write.csv(resmer.sig, file = REScomb.sigfile)
resmer.merEC  <- merECRes_function_v2(resmer, ensEMBL2id, dds.macmer, colnames.comb, REScomb.file)


message("+---  Volcano plot for merge Macrophage res files     --------------+")

resfile   <- read.csv(REScomb.file,header=T)
sig_cut   <- 0.05
logfc_cut <- 1 
title     <- " "  
xrange    <- c(-12, 12, 4)
yrange    <- c(0, 4, 2)
topN      <- 15
xlabel    <- "log2FC (iWD/cWD)"
ylabel    <- "-log10 (padj)"


selmarkers1 <- c("Phactr1", "Slco2b1", "Nrp1","Il10ra","Ednrb", "S1pr1")
selmarkers2 <- c("Cadm1",  "Wdfy4", "Bhlhe40")
selmarkers3 <- c("Vps33a", "Trib1")

## customised volcano plot.

functionPlotDEVolcano_selM <- function(resfile, selmarkers1, selmarkers2, 
                                       selmarkers3, sig_cut, logfc_cut, title,  
                                       xrange, yrange,  xlabel, ylabel) {
  options(ggrepel.max.overlaps = Inf)
  results       <- as.data.frame(resfile)
  results$genes <- results$external_gene_name
  results       <- subset(results, !is.na(padj))
  results       <- results[order(-results$log2FoldChange),]
  #selM          <- results[results$gene%in%selmarkers==T,]
  DWsigN        <- sum(results$log2FoldChange<=-1 & results$padj < 0.05)
  UPsigN        <- sum(results$log2FoldChange>= 1 & results$padj < 0.05)
  results$Markers <- 0
  results$Markers <- ifelse(results$gene%in%selmarkers1==T, 1, results$Markers)
  results$Markers <- ifelse(results$gene%in%selmarkers2==T, 2, results$Markers)
  results$Markers <- ifelse(results$gene%in%selmarkers3==T, 3, results$Markers)
  results$colors  <- "grey"
  results$colors  <- ifelse(results$Markers==1, "darkred", results$colors)
  results$colors  <- ifelse(results$Markers==2, "darkorchid", results$colors)
  results$colors  <- ifelse(results$Markers==3, "black", results$colors)
  
  selM1         <- results[results$gene%in%selmarkers1==T,]
  selM2         <- results[results$gene%in%selmarkers2==T,]
  selM3         <- results[results$gene%in%selmarkers3==T,]
  selM          <- results[results$gene%in%c(selmarkers1, selmarkers2, selmarkers3)==T,]
  selMH         <- results[results$gene%in%c(selmarkers2, selmarkers3)==T,]
  colors        <- selM$colors
  volc.plt <- ggplot(data=results, aes(x=log2FoldChange, y=-log10(padj), label=genes, fontface="bold")) +
    geom_vline(xintercept = logfc_cut,     colour="black", linetype = "dashed", alpha=0.5) +
    geom_vline(xintercept = -(logfc_cut),  colour="black", linetype = "dashed", alpha=0.5) +
    geom_hline(yintercept = -log10(sig_cut), colour="black", linetype = "dashed", alpha=0.5) +
    geom_point(data=subset(results, abs(log2FoldChange) < logfc_cut | padj > sig_cut), alpha=0.75, size=0.3, colour="grey") +
    geom_point(data=subset(results, padj<=sig_cut & log2FoldChange >= logfc_cut),      alpha=0.75, size=0.5, colour="darkgoldenrod1") +
    geom_point(data=subset(results, padj<=sig_cut & log2FoldChange <= -(logfc_cut)),   alpha=0.75, size=0.5, colour="cyan3") +
    geom_point(data=selM3,   alpha=0.75, size=1.5, colour="black") +
    geom_point(data=selM1,   alpha=0.75, size=1.5, colour="darkred") +
    geom_point(data=selM2, alpha=0.75, size=1.5, colour="darkorchid") +
    #geom_text_repel(data= selM1,show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=6) +
    geom_text_repel(data= selM,show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=7) +
    xlab(xlabel) + ylab(ylabel) +
    scale_x_continuous(limits=c(xrange[1],xrange[2]), breaks=seq(xrange[1],xrange[2],xrange[3])) +
    scale_y_continuous(limits=c(yrange[1],yrange[2]), breaks=seq(yrange[1],yrange[2],yrange[3])) +
    theme(aspect.ratio=1) +
    ggtitle(title) +
    theme_update(plot.title = element_text(size=16, face="bold", hjust=0.5),
                 axis.title.x = element_text(size=20, face= "bold"),
                 axis.text.x = element_text(size=20, face="bold"),
                 axis.title.y.left = element_text(size=20, face= "bold"),
                 axis.text.y = element_text(size=20, face="bold")) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "white", colour = NA)) 
  
  
  return(volc.plt)
  
}

message("+------                     Fig 2b                           -------+")

mermacRes_selM.new <- functionPlotDEVolcano_selM(resfile, selmarkers1, selmarkers2,
                                                 selmarkers3, sig_cut, logfc_cut, title, 
                                                 xrange, yrange,  xlabel, ylabel)
pdf("./Figures_Tables/Macrophage/Fig2b-MacMer_VolcanoPlot_selmarkers_05_Jan_2023.pdf", 
    onefile = T)
mermacRes_selM.new
dev.off()

message("+-------------------------------------------------------------------+") 
message("+------  2)    GeneOntology and GSEA analysis                 ------+") 
message("+-------------------------------------------------------------------+") 

## load ensembl file with entrezID
load("./Data/References_Data/GRCm39/Ensembl_GRCm39_entrezID_ensembID_exterName.RData")

## merge res file&prepare GO input data for clusterProfiler.
sel.column.names <-c("entrezgene_id", "ensembl_gene_id.x", "external_gene_name", "log2FoldChange")
resall           <- read.csv(REScomb.files[1], header=T)
resall.sig       <- subset(resall, padj < 0.05 & abs(log2FoldChange)>=1)
resdf.ann        <- merge(resall.sig, ensEMBL2id.GO, by = "external_gene_name", all.x=T)
resdf            <- resdf.ann[,colnames(resdf.ann)%in%sel.column.names==T]
resdf            <- resdf[order(-resdf$log2FoldChange),sel.column.names]
resdf            <- resdf[-which(duplicated(resdf$ensembl_gene_id)==T),]
colnames(resdf)  <- c("ENTREZID", "ENSEMBL", "SYMBOL", "L2FC")  ## 746
resdf            <- resdf[-c(610:611), ] ## duplicated Ndor1 gene ## 744

message("+--------GO over-representation analysis               -------------+")

## GO: BP, CC and MF
resbk     <- subset(resall, !is.na(log2FoldChange))
resbk.sub <- resbk[,c("external_gene_name", "ensembl_gene_id", "log2FoldChange")]
length(unique(resbk.sub$external_gene_name)) ## 33787
resdf.bk <- merge(resbk.sub, ensEMBL2id.GO, by = "ensembl_gene_id")
resdf.bk <- resdf.bk[-which(duplicated(resdf.bk$external_gene_name.x)==T), ]
resdf.bk <- resdf.bk[,c(5,1,2,3)]
colnames(resdf.bk) <- c("ENTREZID", "ENSEMBL", "SYMBOL", "L2FC")
resdf.bk <- resdf.bk[order(-resdf.bk$L2FC),] ## 33786

egoAll <- enrichGO(gene = resdf$SYMBOL, universe = resdf.bk$SYMBOL,
                   keyType = "SYMBOL", OrgDb = org.Mm.eg.db,
                   ont = "All", pAdjustMethod = "BH", pvalueCutoff  = 0.05)
egoAll.dat <- as.data.frame(egoAll)

## KEGG
gene <-  as.character(resdf$ENTREZID[-which(is.na(resdf$ENTREZID))])
bkuni.gene <- as.character(resdf.bk$ENTREZID[-which(is.na(resdf.bk$ENTREZID))])
egoKegg <- enrichKEGG(gene = gene, universe=bkuni.gene, organism = 'mmu', pvalueCutoff = 0.05)
egoKegg <- setReadable(egoKegg, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
egoKegg.dat <- as.data.frame(egoKegg)
egoKegg.dat$ONTOLOGY <- "KEGG"
egoKegg.dat <- egoKegg.dat[,c(10,1:9)]

## merge GO and KEGG
egodat <- rbind(egoAll.dat, egoKegg.dat)

genes.split <- lapply(egodat$geneID, function(x) strsplit(x, split="/"))
upgenes  <- list()
dwgenes  <- list()
dirgenes <- NULL
for(i in 1:length(genes.split)){
  subdat  <- resdf[resdf$SYMBOL%in%genes.split[[i]][[1]], ]
  print(dim(subdat)); print(length(genes.split[[i]][[1]]))
  upgenes[[i]]  <- paste0(subdat[subdat$L2FC>0, "SYMBOL"], collapse="/")
  dwgenes[[i]]  <- paste0(subdat[subdat$L2FC<0, "SYMBOL"], collapse="/")
  dirGsub <- sum(subdat$L2FC>0)-sum(subdat$L2FC<0)
  print(dirGsub)
  dirGsub <- ifelse(dirGsub > 0, "+", "-")
  dirgenes <- c(dirgenes, dirGsub)
}
UPGene     <- as.character(upgenes)
DWGene     <- as.character(dwgenes)
DirPath    <- dirgenes
egodat_Direc <- cbind(egodat, UPGene, DWGene, DirPath)

write.csv(egodat_Direc, 
          file = "./Figures_Tables/Macrophage/GoPathway/Mac_Merge_enrichGO_enrichKEGG_withDirection_30102023.csv")

message("+-------------  GO Gene Set Enrichment Analysis        -------------+")

geneList        <- resdf$L2FC
names(geneList) <- resdf$ENTREZID
geneList        <- geneList[-which(is.na(names(geneList)))]

gseBP  <- gseGO(geneList = geneList, OrgDb = org.Mm.eg.db,
                ont = "BP", minGSSize = 10, maxGSSize = 500, 
                pvalueCutoff = 0.05, verbose = FALSE)
gseBP  <- setReadable(gseBP, OrgDb = org.Mm.eg.db, keyType="ENTREZID")

genes.split <- lapply(gseBP$core_enrichment, function(x) strsplit(x, split="/"))
upgenes  <- list()
dwgenes  <- list()
dirgenes <- NULL
for(i in 1:length(genes.split)){
  subdat  <- resdf[resdf$SYMBOL%in%genes.split[[i]][[1]], ]
  print(dim(subdat)); print(length(genes.split[[i]][[1]]))
  upgenes[[i]]  <- paste0(subdat[subdat$L2FC>0, "SYMBOL"], collapse="/")
  dwgenes[[i]]  <- paste0(subdat[subdat$L2FC<0, "SYMBOL"], collapse="/")
  dirGsub <- sum(subdat$L2FC>0)-sum(subdat$L2FC<0)
  print(dirGsub)
  dirGsub <- ifelse(dirGsub > 0, "+", "-")
  dirgenes <- c(dirgenes, dirGsub)
}
UPGene     <- as.character(upgenes)
DWGene     <- as.character(dwgenes)
DirPath    <- dirgenes
gseBP.dat <- as.data.frame(gseBP)
gseBP.dat <- cbind(gseBP.dat, UPGene, DWGene, DirPath)

write.csv(gseBP.dat, file = "./Figures_Tables/Macrophage/GoPathway/Mac_Merge_gseBP_BP3_30102023.csv")
## no enrichment for CC and MF and KEGG

write.csv(resdf, 
          file = "./Figures_Tables/Macrophage/GoPathway/Mac_Mer_resdf_GO_input_16_Dec_2021.csv")
write.csv(resdf.bk, 
          file = "./Figures_Tables/Macrophage/GoPathway/Mac_Mer_resdf_GO_univers_input_16_Dec_2021.csv")

save(egoAll, egoKegg, gseBP, 
     file ="./Figures_Tables/Macrophage/GoPathway/MacMer_pathway_GO_KEGG_groupUniGSE_31102023.RData")


message("+---- Generate genes and pathways barplot for selected pathways ----+")

load("./Figures_Tables/Macrophage/GoPathway/MacMer_pathway_GO_KEGG_groupUniGSE_Jan_2022.RData")
selBPs <- c("macroautophagy", "process utilizing autophagic mechanism",
            "RNA splicing", "mRNA processing",
            "proteasome-mediated ubiquitin-dependent protein catabolic process",
            "response to virus", "regulation of response to biotic stimulus",
            "regulation of innate immune response",
            "ribonucleoprotein complex biogenesis",
            "ncRNA metabolic process")
resfile   <- read.csv(REScomb.file,header=T)
sig_cut   <- 0.05
logfc_cut <- 1

list_up   <- list()
list_down <- list()
egoALL.dat <- as.data.frame(egoAll)
macmer.sig <- resfile
macmer.sig$gene <- macmer.sig$external_gene_name
ego_BP_selected <- egoALL.dat[egoALL.dat$Description %in% selBPs,]
for (i in 1:length(selBPs)){
  df_tmp <- macmer.sig[macmer.sig$gene %in% unlist(strsplit(ego_BP_selected[i, "geneID"] , split="/")),]
  tmp_up <- length(subset(df_tmp[,"log2FoldChange"], df_tmp[,"log2FoldChange"] > 0))
  list_up[[i]] <- tmp_up
  tmp_down <- length(subset(df_tmp[,"log2FoldChange"], df_tmp[,"log2FoldChange"] < 0))
  list_down[[i]] <- tmp_down
}
ego_BP_selected$genes_UP <- as.numeric(list_up)
ego_BP_selected$genes_DOWN <- -as.numeric(list_down)
ego_BP_selected$Description <- c("process utilizing \nautophagic mechanism",
                                 "macroautophagy",
                                 "proteasome-mediated ubiquitin-\ndependent protein catabolic process",
                                 "response to virus", "RNA splicing",
                                 "ribonucleoprotein complex biogenesis",
                                 "regulation of response \nto biotic stimulus",
                                 "ncRNA metabolic process","mRNA processing",
                                 "regulation of innate immune response")
ego_BP_molten <- melt(ego_BP_selected[,c(3,7,8,9, 10:12)],
                      id.vars=c("Description","p.adjust","qvalue","geneID","Count") )

write.table(ego_BP_molten, 
            file = "./Figures_Tables/Macrophage/Fig2a-MacMer_GOBP_sel_Barplot_Up_DW_Apr_2024.txt")

col_1st       <- "darkorange"
col_2nd       <- "steelblue4"
p_bp_mlt <- ggplot(ego_BP_molten, aes(x=reorder(Description, -p.adjust), y=value, fill=variable)) + 
  geom_bar(stat="identity", aes(alpha = -p.adjust), width = 0.8)+
  coord_flip()+
  scale_fill_manual(values = c(col_1st, col_2nd)) +
  ylab("Gene counts") +
  xlab("") +
  ylim(-max(abs(ego_BP_molten$value)), max(abs(ego_BP_molten$value)))+
  
  theme(axis.text.x = element_text(angle = 0, hjust = 1, size=20, face= "bold"),
        axis.title.x = element_text(size=20, face= "bold"),
        axis.title.y = element_text(size=20, face= "bold"),
        axis.text.y = element_text(size=20, face="bold"),
        legend.title = element_text(color = "blue", size = 16, face="bold"), 
        legend.text = element_text(color = "red",size = 16, face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = "white", colour = NA)) +
  guides(fill=guide_legend("Direction"))

message("+---------  Fig2a, select Biological process for barplot -----------+")

pdf("./Figures_Tables/Macrophage/Fig2a-MacMer_GOBP_sel_Barplot_Up_DW_Apr_2024.pdf", 
    width = 12, height= 5, onefile=T)
p_bp_mlt
dev.off() 


message("+-------------------------------------------------------------------+") 
message("+-------- 3)     Transcription Factor Analysis         -------------+") 
message("+-------------------------------------------------------------------+") 

message("+---Part I: TF motif binding enrichments for SigDEGs&nonSigDEGs-----+")

## a) loading the useful libraries relating to TFs analysis

#TF.packages <- c("CeTF", "enrichTF", "target", "TFBSTools", 
# "RcisTarget","GENIE3", "AUCell")
#BiocManager::install(TF.packages, force = T)
library("devtools")
#devtools::install_github("aertslab/SCENIC") 
library("CeTF")
library("enrichTF")
library("target")
library("TFBSTools")
library("RcisTarget")
library("GENIE3")
library("AUCell")
library("SCENIC")
library("zoo")


data(motifAnnotations_mgi)
motifAnnotations_mgi <- motifAnnotations
table(motifAnnotations_mgi$annotationSource)
## directAnnotation             inferredBy_MotifSimilarity 
## 6206                                  59437 
## inferredBy_MotifSimilarity_n_Orthology                   inferredBy_Orthology 
## 159625                                  26858


## b) load the motif binding enrichments files (TSS+/-10kb &  Gene+500bp, -100bp)

dbFiles <- c("https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather",
             "https://resources.aertslab.org/cistarget/databases/mus_musculus/mm10/refseq_r80/mc9nr/gene_based/mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.genes_vs_motifs.rankings.feather")
# mc9nr: Motif collection version 9: 24k motifs
for(featherURL in dbFiles)
{
  download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
}

motifRankings    <- importRankings("./Data/References_Data/TFs/mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.genes_vs_motifs.rankings.feather")
## Number of genes: 24068 (24068 available in the full DB)
## Number of MOTIFS: 24453

## c) Define the sigDEGs(padj<=0.05&abs(log2FC)>=1) 
##    nonsigDEGs (padj>0.1&abs(log2FC)<1)

resall     <- read.csv(REScomb.file, header=T)
DEs_res  <- as.data.frame(resall)
DEs_dat  <- subset(DEs_res, !is.na(log2FoldChange)) ## 
nonDEs_dat  <- subset(DEs_dat, abs(log2FoldChange) < 1 & padj > 0.1) 
## nonsig 3461
nonDEs      <- unique(nonDEs_dat$external_gene_name) 
## unique gene name 3113.
## duplication gene names with different ensemble id. Total 348.
## 5_8s_rRNA(1), 7SK(329), Aldoa(1), Atp50(1), Ddit3(1), 
## Fam220a(1) Meta20a_SRP(11), Pakap(2), Ptp4a1(1),

sigDEs_dat  <- subset(DEs_dat, abs(log2FoldChange) >=1 & padj <= 0.05) # 746
sigDEs      <- unique(sigDEs_dat$external_gene_name) 
## sig 744, where Nord1 has three Ensembl_gene_ids.

length(unique(nonDEs)); length(unique(sigDEs))

## d) TFs motif binding customised function 

TFs_fn <- function(GroupName, inputGenes, motifRankings,  inputdate){
  # Load gene sets to analyze. e.g.:
  geneListN  <- as.character(inputGenes)
  geneListN  <- geneListN[geneListN %in% colnames(motifRankings@rankings)] 
  geneListsN  <- list(geneListName=unique(geneListN))  
  
  motifs_AUC_NDE <- calcAUC(geneListsN, motifRankings, nCores=1)
  motifEnrichmentTable_NDE <- addMotifAnnotation(motifs_AUC_NDE, nesThreshold=3,
                                                 motifAnnot=motifAnnotations_mgi,
                                                 motifAnnot_highConfCat=c("directAnnotation"), 
                                                 motifAnnot_lowConfCat=c("inferredBy_MotifSimilarity",
                                                                         "inferredBy_MotifSimilarity_n_Orthology",
                                                                         "inferredBy_Orthology"),
                                                 digits=3)
 
  motifEnrichmentTable_wGenes_NDE <- addSignificantGenes(motifEnrichmentTable_NDE,
                                                         rankings=motifRankings, 
                                                         geneSets=geneListsN) 
  submotifEnrichmentTable_wGenes_NDE<- subset(motifEnrichmentTable_wGenes_NDE, TF_highConf!="") 
  
  
  ## AUC histogram plot for promoter
  auc    <- getAUC(motifs_AUC_NDE)[1,];  nes3 <- (3*sd(auc)) + mean(auc)
  aucdat <- as.data.frame(auc); rownames(aucdat) <- names(auc); colnames(aucdat) <- "auc"

  pdf(paste0("./Figures_Tables/Macrophage/TFs/Macmer_", GroupName, "_TFs_NES_AUC_hist_", inputdate, ".pdf"))
  hist1 <- ggplot(aucdat, aes(x=auc)) + 
    geom_histogram(color="blue", fill="white", bins=100) +
    geom_vline(aes(xintercept=nes3), color="red", linetype="dashed", linewidth=1.2)

  print(hist1)
  dev.off()
  
  ## add motif logo into the data sets.
  motifEnrichmentTable_wGenes_wLogo_NDE    <- addLogo(motifEnrichmentTable_wGenes_NDE)

  submotifEnrichmentTable_wGenes_wLogo_NDE<- subset(motifEnrichmentTable_wGenes_wLogo_NDE, TF_highConf!="") 
  return(motifEnrichmentTable_wGenes_wLogo_NDE)
  
}


sigTFs    <- TFs_fn(GroupName="SigDEGs", inputGenes=sigDEs, motifRankings, inputdate="Sep_2022")
nonsigTFs <- TFs_fn(GroupName="NonSigDEGs", inputGenes=nonDEs, motifRankings, inputdate="Sep_2022")


submotifEnrichmentTable_wGenes_wLogo_NDE<- subset(nonsigTFs[[1]], TF_highConf!="") 
submotifEnrichmentTable_wGenes_wLogo_DE<- subset(sigTFs[[1]], TF_highConf!="") 


 ## Find overlap TFs motif binding between sigDEGs and nonsigDEGs

TFssig    <- submotifEnrichmentTable_wGenes_wLogo_DE
TFsnonsig <- submotifEnrichmentTable_wGenes_wLogo_NDE
commonTFs_sigNsig <- TFssig[TFssig$motif%in%TFsnonsig$motif, -10] ## 73
commonTFs_data    <- merge(TFssig, TFsnonsig, by="motif")
commonTFs_data  <- commonTFs_data[,c(1,3,5,6,7,8,9,10,13,14,17,18,19)]
colnames(commonTFs_data) <- c("motif", "logo", "NES_sig", "AUC_sig", 
                              "TF_highConf", "nEnrGenes_sig", "rankAtMax_sig",
                              "enrichedGenes_sig", "NES_Nsig","AUC_Nsig", 
                              "nEnrGenes_Nsig", "rankAtMax_Nsig",
                              "enrichedGenes_Nsig")

## Supplementary Tables relating to TFs results with sig/nonsig DEs analysis

write.xlsx(commonTFs_data, 
           file = "./Figures_Tables/Macrophage/TFs/TFs_enrichment_motif_sigNonsigGene_DirectAnnotation_Sep_2022.xlsx",
           sheetName="sigNonsigTFCommon", append=T)

message("+-Supplement Table of TFs binding motif for sig/nonsig DEs analysis-+")

write.xlsx(submotifEnrichmentTable_wGenes_wLogo_DE, 
           file = "./Figures_Tables/Macrophage/TFs/TFs_enrichment_motif_sigNonsigGene_DirectAnnotation_Sep_2022.xlsx",
           sheetName="sigTFsLogo", append=F)
write.xlsx(submotifEnrichmentTable_wGenes_wLogo_NDE, 
           file = "./Figures_Tables/Macrophage/TFs/TFs_enrichment_motif_sigNonsigGene_DirectAnnotation_Sep_2022.xlsx",
           sheetName="NonsigTFsLogo", append=T)

## e) Perform a two sample proportion ztest between the sig and nonsig 
##    common motifs based on the nEnriched_genes

TFs.overlaps <- commonTFs_data
TFs.pvalues  <- NULL
for(i in 1:dim(TFs.overlaps)[1]){
  propztest <- prop.test(x = c(TFs.overlaps[i,"nEnrGenes_sig"], TFs.overlaps[i,"nEnrGenes_Nsig"]), 
                         n = c(dim(sigDEs_dat)[1], dim(nonDEs_dat)[1]))
  subpv  <- propztest$p.value
  TFs.pvalues <- c(TFs.pvalues, subpv)
  TFs.pvalues
}
TFs.adjp <- p.adjust(TFs.pvalues, method="BH")
TFs.stats <- cbind(TFs.overlaps, TFs.pvalues, TFs.adjp)
write.xlsx(TFs.stats, 
           file = "./Figures_Tables/Macrophage/TFs/TFs_enrichment_motif_sigNonsigGene_DirectAnnotation_Sep_2022.xlsx", 
           sheetName = "OverlapTFs_stats_new", append=T)

## f) Perform a NES(the Normalized Enrichment Score) difference between sig and nonsig common TFs.

TFscommon     <- TFs.stats     
TFscommon_sel <- TFscommon[,-c(2,8,13)]
TFscommon_sel$Group <- gsub("(directAnnotation)."," ", TFscommon_sel$TF_highConf)
TFscommon_sel$Group <- gsub("[( .]", "", TFscommon_sel$Group)
TFscommon_sel$Group[58] <- "Gabpa" 
TFscommon_sel <- TFscommon_sel[order(TFscommon_sel$Group),]
TFscommon_sel$Difference <- TFscommon_sel$NES_sig - TFscommon_sel$NES_Nsig
TFscommon_sel <- TFscommon_sel[order(TFscommon_sel$TFs.adjp),]


meanDiff <- aggregate(TFscommon_sel$Difference, list(TFscommon_sel$Group), FUN=mean) 
padjDiff <- aggregate(TFscommon_sel$TFs.adjp, list(TFscommon_sel$Group), FUN=mean)

ss <- TFscommon_sel %>%
  dplyr::mutate(
    Group = fct_relevel(Group, c("Ehf","Elf1","Elf2","Elf3","Elf4","Elf5","Elk1",
                                 "Elk3","Elk4","Erf","Erg","Ets1","Ets2","Etv1",
                                 "Etv2","Etv3","Etv4","Etv5","Etv6","Fev","Fli1",
                                 "Gabpa","Gabpb1","Gm5454","Spib","Spic","Stat1") 
                        ) )


padjDiff <- padjDiff[order(-padjDiff[,2]),]
ss$Group <- factor(ss$Group, levels=padjDiff[,1])
padjDiff$Ts <- -log10(padjDiff[,2])
padjDplt <- as.matrix(padjDiff[,3], ncol=1)
colnames(padjDplt) <- "padj"
rownames(padjDplt) <- padjDiff[,1]

message("+-Ext Fig6a: TFs NEs difference boxplot with padj heatmap alongside-+")

pvalue_col_fun = colorRamp2(c(0, 2, 10, 20, 60), 
                            c("steelblue", "gray", "gold2",  "orange", "darkorange")) 

ha = HeatmapAnnotation(
  pvalue = anno_simple(padjDplt[,1], col = pvalue_col_fun, pch = 0))
grob.ha = grid.grabExpr(draw(ha))
pvalueHeat = Heatmap(padjDplt, col=pvalue_col_fun,cluster_rows =T, show_row_dend = F,
                     width = unit(0.5, "cm"),
                     heatmap_legend_param =list(title = "-log10(padj)",
                                                title_position = "leftcenter-rot",
                                                title_gp = gpar(fontsize = 10, fontface="bold"),
                                                labels_gp = gpar(fontsize = 10, fontface="bold"),
                                                legend_width = unit(3,"cm"),
                                                #legend_height= unit(2, "cm"),
                                                legend_direction = "vertical"
                     ))
grob.heat = grid.grabExpr(draw(pvalueHeat))

pdf("./Figures_Tables/Macrophage/TFs/Extended-Fig6a-TFs_sig_Nsig_NES_difference_boxplot_01112023_new2.pdf", 
    onefile =T, height=8, width=6)
boxplot_DFs <- ggplot(ss, aes(x = Difference, y = Group)) +
  geom_boxplot(alpha=0.5) +
  geom_vline(xintercept = 0,  color = "red", size=1.2) +
  geom_hline(yintercept=c(1:27), linetype="dashed", color = "lightgrey", size=0.2) +
  labs(x = '', y = '') +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=16, face= "bold"),
        axis.title.x = element_text(size=12, face= "bold"),
        axis.title.y = element_text(size=12, face= "bold"),
        axis.text.y = element_text(size=16, face="bold.italic"),
        legend.position="right") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = "white", colour = NA)) 

plot_grid(boxplot_DFs, grob.heat, rel_widths=c(1.2,1))

dev.off()

## d) TF specific "Spic" pathway analysis, relating to Fig2g, here used the sigDEs Spic TFs relating genes

## no KEGG enrichment pathways identified. 
Spic_genes <- TFscommon[grep("Spic",TFscommon$TF_highConf),]
Spic_genes <- strsplit(Spic_genes$enrichedGenes_sig, split= ";")
Mac_sig    <- read.csv("./Figures_Tables/Macrophage/GoPathway/Mac_Mer_resdf_GO_input_16_Dec_2021.csv", header=T)
Spic_res   <- Mac_sig[Mac_sig$SYMBOL%in%Spic_genes[[1]],] 

Spic_GO    <- enrichGO(gene = Spic_res$SYMBOL, OrgDb = org.Mm.eg.db,
                       keyType = "SYMBOL", ont = "ALL",
                       pAdjustMethod = "BH", pvalueCutoff  = 0.05)
Spic_GL    <- Spic_res$L2FC
names(Spic_GL) <- Spic_res$ENTREZID

save(Spic_GO, 
     file = "./Figures_Tables/Macrophage/GoPathway/TF_Spic_N199_enrichGO_BP79_CC12_MF4_Dec_2022.RData")

write.csv(as.data.frame(Spic_GO), 
          file = "./Figures_Tables/Macrophage/GoPathway/TF_Spic_N199_enrichGO_BP79_CC12_MF4_Dec_2022.csv",
          row.names=F, quote=F)

message("+Ext Fig6b: barplot for selected GO pathways from Biological Process in SuppTab5+")
message("+---- select the GO items with at least 12 genes involved for Spic barplot------+")

SpicgoBar <- barplot(Spic_GO, showCategory=c("ribonucleoprotein complex biogenesis",
                                             "autophagy", "rRNA processing",
                                             "response to endoplasmic reticulum stress",
                                             "proteasome-mediated ubiquitin-dependent protein catabolic process",
                                             "actin filament organization")) 

NewBar <- SpicgoBar + scale_fill_continuous() + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=12, face= "bold"),
        axis.title.x = element_text(size=12, face= "bold"),
        axis.title.y = element_text(size=12, face= "bold"),
        axis.text.y = element_text(size=12, face="bold.italic"),
        legend.title = element_text(size = 10, face="bold")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = "white", colour = NA)) 



pdf("./Figures_Tables/Macrophage/GoPathway/Extended-Fig6b-Spic_selectedPathway_Barplot_Dec_2022.pdf",
    width=5, height=4, onefile=T)
NewBar
dev.off()


message("+-------------------------------------------------------------------+") 
message("+--------------------- 4)     CVD GWAS Analysis    -----------------+") 
message("+-------------------------------------------------------------------+") 

message("+---Part I: extract the GWAS data for atherosclerosis relating------+")

library(readr)
library(data.table)
library(openxlsx)
## a) download gwas catalog from EBI online gwas project
GWAS_all    <- read_tsv("./Data/References_Data/GWAS/gwas_catalog_v1.0.2-associations_e106_r2022-07-09.tsv", 
                        col_names=FALSE)
colnames(GWAS_all) <- GWAS_all[1,]
GWAS_all1   <- GWAS_all[-1,]
Trait_Names <- c("atherosclerosis  EFO_0003914", 
                 "coronary atherosclerosis  MONDO_0021661",
                 "coronary atherosclerosis measurement  EFO_0007938", 
                 "carotid atherosclerosis  EFO_0009783",
                 "cardiovascular disease biomarker measurement  EFO_0005278", 
                 "internal carotid intimal medial thickness  EFO_0005053",
                 "common carotid intimal medial thickness  EFO_0004860", 
                 "carotid plaque build  EFO_0006501",
                 "traffic air pollution measurement  EFO_0007908", 
                 "large artery stroke  EFO_0005524",
                 "coronary artery disease  EFO_0001645", 
                 "peripheral arterial disease  EFO_0004265",
                 "Ischemic stroke  HP_0002140", "age at diagnosis  EFO_0004918", 
                 "carotid artery disease  EFO_0003781")
Trait_Forms <- c("EFO_0003914", "MONDO_0021661", "EFO_0007938", "EFO_0009783", 
                 "EFO_0005278", "EFO_0005053", "EFO_0004860", "EFO_0006501",
                 "EFO_0007908", "EFO_0005524", "EFO_0001645", "EFO_0004265", 
                 "HP_0002140", "EFO_0004918", "EFO_0003781")

GWAS_atheroR <- lapply(Trait_Forms, function(x) GWAS_all1[grep(x,GWAS_all1$MAPPED_TRAIT_URI),])

GWAS_atheroRM  <- as.data.frame(unique(rbindlist(GWAS_atheroR))) ## 2595
newdate  <- strptime(as.character(GWAS_atheroRM[,1]), "%d/%m/%Y")
newdate1 <- format(newdate, "%Y-%m-%d")
GWAS_atheroRM[,1] <- newdate1 

## b) Based on the above one, some of them not complete, 
##    I then Download the gwas trait separately and merge 

gwas_arfiles <- list.files(path= "./Data/References_Data/GWAS/GWAS_atherosclerosis", pattern=".tsv")
## the first one is relating to EFO_0005278 CAD
gwas_dat     <- lapply(gwas_arfiles[-1], function(x) 
  read_tsv(paste0("./Data/References_Data/GWAS/GWAS_atherosclerosis/", x), col_names=T))
gwas_datComb <- unique(rbindlist(gwas_dat)) ## 5210
gwas_dat1    <- read_tsv(paste0("./Data/References_Data/GWAS/GWAS_atherosclerosis/", 
                                gwas_arfiles[1]), col_names=F) 
#EFO_0005278 CAD, 32003
colnames(gwas_dat1) <- colnames(gwas_datComb)
gwas_Comb    <- unique(rbind(gwas_datComb, gwas_dat1))  ## 36918

## c) merge all of the gwas together    

GWAS_Athero <- unique(rbind(as.data.frame(gwas_Comb), GWAS_atheroRM)) ## 39513
GWAS_Athero <- GWAS_Athero[order(GWAS_Athero$MAPPED_TRAIT_URI),]
write.table(GWAS_Athero, 
            file = "./Data/References_Data/GWAS/GWAS_atherosclerosis_N39513_04_Jan_2023.txt", 
            col.names=T, row.names=F)
GWAS_Gene   <- GWAS_Athero[,c('PUBMEDID', 'REGION', 'CHR_ID', 'CHR_POS','CONTEXT',
                              'REPORTED GENE(S)','MAPPED_GENE', 
                              'STRONGEST SNP-RISK ALLELE', 'SNPS')]
RegIndex    <- grep(" - ", GWAS_Gene$MAPPED_GENE)
GWAS_Gene$InterGenes[RegIndex] <-GWAS_Gene$MAPPED_GENE[RegIndex] 
GWAS_Gene$UPSTREM     <- unlist(lapply(GWAS_Gene$'InterGenes', 
                                       function(x) strsplit(x, split=" - ")[[1]][1]))
GWAS_Gene$DOWNSTREM   <- unlist(lapply(GWAS_Gene$'MAPPED_GENE', 
                                       function(x) strsplit(x, split=" - ")[[1]][2]))
GWAS_Gene$MGenes      <- GWAS_Gene$MAPPED_GENE
GWAS_Gene$MGenes[RegIndex] <- "NA"

Reported_Genes <- lapply(GWAS_Gene$`REPORTED GENE(S)`, function(x) strsplit(x, split=", "))
Reported_uniG  <- unique(unlist(Reported_Genes))
Reported_uniGO <- Reported_uniG[order(Reported_uniG)]

CheckGenes <- c("ADAMTS-8", "AP000470.2", "APO cluster", "APO-cluster", 
                "APOA1-C3-A4-A5", "APOA1/C3/A4/A5", "apoB", "APOE/C1/C2/C4",
                 "CCDC92-DNAH10", "CDKN2B x CDKN2A", "CHSY-2", "CTC-340I23.2",
                 "CTSB;DEFB13 6", "DKFZp667P0924", "DRB-DQB", "GRB14-COBLL1", 
                 "HFE2 x HFE2", "HOXA cluster", "HSP90B1 x HSP90B1", "intergenic",
                 "intergenic x intergenic", 
                 "intron/Intron/MODIFIER | nc-RNA/Non coding transcript/MODIFIER",
                 "LINC00837)", "LINC01571;C16orf97", "LITAFy", 
                 "LPL - LOC105379311", "Mar-01", "mir-29b-2", "Mir_1302", 
                 "NYAP2-IRS1","PPP1CB/PLB1", "pseudogene", "RBFOX3/ENPP7", 
                 "RBAK–ZNF890P", "RPL6-PTPN11", "SCN5A-SCN10A", "STK32B x STK32B",
                 "TARiD", "MLIP:MLIP-AS1", "NBEA/MIR548F5", "PMF1-BGLAP", "SKCG-1",
                 "SLMO2-ATP5E", "SPECC1L-ADORA2A", "SYNJ2BP-COX16", 
                 "TBC1D7-LOC100130357", "TMX2-CTNND1", "TNFSF12-TNFSF13", 
                 "TRK-TTT15-1", "VTRNA1-3")

CheckIndex <- lapply(CheckGenes, function(x) grep(x,GWAS_Gene$`REPORTED GENE(S)`))
nCIndex    <- lapply(CheckGenes, function(x) grep(x,Reported_uniGO))
ReplaceGenes <- c("ADAMTS8", "AP000470.2", "APOC1,APOC2,APOC4", "APOC1,APOC2,APOC4", 
                  "APOA1,APOC3,APOA4,APOA5", "APOA1,APOC3,APOA4,APOA5", "APOB",
                  "APOE, APOC1, APOC2, APOC4", "CCDC92, DNAH10", "CDKN2B,CDKN2A", 
                  "CHSY3", "LOC105377673", "CTSB,DEFB136","ANKRD36C", 
                  "HLA-DRB1,HLA-DQB1", "GRB14,COBLL1", "HFE2", "HOXA-AS2, HOXA3", 
                  "HSP90B1", "intergenic", "intergenic", "intron", "LINC00837", 
                  "LINC01571, LINC02911", "LITAF", "LPL,LOC105379311", "MARCH1", 
                  "MIR29B2", "MIR1302-8", "NYAP2,IRS1", "PPP1CB,PLB1", "pseudogene", 
                  "RBFOX3, ENPP7", "RBAK,ZNF890P", "RPL6,PTPN11", "SCN5A,SCN10A", 
                  "STK32B", "TARID", "MLIP,MLIP-AS1", "NEBA,MIR548F5", "PMF1-BGLAP", 
                  "SKCG-1", "SLMO2-ATP5E", "SPECC1L-ADORA2A","SYNJ2BP-COX16",
                  "TBC1D7-LOC100130357","TMX2-CTNND1","TNFSF12-TNFSF13", 
                  "TRK-TTT15-1", "VTRNA1-3")

Reported_AllGenes <- Reported_uniGO ## 4350
for(i in 1:length(nCIndex)){
  Reported_AllGenes[unlist(nCIndex[[i]])] <- ReplaceGenes[i]
}
Reported_AllGenes <- unique(unlist(lapply(Reported_AllGenes, 
                                          function(x) strsplit(x, split=",")))) 
## 4336
RMIndex  <- c("intergenic", "intron", "pseudogene", "NA")
Reported_AllGenesRM <- Reported_AllGenes[-which(Reported_AllGenes%in%RMIndex==T)] 
# 4333
Reported_AllGenesRM <- Reported_AllGenesRM[order(Reported_AllGenesRM)]

## d) : Mapped genes/Reported Genes summary   

Mapped_GenesOnly <- unique(unlist(lapply(GWAS_Gene$MGenes, 
                                         function(x) strsplit(x, split=", "))))
Mapped_GenesOnly <- Mapped_GenesOnly[order(Mapped_GenesOnly)]
CheckMG <- c("CDKN2B-AS1 x CDKN2B-AS1", "LPA; LPA; No mapped genes; SLC22A3",
             "LPA; SLC22A3; LPA; No mapped genes", "No mapped genes x HJV", 
             "No mapped genes x STK32B", "No mapped genes x TTC41P", "NA")

MGIndex <- unlist(lapply(CheckMG, function(x) which(Mapped_GenesOnly==x)))
Mapped_GenesOnly[MGIndex] <- c("CDKN2B-AS1", "LPA,SLC22A3", "LPA", "HJV", 
                               "STK32B", "TTC41P", "NA")
Mapped_GenesOnly <- Mapped_GenesOnly[-c(which(Mapped_GenesOnly=="NA"), 
                                        which(is.na(Mapped_GenesOnly)))]  
## 3970
## remove missing 2284 & 3972
Mapped_UPGenes <- unique(GWAS_Gene$UPSTREM[-c(which(GWAS_Gene$UPSTREM=="NA"), 
                                              which(is.na(GWAS_Gene$UPSTREM)))]) 
## 2477
Mapped_UPGenes <- unique(unlist(lapply(Mapped_UPGenes, 
                                       function(x) strsplit(x, split=", "))))
Mapped_UPGenes <- Mapped_UPGenes[order(Mapped_UPGenes)] ## 2490
Mapped_UPGenes <- Mapped_UPGenes[-which(Mapped_UPGenes=="No mapped genes x N/A")]

Mapped_DWGenes <- unique(GWAS_Gene$DOWNSTREM[-c(which(GWAS_Gene$DOWNSTREM=="NA"), 
                                                which(is.na(GWAS_Gene$DOWNSTREM)))]) 
## 2214
Mapped_DWGenes <- unique(unlist(lapply(Mapped_DWGenes, function(x) 
  strsplit(x, split=", "))))
Mapped_DWGenes <- Mapped_DWGenes[order(Mapped_DWGenes)] ## 2214

save(Reported_AllGenesRM, Mapped_GenesOnly, Mapped_UPGenes, Mapped_DWGenes, 
     file = "./Data/References_Data/GWAS/GWAS_Reported_MappedG_MappedUP_MappedDW_Jan_2023.RData")

write.csv(as.matrix(Reported_AllGenesRM, ncol=1), file = "./Data/References_Data/GWAS/GWAS_Reported_Jan_2023.csv")
write.csv(as.matrix(Mapped_GenesOnly, ncol=1), file = "./Data/References_Data/GWAS/GWAS_Mapped_Jan_2023.csv")
write.csv(as.matrix(Mapped_UPGenes, ncol=1), file = "./Data/References_Data/GWAS/GWAS_Mapped_UP_Jan_2023.csv")
write.csv(as.matrix(Mapped_DWGenes, ncol=1), file = "./Data/References_Data/GWAS/GWAS_Mapped_DOWN_Jan_2023.csv")

message("+-Part II: Merge all of the GWAS, convert Human genes to mouse genes-+")

AllGWAS_Relating <- unique(c(Reported_AllGenesRM, Mapped_GenesOnly, 
                             Mapped_UPGenes, Mapped_DWGenes)) ## 8643
AllGWAS_Relating <- AllGWAS_Relating[order(AllGWAS_Relating)]
write.csv(as.matrix(AllGWAS_Relating, ncol=1), 
          file = "./Data/References_Data/GWAS/GWAS_mergeAll_Jan_2023.csv")

## human and mouse gene convert, Error in biomaRt
# Basic function to convert human to mouse gene names
convertHumanGeneList <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("hgnc_symbol"), filters = "hgnc_symbol", 
                   values = x , mart = human, attributesL = c("mgi_symbol"), 
                   martL = mouse, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}
genes <- convertHumanGeneList(AllGWAS_Relating)

## use g:Profiler,https://biit.cs.ut.ee/gprofiler/orth, g:Orth, orthology search

message("+-Part III: Check overlap between CVD relating GWAS and our sig/nonsig DEGs--+")

HMorth   <- read.csv("./Data/References_Data/GWAS/gProfiler_hsapiens_mmusculus_05-01-2023_15-07-39.csv", header=T)
MouseGW  <- unique(HMorth[-which(HMorth[,5]=="N/A"),5]) ## 6616- unique 6040

colnames(HMorth)[5] <- "external_gene_name"
colnames(HMorth)[2] <- "Human_Symbol"
Mac_Res  <- read.csv("./Figures_Tables/Macrophage/Mac_Mer_Result_res_merE_summary.csv", header=T)
Mac_sel  <- subset(Mac_Res, !is.na(padj))   ## 7832 
Mac_sig  <- subset(Mac_sel, padj <= 0.05 & abs(log2FoldChange) >= 1)  ## 746
Mac_nonsig1 <- subset(Mac_sel,  abs(log2FoldChange) <1 &padj > 0.1) ## 3461

Macsig_overGWAS  <- unique(merge(Mac_sig[,1:8],HMorth[,c(2,5)], by="external_gene_name")) ## 211
Macnsig1_overGWAS <- unique(merge(Mac_nonsig1[,1:8],HMorth[,c(2,5)], by="external_gene_name")) ## 7461

Macsig_GW  <- unique(Macsig_overGWAS$external_gene_name) ## 209
Macnsig1_GW <- unique(Macnsig1_overGWAS$external_gene_name) ## 831

## proportion ztest with alternative "greater", sigProption > NsigProportion overlap
M3 <- as.table(rbind(c(209, 831), c(746-209, 3461-831)))
dimnames(M3) <- list(GWAS = c("Y", "N"),
                     SEGs = c("sig","Nonsig"))
M3xsq   <- chisq.test(M3) ## p-value = 0.02422, two tails chisq test
M3test  <- prop.test(x = c(209, 831), n = c(746, 3461), alternative = "greater") ## 0.01211

message("+---- Fig4a: Barplot to show the proportion overlaps ---------------+")

barplot.dat <- data.frame(Group=c("Nonsig\nDEGs", "Sig\nDEGs"),
                          percentage = c(24.80, 28.02))

pdf("./Figures_Tables/Macrophage/Fig4a-GWAS_sig_Nsig_BarPlot_Jan_2023.pdf", onefile =T)

barplot.gr <- ggplot(barplot.dat, aes(x=reorder(Group, percentage), y=percentage, fill= Group)) + 
  geom_bar(stat="identity")+
  scale_fill_manual(values=c("blue", "red")) +
  xlab("") +
  ylab("% Genes with link to \n cornaryartery disease GWAS") +
  ylim(0,30) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size=16, face= "bold"),
        axis.title.x = element_text(size=20, face= "bold"),
        axis.text.y = element_text(size=20, face="bold"),
        axis.title.y = element_text(size=20, face= "bold"),
        legend.position="none") +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        panel.background = element_rect(fill = "white", colour = NA)) + 
  geom_text(x=1.5, y=30, label="P=0.012", size = 6)
barplot.gr
dev.off()


## write the GWAS overlap with different sig/Nsig datasets
write.csv(unique(Macsig_overGWAS[,c(1,3,7)]), 
          file = "./Data/References_Data/GWAS/GWAS_overlap_sig_N209_TN746_Jan_2023.csv", 
          row.names=F)
write.csv(unique(Macnsig1_overGWAS[,c(1,3,7)]), 
          file = "./Data/References_Data/GWAS/GWAS_overlap_Nsig_p01_l2fc1_N831_TN3461_Jan_2023.csv", 
          row.names=F)

message("+Part IV: Check 209 Overlap GWAS sigDEGs assign to different type of Macrophages+")

## singlecell integrating analysis are performed for both mouse and human are separately from this script.
ScRNA_overlapSig <- openxlsx::read.xlsx("./Figures_Tables/Macrophage/Supplementary_Data_version1_Jan_2023.xlsx", sheet=4, startRow=3)
scSig <- list(ScRNA_overlapSig[1:47,1], ScRNA_overlapSig[1:55,2],ScRNA_overlapSig[1:61,3])
gwas_scsig <- lapply(scSig, function(x) x[x%in%Macsig_overGWAS[,1]])

save(ScRNA_overlapSig, gwas_scsig, 
     file = "./Data/References_Data/GWAS/GWAS_MacSig_scRNA_overlap_Res15_IFN19_Air30_Jan_2023.RData")


message("+-Figrue4b: Pie chart for the GWAS overlap with integrated singlecell RNASeq analysis +")

df.pie <- data.frame(
  group = factor(c("Resident/Adventitia", "Inflammation/IFN",   "Trem2/MacAir"),
                 level=c("Resident/Adventitia", "Inflammation/IFN",   "Trem2/MacAir")),
  value = c(15,19,30)
)
data <- df.pie %>% 
  arrange(desc(group)) %>%
  mutate(prop = value / sum(df.pie$value) *100) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

# Basic piechart
pieplt <- ggplot(data, aes(x="", y=prop, fill=group)) +
  geom_bar(stat="identity", width=1, color="white") +
  coord_polar("y", start=0) +
  theme_void() + 
  theme(legend.position="none") +
  geom_text(aes(y = ypos, label = paste(round(prop, digits=1), "%")), color = "white", size=7) +
  scale_fill_manual(values=c("darkred", "grey", "#4B0082"))

pdf("./Figures_Tables/Macrophage/Fig4b-MacMer_sigDEGs_overlapCelltype_gwas_piechart_group3_Jan_2022.pdf",  onefile=T)
pieplt
dev.off() 

message("+-------------------------------------------------------------------+") 
message("+---5) Overlap with Sym/Nonsym analysis, Fernandez_2020_Human pap---+") 
message("+-------------------------------------------------------------------+") 

### Fernandez et al paper. Single-Cell-Immune-Profiling-of-Atherosclerotic-Plaques, Nature Medicine
### Data https://figshare.com/s/c00d88b1b25ef0c5c788
### used the umi normalised data performed the DEGs.
### merged_plaque_gex-umi-data-mnn-cat_macrophages.txt

param = SnowParam(4, "SOCK", progressbar=TRUE)
ndir <- "/Users/xz289/Documents/Minoru/Project_mt709_0001/Fernandez_2020_Human"
macdat  <- read.delim(paste0(ndir, "/merged_plaque_gex-umi-data-mnn-cat_macrophages.txt"))
rownames(macdat) <- macdat[,1]
macdat  <- macdat[,-1]

metadat <- colnames(macdat)
Barcode <- substr(metadat, 4, 26)
Sample  <- substr(metadat, 39, 44)
Condition <- substr(metadat, 62,73); Condition <- gsub("[.]", "", Condition)
CellType  <- substr(metadat, 88,101); CellType <- gsub("[.]", "", CellType)
BroadCT <- substr(metadat, 115, 129); BroadCT  <- gsub("[.]", "", BroadCT)

meta.data <- data.frame(Barcode, Sample, Condition, BroadCT)
rownames(meta.data) <- metadat
write.table(meta.data, 
            file = paste0(ndir, "/merged_plaque_gex-mnn_macrophages_metadata.txt"))


message("+--- Apply the two-sided Welch’s t-test in R plus BH correction-----+")

Macro_symInd  <- rownames(meta.data)[meta.data$Condition=="Symptomatic"]
Macro_sym <- macdat[,colnames(macdat)%in%Macro_symInd]
Macro_AsymInd  <- rownames(meta.data)[meta.data$Condition=="Asymptomatic"]
Macro_Asym <- macdat[,colnames(macdat)%in%Macro_AsymInd]

MacroSym_mean  <- rowMeans(Macro_sym)
MacroAsym_mean <- rowMeans(Macro_Asym)


Pval <- NULL
for(i in 1:length(MacroSym_mean)){
  subP <- t.test(x=Macro_sym[i,], y=Macro_Asym[i,], "two.sided")$p.value
  Pval <- c(Pval, subP)
  Pval
}

Macro_TtestM   <- data.frame(Sym.mean=MacroSym_mean, Asym.mean=MacroAsym_mean, 
                             log2FoldChange=log2(MacroSym_mean/MacroAsym_mean),
                             pval=Pval)
rownames(Macro_TtestM) <- rownames(macdat)
Macro_sumD <- Macro_TtestM[!is.na(Macro_TtestM$log2FoldChange),] ## 38 NA
Macro_sumD <- Macro_sumD[-which(Macro_sumD$log2FoldChange=="Inf"),] 
## 40 Inf, 214 -Inf
Macro_sumD <- Macro_sumD[-which(Macro_sumD$log2FoldChange=="-Inf"),] ## 9708
Macro_sumD <- Macro_sumD[order(-Macro_sumD$log2FoldChange),]
Macro_sumD$padj <- p.adjust(Macro_sumD$pval, method="BH")

write.csv(Macro_sumD, file = paste0(ndir, "/Macro_SymAsym_DEGs_Jan_2023.csv"))
## check overlap with our DEGs
# R Program to illustrate
# the use of str_to_title function
# Loading Library
library(stringr)
Mac_sigOrth  <- read.csv(paste0(ndir, "/gProfiler_mmusculus_hsapiens_11-01-2023_19-17-35.csv"), header = T)
Mac_sigOrth[,2] <- unlist(lapply(Mac_sigOrth[,2], 
                                 function(x) str_to_title(as.character(x))))
test <- subset(Macro_sumD, padj<=0.05)
test$ortholog_name <- rownames(test)
Mac_merD <- merge(Mac_sigOrth, test, by = "ortholog_name")
Mac_merD <- Mac_merD[,c(1,3,4,6,8,9,10,11,12)]
colnames(Mac_merD) <- c("Human_gene", "Mouse_gene", "Mouse_ID", "Human_ID", 
                        "Sym.mean", "Asym.mean",
                        "Sym_Asym.l2fc", "Sym_Asym.pval", "Sym_Asym.padj")
## Mac_merD with Mac_sig
Mac_Res  <- read.csv("./Figures_Tables/Macrophage/Mac_Mer_Result_res_merE_summary.csv", header=T)
Mac_sel  <- subset(Mac_Res, !is.na(padj))[,c(1,3,7)]   ## 7832 
Mac_sig  <- subset(Mac_sel, padj <= 0.05 & abs(log2FoldChange) >= 1) ## 746
colnames(Mac_sig) <- c("Mouse_gene", "Mac_l2fc", "Mac_padj")

Mac_mersigD <- merge(Mac_merD, Mac_sig, by = "Mouse_gene")
Mac_mersigD_sel <- Mac_mersigD[,c(1,2,7,9,10,11)]
Mac_mersigD_sig <- subset(Mac_mersigD_sel, abs(Sym_Asym.l2fc)>=1 & Sym_Asym.padj <=0.05)
write.csv(Mac_mersigD_sel, 
          file = paste0(ndir, "/MacsigN746_over_MacSymAsymSigN811_N22_Jan_2023.csv"), 
          row.names=F)

message("+---      Fig4c: Volcano plot for Sym/Asym                  --------+")

## generate volcano plot.
summary(Macro_sumD$log2FoldChange)
summary(-log10(Macro_sumD$padj))

## customised volcano function
functionPlotDEVolcano_gwas <- function(resfile, selmarkers1, selmarkers2, 
                                       selmarkers3, sig_cut, logfc_cut, 
                                       xrange, yrange, xlabel, ylabel) {
  options(ggrepel.max.overlaps = Inf)
  results       <- as.data.frame(resfile)
  results$genes <- results$X
  results       <- results[order(-results$log2FoldChange),]
  DWsigN        <- sum(results$log2FoldChange<=-1 & results$padj < 0.05)
  UPsigN        <- sum(results$log2FoldChange>= 1 & results$padj < 0.05)
  results$Markers <- 0
  results$Markers <- ifelse(results$gene%in%selmarkers1==T, 1, results$Markers)
  results$Markers <- ifelse(results$gene%in%selmarkers2==T, 2, results$Markers)
  results$colors  <- "grey"
  results$colors  <- ifelse(results$Markers==1, "darkred", results$colors)
  results$colors  <- ifelse(results$Markers==2, "darkmagenta", results$colors)
  results$colors  <- ifelse(results$Markers==3, "gray28", results$colors)
  
  selM1         <- results[results$gene%in%selmarkers1==T,]
  selM2         <- results[results$gene%in%selmarkers2==T,]
  selM3         <- results[results$gene%in%selmarkers3==T,]
  selM          <- results[results$gene%in%c(selmarkers1, selmarkers2, selmarkers3)==T,]
  colors        <- selM$colors
  selMA         <- results[results$gene%in%selmarkers1[1]==T,]
  selMB         <- results[results$gene%in%c(selmarkers1[2], selmarkers2, selmarkers3)==T,]
  volc.plt <- ggplot(data=results, aes(x=log2FoldChange, y=-log10(padj), label=genes, fontface="bold")) +
    geom_vline(xintercept = logfc_cut,     colour="black", linetype = "dashed", alpha=0.5) +
    geom_vline(xintercept = -(logfc_cut),  colour="black", linetype = "dashed", alpha=0.5) +
    geom_hline(yintercept = -log10(sig_cut), colour="black", linetype = "dashed", alpha=0.5) +
    geom_point(data=subset(results, abs(log2FoldChange) < logfc_cut | padj > sig_cut), alpha=0.75, size=0.3, colour="grey") +
    geom_point(data=subset(results, padj<=sig_cut & log2FoldChange >= logfc_cut),      alpha=0.75, size=0.5, colour="darkgoldenrod1") +
    geom_point(data=subset(results, padj<=sig_cut & log2FoldChange <= -(logfc_cut)),   alpha=0.75, size=0.5, colour="cyan3") +
    geom_point(data=selM1,   alpha=0.75, size=3.0, colour="darkred") +
    geom_point(data=selM2, alpha=0.75, size=3.0, colour="darkmagenta") +
    geom_point(data=selM3, alpha=0.75, size=3.0, colour="gray28") +
    geom_text_repel(data= selMA,show.legend = FALSE, nudge_x=-0.3, nudge_y=0.1, segment.size = 0.25, size=12) +
    geom_text_repel(data= selMB,show.legend = FALSE, nudge_x=0.1, nudge_y=0.1, segment.size = 0.25, size=10) +
    xlab(xlabel) + ylab(ylabel) +
    scale_x_continuous(limits=c(xrange[1],xrange[2]), breaks=seq(xrange[1],xrange[2],xrange[3])) +
    scale_y_continuous(limits=c(yrange[1],yrange[2]), breaks=seq(yrange[1],yrange[2],yrange[3])) +
    theme(aspect.ratio=1) +
    theme_update(axis.title.x = element_text(size=26, face= "bold"),
                 axis.text.x = element_text(size=26, face="bold"),
                 axis.title.y.left = element_text(size=26, face= "bold"),
                 axis.text.y = element_text(size=26, face="bold")) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          panel.background = element_rect(fill = "white", colour = NA)) 
  
  
  return(volc.plt)
  
}

selmarkers1 <- c("NRP1",  "PLEC")
selmarkers3 <- c("EIF4E",  "GEM", "NFKB1")
selmarkers2 <- c("CNN2")

Vol_plt_gwas <- functionPlotDEVolcano_gwas(Macro_sumD, selmarkers1, selmarkers2, 
                                           selmarkers3, sig_cut, logfc_cut, 
                                           xrange, yrange, xlabel, ylabel)

pdf("./Figures_Tables/Macrophage/Fig4c-Fernandez_2020_Human_SymAsym_VolcanoPlot_April_2024_gwas.pdf", 
    onefile=T)
Vol_plt_gwas
dev.off()

write.table(Macro_sumD, 
            file = "./Figures_Tables/Macrophage/Fig4c-Fernandez_2020_Human_SymAsym_VolcanoPlot_April_2024_gwas.txt")

message("+-------------------------------------------------------------------+") 
message("+--------6) GWAS causal genes from Adam Butterworth paper   --------+")
message("+--- Supplementary Table 31 applied here.     ----------------------+")
message("+-------------------------------------------------------------------+") 

## a) Run enrichment of gene ontology pathway Biological Process using clusterProfiler. 

library(clusterProfiler)
load("./Data/References_Data/GRCm39/Ensembl_GRCm39_entrezID_ensembID_exterName.RData")
causalG_dat <- read.csv("./AdamB_NG_paper_N220_Pathwayanalysis/ST31_geneList_NG.csv", header=F)[1:220,1]
eGO_CG <- enrichGO(gene = causalG_dat,
                   keyType = "SYMBOL",
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)

save(eGO_CG, file = "./AdamB_NG_paper_N220_Pathwayanalysis/Biological_Process_clusterProfiler_causalGenes_output.RData")

## total 178 genes out of 220(206) contribute to the BP `go` analysis.

xlsx::write.xlsx(eGO_CG, 
                 file = "./AdamB_NG_paper_N220_Pathwayanalysis/Biological_Process_clusterProfiler_causalGenes_output.xlsx",
                 append = F, sheetName = "BP_all220")
eGO_CG_act <- eGO_CG[grep("actin",eGO_CG$Description ),]
xlsx::write.xlsx(eGO_CG_act, 
                 file = "./AdamB_NG_paper_N220_Pathwayanalysis/Biological_Process_clusterProfiler_causalGenes_output.xlsx",
                 append = T, sheetName="Actin_Rel")



categorys <- c("neutral lipid metabolic process", "ameboidal-type cell migration",
  "renal system vasculature development", "actin filament organization")

CGs <- read.csv("./AdamB_NG_paper_N220_Pathwayanalysis/CausalGenes_gProfiler_hsapiens_mmusculus_29-11-2023_12-39-56.csv", header = T)
Mac_Res  <- read.csv("./Figures_Tables/Macrophage/Mac_Mer_Result_res_merE_summary.csv", header=T)
Mac_sel  <- subset(Mac_Res, !is.na(padj))   ## 7832 
Mac_sig  <- subset(Mac_sel, padj <= 0.05 & abs(log2FoldChange) >= 1)  ## 746

selGOs   <- eGO_CG_dat[eGO_CG_dat$Description%in%categorys, "geneID"]
selGOs_genes <- lapply(selGOs, function(x) strsplit(x, split="/")[[1]])
selGOs_ortho <- lapply(selGOs_genes, function(x) CGs[CGs$initial_alias%in%x, "ortholog_name"])
overGO_genes <- lapply(selGOs_ortho, function(x) Mac_sig[Mac_sig$external_gene_name%in%x, c(1,3,6,7)])
overGO_genesNS <- lapply(selGOs_ortho, function(x) Mac_sel[Mac_sel$external_gene_name%in%x, c(1,3,6,7)])

message("+-----    Extended Fig7: network plot for selected pathways --------+")

pdf("./AdamB_NG_paper_N220_Pathwayanalysis/Extended-Fig7-CVD_CasusalGenes_selPathway_cnetplot_112023.pdf",
    width = 10, height =8)
cnetplot(eGO_CG, showCategory = categorys, color_category='firebrick', color_gene='steelblue')
## To make the GO enrichment more clear, we give the cut-off of padj is 0.01
dev.off()


##------------------ Finished after revised on 19/04/2024 -------------------##




