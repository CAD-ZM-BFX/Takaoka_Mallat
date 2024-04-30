#!/usr/local/bin/Rscript
# Four public mouse athero singlecell RNASeq study analysis, Ldlr-/- mouse, male
# chow diet,  HFD different time
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
message("+ General settings and libraries calling, load scRNA objects.       +")
message("+-------------------------------------------------------------------+")

setwd('/Users/xz289/Documents/Minoru/Seurat_v3/Athero_human/')
library(Seurat)
library(ggplot2)
library(viridis)
library(RColorBrewer)
library(dplyr)
library(cowplot)
library(tidyr)
library(tidyverse)
library(xlsx)

load(file="./Alsaigh_ComBio/Alsaigh_merged.Robj")
load(file="./Chowdhury_Circu/Chowdhury_merged.Robj")
load(file="./plaquemac_unfiltered.Robj")
load(file="./Fernandez_NatMed/CITE_seq_Plaque/hplaque_cite_mac.Robj")
load(file="./Li_Circu/Li_merged.Robj")
load(file="./Pan_Circu/pan_merged.Robj")
load(file="./Wirka_NatMed/coronary.Robj")

message("+---merge the 7 data set and recluster to identify the Macrophage --+")

human_merge<- merge(x=Alsaigh_merged, y=list(Chowdhury_merged, 
                                             plaquemac_unfiltered, 
                                             hplaque_cite_mac, Li_merged, 
                                             pan_merged, coronary))

human.merge.list[["RNA"]] <- split(human_merge[["RNA"]], f = human_merge$Protocol)

human.merge.list <- human.merge.list[c("Alsaigh", "Chowdhury", "Fernandez", 
                                   "Fernandez_CITE_MAC", "Li", "Pan", "Wikra")]

head(human.merge.list)

message("+--------      Perform analysis without integration     ------------+")

human.merge.list <- NormalizeData(Ahuman.merge.list, verbose = FALSE)
human.merge.list <- FindVariableFeatures(human.merge.list, 
                                          selection.method = "vst", 
                                          nfeatures = 2000, verbose = FALSE)
human.merge.list <- ScaleData(human.merge.list)
human.merge.list <- RunPCA(human.merge.list)

human.merge.list <- FindNeighbors(human.merge.list, dims = 1:30, 
                                   reduction = "pca")


message("+Using clustree to decide the clusters with different resolutions.+")

human.merge.rlist <- human.merge.list
human.merge.rlist <- FindNeighbors(human.merge.rlist, dims = 1:30)
human.merge.rlist <- FindClusters(human.merge.rlist, resolution = 0)

for(res in c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
{
  message(paste0("+--Method 0 : Finding Clusters at ", res, " resolution ---+"))
  human.merge.rlist <- FindClusters(human.merge.rlist, 
                                     resolution = res, 
                                     print.output = FALSE) 
}

library(clustree)
mat.convert  <- human.merge.rlist@meta.data
mat.convert.numeric <- apply(mat.convert[,grep("integrated_snn_res.", colnames(mat.convert))], 
                             2, as.numeric)
clust.input <- as.data.frame(mat.convert.numeric)

pdf("./Human-Clustree_genecut200_cellcut3_res0_1.pdf", width = 15, height = 20)
clustree(clust.input, prefix ="integrated_snn_res.")
dev.off()

message("+-------- Final decision of res is going to use 0.3           ------+")

human.merge.list <- FindClusters(human.merge.list, resolution = 0.3, 
                                  cluster.name = "unintegrated_clusters")  

human.merge.list <- RunUMAP(human.merge.list, dims = 1:30, 
                             reduction = "pca", 
                             reduction.name = "umap.unintegrated")

DimPlot(human.merge.list, reduction = "umap.unintegrated", 
        group.by = c("Protocol", "seurat_clusters"))


message("+-------------------------------------------------------------------+")
message("+---------     Perform integration analysis                ---------+")
message("+-------------------------------------------------------------------+") 

human.merge.integrated <- IntegrateLayers(object = human.merge.list, 
                                           method = scVIIntegration, 
                                           orig.reduction = "pca", 
                                           new.reduction = "integrated.scvi",
                                           verbose = FALSE)

# re-join layers after integration
human.merge.integrated[["RNA"]] <- JoinLayers(human.merge.integrated[["RNA"]])

human.merge.integrated <- FindNeighbors(human.merge.integrated, 
                                         reduction = "integrated.scvi", 
                                         dims = 1:30)
human.merge.integrated <- FindClusters(human.merge.integrated, 
                                        resolution = 0.3,
                                        cluster.name = "scvi_clusters")

human.merge.integrated <- RunUMAP(human.merge.integrated, 
                                   reduction = "integrated.scvi", 
                                   dims = 1:30, reduction.name = "umap.scvi")
p1 <- DimPlot(
  human.merge.integrated,
  reduction = "umap.scvi",
  group.by = c("Protocol", "scvi_clusters"),
  combine = FALSE, label.size = 2
)

save(human.merge.integrated, 
     file="./human_merge_integrated_scvi_res03.Robj")

message("+-------------------------------------------------------------------+")
message("+---------     Perform annotation analysis                 ---------+")
message("+-------------------------------------------------------------------+")


message("+--------  Find Markers at res0.3 clustering,             ----------+") 

human.merge.integrated.markers<- FindAllMarkers(human.merge.integrated, 
                                                 min.pct = 0.3, 
                                                 logfc.threshold = 0.2, 
                                                 max.cells.per.ident = 250, 
                                                 only.pos = T)##20 clusters

human.merge_scvi_Average <- AverageExpression(human.merge.integrated,
                                              return.seurat = T)

write.csv(human.merge.integrated.markers,
          file="./human.merge.integrated.markers.res03.csv")
save(human.merge_scvi_Average,file="./human.merge_scvi_Average_res03.Robj")

message("+--- Apply mac surface markers to identify macrophage clusters -----+")

human.scvi.integrated <- human.merge.integrated
head(human.scvi.integrated@meta.data)
Idents(human.scvi.integrated)=human.scvi.integrated[["scvi_clusters"]]
levels(human.scvi.integrated)


scvi_HUMAN_MARKERS=FindAllMarkers(human.scvi.integrated)
save(scvi_HUMAN_MARKERS,file="./scvi_HUMAN_MARKERS.Robj")

## Macrophage surface markers list

MAC_genes <- c("CXCL8", "CCL4", "CXCL2", "IER3", "CLEC4E", "IL1B", "CD74", 
               "HLA-DRA", "FABP4", "APOC1", "GPNMB","LGALS3","TREM2", "CD9",
               "APOE", "CTSD", "SPP1", "MARCO", "F13A1", "LGMN", "FOLR2", 
               "CD163", "TIMD4", "LYVE1",  "NRP1", "IFI44L", "XIST", "FCGR3A", 
               "CX3CR1", "IFI6", "EGR1")

test.scvi.integrated <- human.scvi.integrated
DefaultAssay(test.scvi.integrated) <- "RNA"
MAC_marker_list <- list(MAC_genes)
MAC_object <- AddModuleScore(object = test.scvi.integrated, 
                             features = MAC_marker_list, 
                             name = "MAC_marker_score")

pdf("./AFig2a-Human_integrated_all_UMAP_MacGenes.pdf", width=9, height=4, onefile = T)

plt_umap <- UMAPPlot(human.scvi.integrated, label=T) + theme(legend.position="none")
plt_mac  <- FeaturePlot(object = MAC_object, features = "MAC_marker_score1") +
     scale_size(range = c(2, 7)) + theme(plot.title = element_blank())+
     scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
plot_grid(plt_umap, plt_mac, ncol = 2)

dev.off()

## plot the umap split by patient and see how it distributed differently. 
UMAPPlot(human.scvi.integrated, group.by="Patient", split.by="Patient", ncol=4) + 
  scale_color_brewer(palette = "Paired") + 
  NoLegend()

Idents(human.scvi.integrated)="Patient"

## Repartition of clusters across patients using dittoSeq

library(dittoSeq)
dittoBarPlot(human.scvi.integrated, var = Idents(human.scvi.integrated), 
             group.by = "Patients", scale = "percent") + 
  scale_fill_manual(values=c("#1B9E77","#D95F02" ,"#7570B3","#E7298A", "#66A61E",
                             "#E6AB02","#A6761D", "#666666","#2166AC","#810F7C"))

## Combine Zernecke's paper Supplementary Table 3 human MPC markers, 
## Plot the subtype of Macrophages umap plot to select the clusters of macrophage
## from our integrated object.


MAC_genes_Data <- read.xlsx("./cvac161_supplementary_data/Supp_xls_3_Human_MPC_Marker_Genes.xlsx", 
                            sheetIndex = 1, header = T)
colnames(MAC_genes_Data) <- MAC_genes_Data[1,]
MAC_genes_Data <-MAC_genes_Data[-1,]

ResMAC_genes     <- MAC_genes_Data[MAC_genes_Data$cluster=="LYVE1_Mac", "gene"]
FoamyMAC_genes   <- MAC_genes_Data[MAC_genes_Data$cluster=="Foamy_Mac", "gene"] 
C3MAC_genes      <- MAC_genes_Data[MAC_genes_Data$cluster=="C3_Mac", "gene"]
IFNICMAC_genes   <- MAC_genes_Data[MAC_genes_Data$cluster=="IFNIC_Mac", "gene"]
InflamMAC_genes  <- MAC_genes_Data[MAC_genes_Data$cluster=="Inflammatory_Mac", "gene"]

## check the macrophage clusters first
MAC_ALL <- list(unique(c(ResMAC_genes, FoamyMAC_genes, C3MAC_genes, IFNICMAC_genes, InflamMAC_genes)))
MACALL_object <- AddModuleScore(object = test.scvi.integrated, features = MAC_ALL, 
                                name = "MACALL_marker_score")
plt_1 <- FeaturePlot(object = MACALL_object, 
                     features = "MACALL_marker_score1", raster=F) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


## Res-like, LYVE1
MAC_Res <- list(ResMAC_genes)
MACRes_object <- AddModuleScore(object = test.scvi.integrated, features = MAC_Res, 
                                name = "MACRes_marker_score")
plt_Res <- FeaturePlot(object = MACRes_object, 
                       features = "MACRes_marker_score1", raster=F) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


##  Foamy like, TREM2
MAC_Foam <- list(FoamyMAC_genes)
MACFoam_object <- AddModuleScore(object = test.scvi.integrated, features = MAC_Foam, 
                                 name = "MACFoam_marker_score")
plt_Foam <- FeaturePlot(object = MACFoam_object, 
                        features = "MACFoam_marker_score1", raster=F) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

## C3 inflam 
MAC_C3 <- list(C3MAC_genes)
MACC3_object <- AddModuleScore(object = test.scvi.integrated, features = MAC_C3, 
                               name = "MACC3_marker_score")
plt_C3 <- FeaturePlot(object = MACC3_object, 
                      features = "MACC3_marker_score1", raster=F) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


## IFNIC mac
MAC_IFNIC <- list(IFNICMAC_genes)
MACIFNIC_object <- AddModuleScore(object = test.scvi.integrated, 
                                  features = MAC_IFNIC, 
                                  name = "MACIFNIC_marker_score")
plt_IFNIC <- FeaturePlot(object = MACIFNIC_object, 
                         features = "MACIFNIC_marker_score1", raster=F) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))


## inflammation
MAC_Inflam <- list(InflamMAC_genes)
MACInflam_object <- AddModuleScore(object = test.scvi.integrated, features = MAC_Inflam, 
                                   name = "MACInflam_marker_score")
plt_Inflam <- FeaturePlot(object = MACInflam_object, 
                          features = "MACInflam_marker_score1", raster=F) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

## Joint plot
pdf("Human_Integrated_UMAP_highlight_MACFeatures_MPCs.pdf", width = 10, height = 4)
plt_2 + plt_Res
plt_2 + plt_Foam
plt_2 + plt_C3
plt_2 + plt_IFNIC
plt_2 + plt_Inflam
dev.off()


top10 <-scvi_HUMAN_MARKERS %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)

heat<- DoHeatmap(human.merge_scvi_Average, features = top10$gene, 
                 draw.lines = F,group.bar = F,size = 4, angle = 90, raster=F) + 
       NoLegend()

heat


message("+----    Select the macrophage clusters and recluster.      --------+")

FeaturePlot(object = human.scvi.integrated, features = c("CSF1R"), raster=F) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
MAC_ONLY <- subset(human.scvi.integrated, idents= c(3,5,12,18,14,15)) 
## 28798 features, 25133 cells

UMAPPlot(MAC_ONLY)
DefaultAssay(MAC_ONLY) <- "integrated.scvi"
MAC_Reclustered <- MAC_ONLY
mygenes <- rownames(GetAssay(MAC_Reclustered))

MAC_Reclustered=FindVariableFeatures(MAC_Reclustered)
MAC_Reclustered<- ScaleData(MAC_Reclustered, features= mygenes, verbose = FALSE)
MAC_Reclustered<- RunPCA(MAC_Reclustered, npcs = 30, verbose = FALSE)
MAC_Reclustered<- RunUMAP(MAC_Reclustered, reduction = "pca", dims = 1:30)

## here we use default resolution res = 0.8
MAC_Reclustered <- FindNeighbors(MAC_Reclustered, dims = 1:30)
MAC_Reclustered <- FindClusters(MAC_Reclustered, resolution = 0.8)

## Reclustered UMAP
Reclustered_UMAP=UMAPPlot(MAC_Reclustered, label = T, raster = F)
Reclustered_UMAP

Reclustered_Markers=FindAllMarkers(MAC_Reclustered,  only.pos = T)

## save objects
save(Reclustered_Markers, file="/./Reclustered_Markers.Robj")
save(MAC_Reclustered, file="./MAC_Reclustered.Robj")

message("+--- Find top markers for each cluster and prepare annotation  -----+")

top10 <- Reclustered_Markers %>% group_by(cluster) %>% top_n(n = 10, wt =avg_log2FC)
top10

## mannually annnotated the clusters and assign below

new.cluster.ids <- c("IFNIC", "Foamy", "Res", "IFNIC", "Res",  "Inflam",
                     "Inflam", "Inflam", "Res", "Res", "Inflam", "MixMac",
                     "MixMac", "Inflam",  "Res", "Foamy", "IFNIC", "MixMac",
                     "Inflam", "MixMac", "Foamy", "MixMac", "Res", "MixMac")
newReclustered_Mac <- MAC_Reclustered
names(new.cluster.ids) <- levels(newReclustered_Mac)
newReclustered_Mac <- RenameIdents(newReclustered_Mac, new.cluster.ids)


pdf("./AFig2b-Human_Macrophage_reclustered_UMAP_all_Feb_2023.pdf",
    width=8, height = 4, onefile = T)
umap_all  <- UMAPPlot(newReclustered_Mac, label=T) + 
  scale_color_manual(values=c("darkred",  "darkgray",  "darkorchid","cyan")) +
  theme(legend.position = "none")
umap_all1 <- UMAPPlot(MAC_Reclustered, label=T) + theme(legend.position = "none")
plot_grid(umap_all1, umap_all, ncol = 2)
dev.off()

pdf("./AFig2c-Human_Macrophage_reclustered_UMAP_bystudy_Feb_2023.pdf", width = 9, 
    height = 9, onefile = T)
umap_split <- UMAPPlot(newReclustered_Mac, label=F, split.by= "Protocol", ncol = 3) + 
  scale_color_manual(values=c("darkred",  "darkgray",  "darkorchid","cyan")) +
  theme(legend.position = "none")
umap_split
dev.off()


DefaultAssay(newReclustered_Mac) <- "RNA"
pdf("Human_integrated_selMac_dotplot_selMarkers_01_02_2023_N1.pdf", 
    width=10, height=3.5)
plt_dot<- DotPlot(newReclustered_Mac, features=MAC_genes, 
                  idents=c("Res", "IFNIC", "Inflam", "Foamy")) + 
          RotatedAxis() + scale_color_viridis() +
          scale_size(range = c(2, 7)) +
          theme(legend.position="bottom", legend.box = "horizontal",
                legend.key.width= unit(2, 'cm'))
plt_dot
dev.off()

message("+--Fig 4e: Merge Inflam/IFNIC to produce selected Markers dotplot---+")
## based on the previous dot plot, order selected genes
newReclustered_Mac@meta.data$FineCluster <- Idents(newReclustered_Mac)
newReclustered_Mac@meta.data$FineCluster <- factor(newReclustered_Mac@meta.data$FineCluster, 
                                                  levels = c("Res", "Foamy", "IFNIC", "Inflam", "MixMac"))
newReclustered_Mac[["old.ident"]] <- Idents(object =newReclustered_Mac)
new.cluster.ids <- c("IFN/Inflam", "Foamy", "Res", "IFN/Inflam", "MixMac")
names(new.cluster.ids) <- levels(newReclustered_Mac)

newReclustered_Mac <- RenameIdents(newReclustered_Mac, new.cluster.ids)
newReclustered_Mac[["new.ident"]] <- Idents(object =newReclustered_Mac)
DefaultAssay(newReclustered_Mac) <- "RNA"

MAC_genes_ord <- c("CD163", "NRP1", "LYVE1", "TIMD4", "FOLR2", "LGMN", "F13A1",   
                   "LGALS3", "CD9", "MARCO", "FABP4", "CTSD", "SPP1", "GPNMB", 
                   "APOC1", "TREM2", "IFI6", "APOE", "IL1B", "CXCL8",  "IER3", 
                   "FCGR3A", "CXCL2", "CCL4", "XIST", "IFI44L", "HLA-DRA", "CD74", 
                   "EGR1", "CLEC4E" )

pdf("./Fig4e-Human_integrated_selMac_dotplot_selMarkers_Apr_2024_1.pdf", width=12, height=4)

plt_dot<- DotPlot(newReclustered_Mac, features=MAC_genes_ord, scale = T, 
                  idents=c("Res",  "Foamy", "IFN/Inflam"))+
  RotatedAxis()+scale_color_viridis() +
  scale_size(range = c(2, 7)) + xlab("") + ylab("")  +
  theme(legend.position="top",legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=16, face = "bold"), #change legend title font size
        legend.text = element_text(size=16, face = "bold"),
        axis.text.x = element_text(size=20, face = "italic"),
        axis.text.y = element_text(size=20,face = "bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.border = element_rect(colour = "black", fill=NA, size = 1))
plt_dot
dev.off()


message("+-----        NRP1 high/low list for EnrichR GO analysis       -----+")

library(enrichR)

dbs <- c("GO_Biological_Process_2021")
NRP1_markers <- FindMarkers(newReclustered_Mac, ident.1= "Res", 
                            ident.2=c("IFN/Inflam", "Foamy"), only.pos = T)
NRP1_GO <- enrichr(rownames(NRP1_markers), dbs)
NRP1_GO <- NRP1_GO[["GO_Biological_Process_2021"]]
NRP1_GO <- subset(NRP1_GO, Adjusted.P.value < 0.05)
write.xlsx(NRP1_GO, file = "Macphage_subGO_BiologicalProcess_summary.xlsx", 
           append = F, sheetName = "NRP1_highvslow GO")

message("+-------    Fig 4f: selGO NRP1 highvslow              --------------+")

selBPs <- c("receptor-mediated endocytosis (GO:0006898)", 
            "complement activation, classical pathway (GO:0006958)",
            "negative regulation of leukocyte apoptotic process (GO:2000107)",
            "transmembrane receptor protein tyrosine kinase signaling pathway (GO:0007169)",
            "positive regulation of leukocyte chemotaxis (GO:0002690)", 
            "regulation of cell migration (GO:0030334)", 
            "regulation of actin cytoskeleton reorganization (GO:2000249)",
            "axon guidance (GO:0007411)")
selBPs.1 <- c("receptor-mediated endocytosis", 
              "complement activation, \nclassical pathway",
              "negative regulation of \nleukocyte apoptotic process",
              "transmembrane receptor \nprotein tyrosine kinase \nsignaling pathway",
              "positive regulation of \nleukocyte chemotaxis", 
              "regulation of cell migration", 
              "regulation of actin \ncytoskeleton reorganization",
              "axon guidance")

BPplt_dat <- NRP1_GO[NRP1_GO$Term%in%selBPs,]
BPplt_dat$Description <-selBPs.1
BPplt_dat$Count <- as.numeric(unlist(lapply(BPplt_dat$Overlap, 
                                            function(x) strsplit(x, split="/")[[1]][1])))


pdf("./Fig4f-SelBP_scRNAHuman_NRP1_highvslow_Apr_2024-rm3.pdf", 
    width= 13, height = 7, onefile = T)
BPplt <- ggplot(BPplt_dat, aes(x=reorder(Description, -Adjusted.P.value),
                               y=Count, color=Adjusted.P.value, fill=Adjusted.P.value)) + 
  geom_bar(stat="identity", aes(alpha = -Adjusted.P.value), color="white")+
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

BPplt
dev.off()

##-------------      Finish revised on 29/04/2024                -------------## 
