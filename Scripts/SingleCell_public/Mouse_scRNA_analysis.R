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

setwd('/Users/xz289/Documents/Minoru/Seurat_v3/Athero_mouse/')

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(viridis)
  library(RColorBrewer)
  library(dplyr)
  library(cowplot)
  library(tidyr)
  library(tidyverse)
  library(xlsx)
  library(readxl)
  library(patchwork)
})


load(file="./Cochain_Chow.Robj")
load(file="./Cochain_11weeks.Robj")
load(file="./Cochain_20weeks.Robj")
load(file="./Kim_Ldlr.Robj") 
load(file="./Cochain_A11A12.Robj") ## Vafadarnejad et al. (2020) 
load(file="./Williams_21days_HFD.Robj")
load(file="./Williams_C57.Robj")


message("+-------------------------------------------------------------------+")
message("+---------            merge seven data sets.               ---------+")
message("+-------------------------------------------------------------------+")

cellids <- c("Cochain_Ldlr_11_weeks_HFD", "Cochain_Ldlr_20_weeks_HFD", 
             "Cochain_A11A12", "Cochain_Chow",  "Kim_Ldlr_10_weeks_HFD", 
             "Williams_21days_HFD", "Williams_C57")

AtheroUpdate.list <- merge(x=Cochain_11weeks, y=list(Cochain_20weeks, 
                                                     Cochain_A11A12, 
                                                     Cochain_Chow,  
                                                     Kim_Ldlr, 
                                                     Williams_21days_HFD, 
                                                     Williams_C57), 
               add.cell.ids=cellids)

AtheroUpdate.list[["RNA"]] <- split(AtheroUpdate.list[["RNA"]], 
                                    f = AtheroUpdate.list$Protocol)

AtheroUpdate.list <- AtheroUpdate.list[cellids]
head(AtheroUpdate.list)

message("+--------      Perform analysis without integration     ------------+")

AtheroUpdate.list <- NormalizeData(AtheroUpdate.list, verbose = FALSE)
AtheroUpdate.list <- FindVariableFeatures(AtheroUpdate.list, 
                                               selection.method = "vst", 
                                               nfeatures = 2000, verbose = FALSE)
AtheroUpdate.list <- ScaleData(AtheroUpdate.list)
AtheroUpdate.list <- RunPCA(AtheroUpdate.list)

AtheroUpdate.list <- FindNeighbors(AtheroUpdate.list, dims = 1:30, 
                                   reduction = "pca")


message("+Using clustree to decide the clusters with different resolutions.+")

AtheroUpdate.rlist <- AtheroUpdate.list
AtheroUpdate.rlist <- FindNeighbors(AtheroUpdate.rlist, dims = 1:30)
AtheroUpdate.rlist <- FindClusters(AtheroUpdate.rlist, resolution = 0)

for(res in c(0.3,0.4,0.5,0.6,0.7,0.8,0.9,1))
{
  message(paste0("+--Method 0 : Finding Clusters at ", res, " resolution ---+"))
  AtheroUpdate.rlist <- FindClusters(AtheroUpdate.rlist, 
                                     resolution = res, 
                                     print.output = FALSE) 
}

library(clustree)
mat.convert  <- AtheroUpdate.rlist@meta.data
mat.convert.numeric <- apply(mat.convert[,grep("integrated_snn_res.", colnames(mat.convert))], 
                             2, as.numeric)
clust.input <- as.data.frame(mat.convert.numeric)

pdf("./Athero_mouse_LDLR_new/Clustree_genecut200_cellcut3_res0_1.pdf", width = 15, height = 20)
clustree(clust.input, prefix ="integrated_snn_res.")
dev.off()

message("+-------- Final decision of res is going to use 0.3  ---------------+")

AtheroUpdate.list <- FindClusters(AtheroUpdate.list, resolution = 0.3, 
                                  cluster.name = "unintegrated_clusters")  

AtheroUpdate.list <- RunUMAP(AtheroUpdate.list, dims = 1:30, 
                             reduction = "pca", 
                             reduction.name = "umap.unintegrated")

DimPlot(AtheroUpdate.list, reduction = "umap.unintegrated", 
        group.by = c("Protocol", "seurat_clusters"))


message("+-------------------------------------------------------------------+")
message("+---------     Perform integration analysis                ---------+")
message("+-------------------------------------------------------------------+") 

AtheroUpdate.integrated <- IntegrateLayers(object = AtheroUpdate.list, 
                                           method = scVIIntegration, 
                                           orig.reduction = "pca", 
                                           new.reduction = "integrated.scvi",
                                           verbose = FALSE)

# re-join layers after integration
AtheroUpdate.integrated[["RNA"]] <- JoinLayers(AtheroUpdate.integrated[["RNA"]])

AtheroUpdate.integrated <- FindNeighbors(AtheroUpdate.integrated, 
                                         reduction = "integrated.scvi", 
                                         dims = 1:30)
AtheroUpdate.integrated <- FindClusters(AtheroUpdate.integrated, 
                                        resolution = 0.3,
                                        cluster.name = "scvi_clusters")

AtheroUpdate.integrated <- RunUMAP(AtheroUpdate.integrated, 
                                   reduction = "integrated.scvi", 
                                   dims = 1:30, reduction.name = "umap.scvi")
p1 <- DimPlot(
  AtheroUpdate.integrated,
  reduction = "umap.scvi",
  group.by = c("Protocol", "scvi_clusters"),
  combine = FALSE, label.size = 2
)

save(AtheroUpdate.integrated, 
     file="./Athero_mouse_LDLR_new/AtheroUpdate_integrated_scvi_res03.Robj")

message("+-------------------------------------------------------------------+")
message("+---------     Perform annotation analysis                 ---------+")
message("+-------------------------------------------------------------------+")


message("+--------  Find Markers at res0.3 clustering,             ----------+") 

AtheroUpdate.integrated.markers<- FindAllMarkers(AtheroUpdate.integrated, 
                                                 min.pct = 0.3, 
                                                 logfc.threshold = 0.2, 
                                                 max.cells.per.ident = 250, 
                                                 only.pos = T)

AtheroUpdate_scvi_Average <- AverageExpression(AtheroUpdate.integrated, return.seurat = T)

write.csv(AtheroUpdate.integrated.markers,
          file="./Athero_mouse_LDLR_new/AtheroUpdate.integrated.markers.res03.csv")
save(AtheroUpdate_scvi_Average,file="./Athero_mouse_LDLR_new/AtheroUpdate_scvi_Average_res03.Robj")

## select top markers generate heatmap

Markers <- AtheroUpdate.integrated.markers %>% filter(pct.1 > 0.3)
top5= AtheroUpdate.integrated.markers%>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
heat<- DoHeatmap(AtheroUpdate_scvi_Average, features = top5$gene, group.bar = T, 
                 size = 4, angle = 90,draw.lines=F, raster = F) + NoLegend()

heat
message("+--------    Add cell cycle score and IEG score           ----------+")

## load cell cycle genes, convert to mouse gene names.  
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

s.genes <-tolower(s.genes)
g2m.genes <- tolower(g2m.genes)

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

s.genes <- sapply(s.genes, simpleCap)
g2m.genes <-sapply(g2m.genes, simpleCap)

AtheroUpdate.integrated <- CellCycleScoring(AtheroUpdate.integrated, 
                                            s.features = s.genes, 
                                            g2m.features = g2m.genes, 
                                            set.ident = TRUE)

UMAPPlot(AtheroUpdate.integrated, label=T)

## UMAP for cell cycle score.

g2mphase <- FeaturePlot(AtheroUpdate.integrated, features=c("G2M.Score")) + 
  scale_color_viridis(option="inferno")

sphase <-FeaturePlot(AtheroUpdate.integrated, features=c("S.Score")) + 
  scale_color_viridis(option="inferno")

plot_grid(sphase, g2mphase)

##----- cell cycle did not make effects  --------##

## Check dissociation induced immediate early genes**

ieg.genes <-list( c("Jun", "Junb", "Jund", "Fos", "Fosb", "Atf3", "Hsp90aa1", 
                    "Hsp90ab1","Egr1", "Hspa1a", "Hspa1b", "Hspa8", "Zfp36",
                    "Cebpb", "Cebpd", "Socs3", "Hspe1", "Dusp1"))

## Add the IEG score to the metadata.

AtheroUpdate.integrated <- AddModuleScore(AtheroUpdate.integrated, 
                                          features = ieg.genes, name="IEG_Score", ctrl = 10)


head(AtheroUpdate.integrated@meta.data)
FeaturePlot(AtheroUpdate.integrated, features=c("IEG_Score1"), 
            min.cutoff = "q5", 
            max.cutoff = "q95") + 
  scale_color_viridis(option="inferno") + 
  labs(title="Immediate Early Gene expression score")

message("+--------      Identify macrophage clusters                ---------+")

head(AtheroUpdate.integrated@meta.data)
Idents(AtheroUpdate.integrated)=AtheroUpdate.integrated[["scvi_clusters"]]
levels(AtheroUpdate.integrated)

save(AtheroUpdate.integrated,file="./Athero_mouse_LDLR_new/AtheroUpdate.integrated.res03.Robj")

## apply a group of macrophage surface markers to identify the macrophage cluster.

Mac_genes <- c("Sepp1", "Pf4", "F13a1", "Folr2", "Nrp1", "Lyve1", "Cd163", "Ednrb",
               "Timd4", "Phactr1", "Egr1"," Lgals3", "Ctsd", "Cxcl2", "Ccl4",
               "Trem2", "Spp1", "Cd9", "Gpnmb", "Adgre1", "Il1b", "Cd74", "Stat1", 
               "Ifit3", "Mnda", "Irf7", "Isg15", "Ifit2" , "Fcgr4", "Fcgr1")

pdf("./Athero_mouse_LDLR_new/AFig1a-Mouse_integrated_all_UMAP_MacGenes.pdf", 
    width = 9, height = 4, onefile = T)
umap_inte <- UMAPPlot(AtheroUpdate.integrated, label=T) + theme(legend.position = "none")
mac_int   <- list(Mac_genes)
object    <- AtheroUpdate.integrated
DefaultAssay(object) <- "RNA"
object    <- AddModuleScore(object = object, features = mac_int, name = "Mac_Genes")
macF      <- FeaturePlot(object = object, features = "Mac_Genes1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
plot_grid(umap_inte, macF, ncol = 2)
dev.off()

## Key macrophages featurePlot to identify the Macrophage clusters.

pdf("./Athero_mouse_LDLR_new/Mouse_integrated_FeaturePlot_Adgre1_Fcgr1_Csf1r_Pf4.pdf", onefile = T)

macAdgre1 <- FeaturePlot(object = AtheroUpdate.integrated, features = "Adgre1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
macFcgr1 <- FeaturePlot(object = AtheroUpdate.integrated, features = "Fcgr1") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
macCsf1r <- FeaturePlot(object = AtheroUpdate.integrated, features = "Csf1r") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
macPf4 <- FeaturePlot(object = AtheroUpdate.integrated, features = "Pf4") +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
plot_grid(macAdgre1, macFcgr1, macCsf1r, macPf4, ncol = 2)
dev.off()

message("+--------  Reclustering the macrophage clusters            ---------+") 

AtheroUpdate_MAC <- subset(AtheroUpdate.integrated, idents=c(0,1,2,7,3,6,14))
save(AtheroUpdate_MAC, file="./Athero_mouse_LDLR_new/AtheroUpdate.integrated.res03.addc7.Macrophage.Robj")

MAC_Reclustered <- AtheroUpdate_MAC
DefaultAssay(MAC_Reclustered) <- "integrated.scvi"
MAC_Reclustered=FindVariableFeatures(MAC_Reclustered)
mygenes <- rownames(GetAssay(MAC_Reclustered))
MAC_Reclustered <- ScaleData(MAC_Reclustered, features= mygenes, verbose = FALSE)
MAC_Reclustered <- RunPCA(MAC_Reclustered, npcs = 30, verbose = FALSE)
MAC_Reclustered <- RunUMAP(MAC_Reclustered, reduction = "pca", dims = 1:30)
MAC_Reclustered <- FindNeighbors(MAC_Reclustered, dims = 1:30)
MAC_Reclustered <- FindClusters(MAC_Reclustered, resolution = 1.2)

pdf("./Athero_mouse_LDLR_new/Macrophage_c01236714_Ldlr_ChowHFD_res12_c19_UMAP.pdf", width=10, height=4)
umap_all <- UMAPPlot(MAC_Reclustered, label=T)
macair_highlight <- rownames(MAC_Reclustered@meta.data[MAC_Reclustered@meta.data$Williams_C57_ID=="Mac_AIR",])
umap_macair<- UMAPPlot(MAC_Reclustered, label=F, cells.highlight=list(macair_highlight))
plot_grid(umap_all, umap_macair, ncol = 2)
dev.off()

MAC.integrated.markers <- FindAllMarkers(MAC_Reclustered, min.pct = 0.3, 
                                         logfc.threshold = 0.20, only.pos = T)

MAC.top20 = MAC.integrated.markers%>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

## add meta information of treatment, Ctrl and HFD

MAC_Reclustered@meta.data$Treatment <- ifelse(MAC_Reclustered@meta.data$Protocol == "Cochain_Chow"|MAC_Reclustered@meta.data$Protocol == "Williams_C57", 
                                              "Ctrl", "HFD")

## find conserve markers for some clusters between treatment

MAC.conserve.markers9 <- FindConservedMarkers(MAC_Reclustered,
                                              ident.1 = 9,
                                              grouping.var = "Treatment",
                                              only.pos = TRUE)
MAC.conserve.markers7 <- FindConservedMarkers(MAC_Reclustered,
                                              ident.1 = 7,
                                              grouping.var = "Treatment",
                                              only.pos = TRUE)
MAC.conserve.markers13 <- FindConservedMarkers(MAC_Reclustered,
                                               ident.1 = 13,
                                               grouping.var = "Treatment",
                                               only.pos = TRUE)
MAC.conserve.markers10 <- FindConservedMarkers(MAC_Reclustered,
                                               ident.1 = 10,
                                               grouping.var = "Treatment",
                                               only.pos = TRUE)
MAC.conserve.markers1 <- FindConservedMarkers(MAC_Reclustered,
                                              ident.1 = 1,
                                              grouping.var = "Treatment",
                                              only.pos = TRUE)
MAC.conserve.markers12 <- FindConservedMarkers(MAC_Reclustered,
                                               ident.1 = 12,
                                               grouping.var = "Treatment",
                                               only.pos = TRUE)
MAC.conserve.markers2 <- FindConservedMarkers(MAC_Reclustered,
                                              ident.1 = 2,
                                              grouping.var = "Treatment",
                                              only.pos = TRUE)
MAC.conserve.markers6 <- FindConservedMarkers(MAC_Reclustered,
                                              ident.1 = 6,
                                              grouping.var = "Treatment",
                                              only.pos = TRUE)
MAC.conserve.markers0 <- FindConservedMarkers(MAC_Reclustered,
                                              ident.1 = 0,
                                              grouping.var = "Treatment",
                                              only.pos = TRUE)
MAC.conserve.markers3 <- FindConservedMarkers(MAC_Reclustered,
                                              ident.1 = 3,
                                              grouping.var = "Treatment",
                                              only.pos = TRUE)

## average expression for reclustering markers.

MAC_Reclustered_scvi_Average <- AverageExpression(MAC_Reclustered, return.seurat = T)
MAC.top5 = MAC.integrated.markers%>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
heatmac <- DoHeatmap(MAC_Reclustered_scvi_Average, features = MAC.top5$gene, group.bar = T, size = 4, angle = 90,draw.lines=F, raster = F) + NoLegend()

## top markers heatmap
heatmac

message("+---   Combine SuppTable for MPC markers to assist clustering ------+")

MPCs <- read.xlsx("./cvac161_supplementary_data/Supp_xls_1_Mouse_MPC_Marker_Genes.xlsx", sheetIndex=1)
MPCs <- MPCs[-1,]
MPC_clust <- unique(MPCs[-1,7])
addM_MAC  <- MAC_Reclustered
DefaultAssay(addM_MAC) <- "RNA"

Fplt_fn <- function(seobj, Minput, clustM, clustCol, geneCol){
  celltype_marker_gene_list <- list(Minput[Minput[,clustCol]==clustM, geneCol])
  newobject <- AddModuleScore(object = seobj, features = celltype_marker_gene_list, 
                              name = paste0(clustM, "_score"))
  Fplt_M    <- FeaturePlot(object = newobject, features = paste0(clustM, "_score1")) +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
    
  Fplt_M
}
MclustPlt <- list()
for(i in c(1:6,8:15)){
  print(i)
  MclustPlt[[i]] <- Fplt_fn(addM_MAC, MPCs, MPC_clust[i], 7, 8)
  MclustPlt
}

i = 7
celltype_marker_gene_list <- list(MPCs[MPCs[,7]=="CCR2intMHCII+", 8])
newobject <- AddModuleScore(object = addM_MAC, features = celltype_marker_gene_list, 
                              name = "CCR2intMHCII_score")
Fplt_M    <- FeaturePlot(object = newobject, features = "CCR2intMHCII_score1") +
    scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))
MclustPlt[[7]] <-  Fplt_M  

pdf("./Athero_mouse_LDLR_new/MPCs_FeaturePlot_clustersMarkers_Res12_c19_June_2023.pdf", 
    width = 16, height = 16)
umapcplt <- DimPlot(MAC_Reclustered, label=T)
plot_grid(MclustPlt[[1]], MclustPlt[[2]], MclustPlt[[3]], MclustPlt[[4]],
          MclustPlt[[5]], MclustPlt[[6]], MclustPlt[[7]], MclustPlt[[8]],
          MclustPlt[[9]], MclustPlt[[10]], MclustPlt[[11]], MclustPlt[[12]],
          MclustPlt[[13]], MclustPlt[[14]], MclustPlt[[15]], umapcplt, ncol = 4, byrow=T)
dev.off()

message("+---- Refine the clusters annotations based on above process ------+")

new.cluster.ids <- c("Adv/Res", "Adv/Res", "CCR2-MHCII", "Adv/Res",
                     "Inflam", "Mac_AIR/Trem2", "Inflam", "Trem2-Slamf9",
                     "Adv/Res", "IFNIC", "Inflam", "Adv/Res", "Prolif_G2M",
                     "Inflam", "Adv/Res", "Prolif_G2M", "Monocyte", 
                     "Mac_AIR/Trem2", "Mono/DCs")
new.cluster.idsm <- c("Adv/Res", "Adv/Res", "CCR2-MHCII", "Adv/Res",
                     "IFN/Inflam", "Mac_AIR/Trem2", "IFN/Inflam", "Mac_AIR/Trem2",
                     "Adv/Res", "IFN/Inflam", "IFN/Inflam", "Adv/Res", "Prolif_G2M",
                     "IFN/Inflam", "Adv/Res", "Prolif_G2M", "Monocyte", 
                     "Mac_AIR/Trem2", "Mono/DCs")

message("+-------Add the cluster celltype annotations-------+")

## global annotation with 3 main subtype of Macrophages
Idents(MAC_Reclustered) <- MAC_Reclustered[["seurat_clusters"]]
names(new.cluster.idsm) <- levels(MAC_Reclustered)
MAC_Reclustered <- RenameIdents(MAC_Reclustered, new.cluster.idsm)
MAC_Reclustered[["Global_scvi_Annotation"]] <- Idents(object = MAC_Reclustered)

## Fine annotation 6 subcluster of macrophage
Idents(MAC_Reclustered) <- MAC_Reclustered[["seurat_clusters"]]
names(new.cluster.ids) <- levels(MAC_Reclustered)
MAC_Reclustered <- RenameIdents(MAC_Reclustered, new.cluster.ids)
MAC_Reclustered[["Fine_scvi_Annotation"]] <- Idents(object = MAC_Reclustered)

pdf("./Athero_mouse_LDLR_new/MouseMacrophage_reclustered_Ldlr_ChowHFD_UMAP_MacAir_Res12_C19.pdf",
    width = 12, height = 10, onefile = T)
Idents(MAC_Reclustered) <- MAC_Reclustered[["Global_scvi_Annotation"]]
umap_allG  <- UMAPPlot(MAC_Reclustered, label=T) 
Idents(MAC_Reclustered) <- MAC_Reclustered[["Fine_scvi_Annotation"]]
umap_allf   <- UMAPPlot(MAC_Reclustered, label=T)
macair_highlight <- rownames(MAC_Reclustered@meta.data[MAC_Reclustered@meta.data$Williams_C57_ID=="Mac_AIR",])
umap_macair<- UMAPPlot(MAC_Reclustered, label=F, cells.highlight=list(macair_highlight)) 
plot_grid(umap_allG, umap_macair, umap_all, umap_allf, ncol = 2)
dev.off()

## Make fine plot with consistent macrophage subset colors.

pdf("./Athero_mouse_LDLR_new/Macrophage_Ldlr_ChowHFD_res12_c19_seldefineColor_UMAP_Global_GlobalF.pdf", 
    width = 5, height = 4)
Idents(MAC_Reclustered) <- MAC_Reclustered[["Global_scvi_Annotation"]]
UMAP= UMAPPlot(MAC_Reclustered, label=F)+
  scale_color_manual(values=c("darkred", "darkorange", "darkgray",  "darkorchid",
                              "cyan", "pink", "green"))
LabelClusters(UMAP, id = "ident",  fontface = "bold", size=5)

Idents(MAC_Reclustered) <- MAC_Reclustered[["seurat_clusters"]]
## merge the MacAir and Trem2, split inflam and IFNIC
new.cluster.ids1 <- c("Adv/Res", "Adv/Res", "CCR2-MHCII", "Adv/Res",
                      "Inflam", "Mac_AIR/Trem2", "Inflam", "Mac_AIR/Trem2",
                      "Adv/Res", "IFNIC", "Inflam", "Adv/Res", "Prolif_G2M",
                      "Inflam", "Adv/Res", "Prolif_G2M", "Monocyte", 
                      "Mac_AIR/Trem2", "Mono/DCs")
names(new.cluster.ids1) <- levels(MAC_Reclustered)
MAC_Reclustered <- RenameIdents(MAC_Reclustered, new.cluster.ids1)
MAC_Reclustered[["GlobalF_scvi_Annotation"]] <- Idents(object = MAC_Reclustered)
Idents(MAC_Reclustered) <- MAC_Reclustered[["GlobalF_scvi_Annotation"]]
UMAPF= UMAPPlot(MAC_Reclustered, label=F)+
  scale_color_manual(values=c("darkred", "darkorange", "darkgray", "darkorchid",
                              "gray85", "cyan", "pink", "green"))
LabelClusters(UMAPF, id = "ident",  fontface = "bold", size=5)

dev.off()


## merge MacAIR/Trem2, IFN/Inflam UMAP
pdf("./Athero_mouse_LDLR_new/AFig1b-Macrophage_Ldlr_ChowHFD_res12_c19_seldefineColor_UMAP.pdf", width = 9, height = 3.5)
UMAPF1 = UMAPPlot(MAC_Reclustered, label=T)+
  scale_color_manual(values=c("darkred", "darkorange", "darkgray",  "darkorchid",
                              "gray85", "cyan", "pink", "green")) 

plot_grid(UMAPF1, umap_macair, ncol = 2)
dev.off()

## split by study

pdf("./Athero_mouse_LDLR_new/AFig1c_Macrophage_Ldlr_ChowHFD_res12_c19_seldefineColor_UMAP_GlobalF.pdf", width = 9, height = 9)
UMAP_study = UMAPPlot(MAC_Reclustered, label=F, split.by = "Protocol", ncol = 3)+
  scale_color_manual(values=c("darkred", "darkorange", "darkgray", "darkorchid", 
                              "gray85", "cyan","pink", "green")) + 
  theme(legend.position="none")

UMAP_study
dev.off()

message("+---  Dotplot for selectearkers of each subtype of Macrophage  -----+")
message("+--------               Fig 4d                          ------------+")

DefaultAssay(MAC_Reclustered) <- "RNA"
Mac_genes <- c("Pf4","F13a1","Lyve1","Csf1r","Cd163","Ednrb","Folr2","Mrc1","Sepp1",
               "Nrp1","Timd4","Cd209f", "Cd209d","Phactr1",
               "Cxcl2","Nfkbiz","Nfkbia","Nlrp3","Cd14","Cebpb","Il1b","Il1rn", 
               "Ifit3", "Stat1", "Mnda", "Irf7","Isg15",  "Ifit2" ,  "Ccl12", "Fcgr4",
               "Lgals3", "Trem2","Cd9","Ctsb","Gngt2","Spp1","Acp5","Fabp5","Cadm1",
               "Mmp12", "Itgax","Gpnmb")

Idents(MAC_Reclustered) <- MAC_Reclustered[["Global_scvi_Annotation"]]
MAC_Reclustered_sub <- subset(MAC_Reclustered, idents=c("Adv/Res", "IFN/Inflam", "Mac_AIR/Trem2"))
levels(MAC_Reclustered_sub) <- c("Adv/Res", "IFN/Inflam", "Mac_AIR/Trem2")

pdf("./Athero_mouse_LDLR_new/Fig4d-Mouse_integrated_selMac_dotplot_selMarkers_MerCluster_Res12_C19_MC3.pdf", width=10, height=3.5)
plt_dotG <- DotPlot(MAC_Reclustered_sub, features=Mac_genes, scale = T)+
  RotatedAxis()+scale_color_viridis() +
  scale_size(range = c(2, 7)) + xlab("") + ylab("") +
  theme(legend.position="top",legend.key.size = unit(1, 'cm'), #change legend key size
        legend.key.height = unit(0.5, 'cm'), #change legend key height
        legend.key.width = unit(1.2, 'cm'), #change legend key width
        legend.title = element_text(size=12, face = "bold"), #change legend title font size
        legend.text = element_text(size=10, face = "bold"),
        axis.text.x = element_text(size=12, face = "italic"),
        axis.text.y = element_text(size=12,face = "bold"))
plt_dotG
dev.off()

message("+Generate the key markers list for each subtype of Macrophage and overlap with iWDvscWD sig DEGs+")

MacList <- FindAllMarkers(MAC_Reclustered_sub, only.pos = T)
save(MacList, file = "./Athero_mouse_LDLR_new/Allmarkers_Mactypes_list.RData")
## Supplementary Table S3 with overlap list.


message("+different subtypes of Macrophage under different time DEGs analysis+")

##Prepare the metadata with different conditions as compared
## 1) Ctrl HFD, Condition
testall <- MAC_Reclustered
Idents(testall) <- testall@meta.data$Protocol
testall@meta.data$Condition <- testall@meta.data$Treatment

## 2) assign different Conditions with cell type specific
## global condition
testall@meta.data$GlobalCondition <- paste0(testall@meta.data$Global_scvi_Annotation, "_",
                                            testall@meta.data$Condition)
table(testall@meta.data$GlobalCondition)

## 3) split by time (3wk, 10+11wks, 20 week for HFD) and add celltype
timeCond <- ifelse(testall@meta.data$Protocol=="Cochain_Ldlr_20_weeks_HFD", "wk20_HFD",
                   "wk11_HFD")
timeCond <- ifelse(testall@meta.data$Protocol=="Williams_C57" | testall@meta.data$Protocol=="Cochain_Chow", "Ctrl", timeCond)
timeCond <- ifelse(testall@meta.data$Protocol=="Williams_21days_HFD", "wk3_HFD", timeCond)
table(timeCond)
testall@meta.data$TimeCond <- timeCond

testall@meta.data$NewTimeCond <- paste0(testall@meta.data$Global_scvi_Annotation, "_",timeCond)
table(testall@meta.data$NewTimeCond)

MAC_Reclustered <- testall

save(MAC_Reclustered, 
     file="./Athero_mouse_LDLR_new/NewMAC_Reclustered_withConditions_meta_Robj_Res12_C19_GC3_GFC4_June_2023.RData")

message("+-----      Nrp1 highvslow gene list identification         --------+")

Nrp1_markers <- FindMarkers(MAC_Reclustered, ident.1= "Adv/Res", 
                            ident.2=c("IFN/Inflam", "Mac_AIR/Trem2"), 
                            only.pos = T)
## gene ontology applied to enrichr due to small number of genes identified.
## results see Supplementary Table S10
## convert mouse genes to human genes to apply enrichr input.

library(Homo.sapiens)
library(TxDb.Hsapiens.UCSC.hg38.refGene)
TxDb(Homo.sapiens) <- TxDb.Hsapiens.UCSC.hg38.refGene
library(Mus.musculus)
library(Orthology.eg.db)

## cob together a function
mapFun <- function(genes) {
  mouseg <- mapIds(org.Mm.eg.db, genes, "ENTREZID","SYMBOL")
  humeg <- mapIds(Orthology.eg.db, mouseg, "Homo.sapiens", "Mus.musculus")
  humdat <- select(Homo.sapiens, humeg, c("SYMBOL","CDSCHROM"), "GENEID")
  mousedat <- select(Mus.musculus, mouseg,  c("SYMBOL","CDSCHROM"), "GENEID")
  return(list(humdat = humdat, mousedat = mousedat))
}

Ortho_Res  <- lapply(Nrp1_markers, function(x) mapFun(x))
Ortho_Resh <- lapply(Ortho_Res, function(x) x$humdat)

Nrp1_hmarkers <- unique(unlist(lapply(Ortho_Resh, function(x) x$SYMBOL)))
dbs <- c("GO_Biological_Process_2021")

Nrp1_GO <- enrichr(rownames(Nrp1_hmarkers), dbs)
Nrp1_GO <- Nrp1_GO[["GO_Biological_Process_2021"]]
Nrp1_GO <- subset(Nrp1_GO, Adjusted.P.value < 0.05)
write.xlsx(Nrp1_GO, file = "Mouse-Macphage_subGO_BiologicalProcess_summary.xlsx", 
           append = F, sheetName = "NRP1_highvslow GO")

message("+----    Identify DEGs for the following comparisons interest ------+")

## HFDvsCtrl Adv/Res, IFN, Inflam, Mac_AIR/Trem2;
## HFDvsCtrl across time points 
## From the paper "Bias, robustness and scalability in single-cell differential expression analysis" 
## by Soneson and Robinson (2018, https://www.nature.com/articles/nmeth.4612, 
## DOI:: http://dx.doi.org/10.1038/nmeth.4612), I will use Figure 5 summary, 
## choose the "MAST" method which shows the better performance across different method and with flexible modelling. 

DefaultAssay(MAC_Reclustered) <- "RNA"
Idents(MAC_Reclustered) <- MAC_Reclustered[["GlobalCondition"]]

AdvR_MAST   <- FindMarkers(MAC_Reclustered, assay = "RNA", ident.1 = "Adv/Res_HFD", 
                           ident.2 = "Adv/Res_Chow", verbose = FALSE, 
                           logfc.threshold = 0.00, test.use = "MAST")

IFNI_MAST   <- FindMarkers(MAC_Reclustered, assay = "RNA", ident.1 = "IFN/Inflam_HFD", 
                           ident.2 = "IFN/Inflam_Chow", verbose = FALSE, 
                           logfc.threshold = 0.00, test.use = "MAST")

AIRT_MAST   <- FindMarkers(MAC_Reclustered, assay = "RNA", ident.1 = "Mac_AIR/Trem2_HFD", 
                           ident.2 ="Mac_AIR/Trem2_Chow", verbose = FALSE, 
                           logfc.threshold = 0.00, test.use = "MAST")

write.xlsx(AdvR_MAST, file = "./Athero_mouse_LDLR_new/DEGs_MACs_HFDvsCtrl_MAST_summaryTable_June_2023.xlsx",
           append = F, sheetName = "Res_MAST")
write.xlsx(IFNI_MAST, file = "./Athero_mouse_LDLR_new/DEGs_MACs_HFDvsCtrl_MAST_summaryTable_June_2023.xlsx", 
           append = T, sheetName = "IFN_Inflam_MAST")
write.xlsx(AIRT_MAST, file = "./Athero_mouse_LDLR_new/DDEGs_MACs_HFDvsCtrl_MAST_summaryTable_June_2023.xlsx", 
           append = T, sheetName = "MacAIR_Trem2_MAST")

## Perform the DEGs analysis across the time for global annotation 
## with 3 clusters for Macrophages.
## 10/11 weeks will put together. Named as wk11

Idents(MAC_Reclustered) <- MAC_Reclustered[["NewTimeCond"]]
Timelab  <- c("wk3", "wk11", "wk20")

AdvR_labs <- cbind(paste0("Adv/Res_",Timelab,"_HFD"), 
                   rep("Adv/Res_Ctrl",3))
IFNI_labs <- cbind(paste0("IFN/Inflam_",Timelab,"_HFD"),
                   rep("IFN/Inflam_Ctrl",3))
AIRT_labs <- cbind(paste0("Mac_AIR/Trem2_",Timelab,"_HFD"), 
                   rep("Mac_AIR/Trem2_Ctrl",3))

AdvR_Time_MAST   <- apply(AdvR_labs, 1, 
                          function(x) FindMarkers(MAC_Reclustered, assay = "RNA", 
                                                  ident.1 = x[1], ident.2 = x[2],
                                                  verbose = FALSE, 
                                                  logfc.threshold = 0.00,
                                                  test.use = "MAST"))
IFNI_Time_MAST   <- apply(IFNI_labs, 1, 
                          function(x) FindMarkers(MAC_Reclustered, assay = "RNA", 
                                                  ident.1 = x[1], ident.2 = x[2],
                                                  verbose = FALSE, 
                                                  logfc.threshold = 0.00,
                                                  test.use = "MAST"))

AIRT_Time_MAST   <- apply(AIRT_labs, 1, 
                          function(x) FindMarkers(MAC_Reclustered, assay = "RNA", 
                                                  ident.1 = x[1], ident.2 = x[2],
                                                  verbose = FALSE, 
                                                  logfc.threshold = 0.00,
                                                  test.use = "MAST"))

save(AdvR_MAST, AdvR_Time_MAST, AIRT_MAST, AIRT_Time_MAST, IFNI_MAST, IFNI_Time_MAST,
    file = "./Athero_mouse_LDLR_new/DEGs_methodCompare_MACs_ResAdv_HFDvsCtrl_MAST_summaryTable_June_2023.RData")


message("+- genes expression changes across time between treatment in different subtype Mac----+")

## heatmap for selected MacGenes for three macrophage types with log2FoldChange 
## across time wk3, wk11 and wk20

macheatGenes <- c(Mac_genes, "Igfbp4", "Igf1", "Igf1r")
Res_Ghmat1 <- Res_MAST[rownames(Res_MAST)%in%macheatGenes, c("avg_log2FC", "p_val_adj")]
Res_Ghmat1$gene <- rownames(Res_Ghmat1)
IFN_Ghmat1 <- IFNI_MAST[rownames(IFNI_MAST)%in%macheatGenes, c("avg_log2FC", "p_val_adj")]
IFN_Ghmat1$gene <- rownames(IFN_Ghmat1)
AIR_Ghmat1 <- AIRT_MAST[rownames(AIRT_MAST)%in%macheatGenes, c("avg_log2FC", "p_val_adj")]
AIR_Ghmat1$gene <- rownames(AIR_Ghmat1)
macheat_dat <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all=T),
                      list(Res_Ghmat1, IFN_Ghmat1, AIR_Ghmat1))

macheat_pltd <- macheat_dat[, c(2,4,6)]
rownames(macheat_pltd) <- macheat_dat$gene
colnames(macheat_pltd) <- c("Adv_Res", "IFN_Inflam", "MacAIR_Trem2")


## prepare heatmap data 
library(circlize)
library(grid)

## heatmap datasort function
hpltm_datasort <- function(datinput, mgenes, selcolName){
  subd <- datinput[rownames(datinput)%in%mgenes, selcolName]
  subd$gene <- rownames(subd)
  subd
}

## selected key markers
Res_markers <- c("Pf4","F13a1","Lyve1","Csf1r","Cd163","Ednrb","Folr2","Mrc1","Sepp1",
                 "Nrp1","Timd4","Cd209f", "Cd209d", "Igfbp4", "Igf1", "Igf1r")
IFN_markers <- c("Cxcl2","Nfkbiz","Nfkbia","Nlrp3","Cd14","Cebpb","Il1b","Il1rn", 
                 "Stat1", "Mnda", "Irf7","Isg15",  "Ccl12", "Fcgr4")
AIR_markers <- c("Lgals3", "Trem2","Cd9","Ctsb","Spp1","Fabp5","Slamf9", "Maf", "Acp5",
                 "Gngt2", "Mpeg1", "Cd72", "Gpnmb", "Ctsd", "H2-Eb1", "H2-Aa", "Cadm1", "Mmp12")

selcolName <- c("avg_log2FC", "p_val_adj")

Res_selhmatT <- lapply(AdvR_Time_MAST,  function(x) hpltm_datasort(x, Res_markers, selcolName))
IFN_selhmatT <- lapply(IFNI_Time_MAST,  function(x) hpltm_datasort(x, IFN_markers, selcolName))
AIR_selhmatT <- lapply(AIRT_Time_MAST,  function(x) hpltm_datasort(x, AIR_markers, selcolName))

Resheat_seldatT <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all=T), Res_selhmatT)[,c(1,2,4,6)]
IFNheat_seldatT <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all=T), IFN_selhmatT)[,c(1,2,4,6)]
AIRheat_seldatT <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all=T), AIR_selhmatT)[,c(1,2,4,6)]

colnames(Resheat_seldatT)[2:4] <- paste0(c(3,11,20), "wks") 
colnames(IFNheat_seldatT)[2:4] <- paste0(c(3,11,20), "wks")
colnames(AIRheat_seldatT)[2:4] <- paste0(c(3,11,20), "wks")

Resheat_seldatT1 <- Resheat_seldatT[,-1]; rownames(Resheat_seldatT1) <- Resheat_seldatT$gene
IFNheat_seldatT1 <- IFNheat_seldatT[,-1]; rownames(IFNheat_seldatT1) <- IFNheat_seldatT$gene
AIRheat_seldatT1 <- AIRheat_seldatT[,-1]; rownames(AIRheat_seldatT1) <- AIRheat_seldatT$gene

all_seldatT1  <- rbind(Resheat_seldatT1, IFNheat_seldatT1, AIRheat_seldatT1)
all_seldatT1$Time <- rep(c("Res", "IFN", "AIR"), c(dim(Resheat_seldatT1)[1],
                                                   dim(IFNheat_seldatT1)[1],
                                                   dim(AIRheat_seldatT1)[1]))



## merge all data and sort by sub macrophage type. 

allhpt_seldat <- allhpt_seldat1[order(allhpt_seldat1$Time),]
allhpt_selplt <- allhpt_seldat[,-c(1,5,7)]
colnames(allhpt_selplt) <- c("wkA", "wkB", "wkC")
rownames(allhpt_selplt) <- allhpt_seldat[,1] 

message("+-------                        Extended Fig 7a           ----------+")

colOrd <- factor(colnames(Resheat_seldatT1), levels = c("3wks",  "11wks", "20wks"))
hplt_res_all <- ComplexHeatmap::Heatmap(as.matrix(allhpt_selplt), na_col ="black",
                                       name = "macheatTimeres",  
                                       row_title = "", column_title = NULL, 
                                       show_row_names = TRUE, show_column_names = TRUE,
                                       heatmap_legend_param = list(title = "log2FoldChange", 
                                                                   legend_height = unit(6, "cm"), 
                                                                   title_position = "leftcenter-rot"), 
                                       cluster_columns = FALSE, row_names_side ="right", 
                                       cluster_rows = TRUE, column_split = c("wkA", "wkB", "wkC"),
                                       
                                       row_split = rep(c("IFN/Inflam", "Mac_AIR/Trem2","Adv/Res"), 
                                                       c(dim(AIRheat_seldatT1)[1],dim(IFNheat_seldatT1)[1],
                                                         dim(Resheat_seldatT1)[1])),
                                       cluster_row_slices = T,
                                       column_labels = c("3 weeks", "11 weeks", "20 weeks"),
                                       column_order = c(1,2,3,4),
                                       column_names_rot = 0, column_names_side = "bottom", 
                                       column_title_rot = 0, column_names_centered = TRUE,
                                       column_names_gp = gpar(fontsize = 10, fontface = "bold"),
                                       row_names_gp = gpar(fontsize = 10, fontface = "bold.italic"),
                                       column_title_gp = gpar(fontsize = 12, fontface = "bold"),
                                       row_title_gp = gpar(fontsize = 12, fontface = "bold"),
                                       row_title_rot = 90, row_gap = unit(3, "mm"), 
                                       column_gap = unit(2, "mm"))


pdf("./Athero_mouse_LDLR_new/ExtFig7a-Heatmap_DEGs_MAST_GlobalC3_acrossTime_selMarkers_July_2023.pdf",
    height = 7, width = 5)
  
hplt_res_all

dev.off()

message("+-------                        Extended Fig 7b           ----------+")
## simply just plot the proportion of each subtype of macrophage following time
## Res, IFN/Inflam, MacAIR_Trem2 at wk3/10-11/20 HFD

load("Athero_mouse_LDLR_new/NewMAC_Reclustered_withConditions_meta_Robj_Res12_C19_GC3_GFC4_June_2023.RData")
newmeta <- MAC_Reclustered@meta.data
newTab  <- as.data.frame(table(newmeta$NewTimeCond))
newTab_HFD <- newTab[grep("HFD", newTab$Var1),]
newTab_HFDmac <- newTab_HFD[c(1:3,7:12),]
newTab_HFDmac$CellType <- rep(c("Res", "IFN/Inflam", "Mac_AIR/Trem2"), each = 3)
newTab_HFDmac$Time <- rep(c("wk11", "wk20", "wk3"), length=9)
newTab_HFDmac$Treatment <- paste0(newTab_HFDmac$Time, "_HFD")
rownames(newTab_HFDmac) <- newTab_HFDmac[,1]
newTab_HFDmac <- newTab_HFDmac[,-1]
newTab_HFDmac <- newTab_HFDmac[,c("CellType", "Time", "Treatment", "Freq")]
newTab_HFDmac$Time <- factor(newTab_HFDmac$Time, 
                             levels = c("wk20", "wk11", "wk3"))
newTab_HFDmac$CellType <- factor(newTab_HFDmac$CellType, 
                                 levels = c("IFN/Inflam", "Mac_AIR/Trem2", "Res"))

library(dplyr)
newTab_HFDmac_new <- newTab_HFDmac %>%
                     group_by(Time) %>%
                     arrange(Treatment, desc(Time)) 

Proportion <- c(newTab_HFDmac_new[1:3, "Freq"]/sum(newTab_HFDmac_new[1:3, "Freq"])*100,
                newTab_HFDmac_new[4:6, "Freq"]/sum(newTab_HFDmac_new[4:6, "Freq"])*100,
                newTab_HFDmac_new[7:9, "Freq"]/sum(newTab_HFDmac_new[7:9, "Freq"])*100)
newTab_HFDmac_new$Proportion <- unlist(Proportion) 
pdf("./Athero_mouse_LDLR_new/ExtFig7b-Barplot_DEGs_MAST_GlobalC3_acrossTime_selMarkers_July_2023.pdf")
p <- ggplot(data = newTab_HFDmac_new, aes(x = Time, y = Proportion)) +
  geom_col(aes(fill = CellType), width = 0.7) + coord_flip() +
  scale_fill_manual(values = c( "darkgray", "darkorchid","darkred"))
p
dev.off() 

##--------------------------- Finish. revised on 25/04/2024---------------------##
