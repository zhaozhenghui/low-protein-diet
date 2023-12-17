

rm(list = ls())

library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(patchwork)

M2_2<-Seurat::Load10X_Spatial("D:\\spatial\\spatial\\1M")


plot1 <- VlnPlot(M2_2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(M2_2, features = "nCount_Spatial", pt.size.factor = 6) + theme(legend.position = "right")
plot_grid(plot1, plot2)



#2M_2归一化

ovary_2M_2 <- SCTransform(M2_2, assay = "Spatial", return.only.var.genes = FALSE, verbose = FALSE)

c <- SpatialFeaturePlot(ovary_2M_2, features = c("Inhbb", "Grb14","Ube2c","Prss23"),ncol = 2,pt.size.factor = 5, alpha = c(0.1, 1))
ggsave("PD14-2-Grb14ming.tiff",plot = c, device = "tiff", width = 10, height = 10, dpi = 300)

SpatialFeaturePlot(ovary_2M_2, features = c("Onecut2", "Akr1c18"),pt.size.factor = 6)
SpatialFeaturePlot(ovary_2M_2, features = c("Prss23", "Grb14","Ube2c","Top2a"),ncol = 2,pt.size.factor = 6)
ggsave("1M-Grb14.tiff",plot = c, device = "tiff", width = 10, height = 10, dpi = 300)


DefaultAssay(ovary_2M_2) <- "SCT"
VariableFeatures(ovary_2M_2) <- VariableFeatures(ovary_2M_2)



#降维、聚类、可视化
ovary_2M_2 <- RunPCA(ovary_2M_2, assay = "SCT", verbose = FALSE)
ovary_2M_2 <- FindNeighbors(ovary_2M_2, reduction = "pca", dims = 1:30)
ovary_2M_2 <- FindClusters(ovary_2M_2, verbose = FALSE)
ovary_2M_2 <- RunUMAP(ovary_2M_2, reduction = "pca",dims = 1:30)

p1 <- DimPlot(ovary_2M_2, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(ovary_2M_2, label = TRUE, label.size = 3, pt.size.factor =6)
plot_grid(p1, p2)


#单细胞整合
load("D:/mousemodel/LPD/LPD10/LPD10FINAL_rename_Somatic.rData")
LPD10E <- subset(x=SomaticCell3, stage =="LPD10_EXP")
DimPlot(LPD10E, label = TRUE)
PD90_sc <- LPD10E
table(PD90_sc@active.ident)
PD90_sc$subclass <- PD90_sc@active.ident
DimPlot(PD90_sc, group.by = "subclass", label = TRUE)
allen_reference <- SCTransform(PD90_sc, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)
DimPlot(allen_reference, group.by = "subclass", label = TRUE)
ovary_2M_2 <- SCTransform(ovary_2M_2, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE)

anchors <- FindTransferAnchors(reference = allen_reference, query = ovary_2M_2, normalization.method = "SCT")

predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE,
                                  weight.reduction = ovary_2M_2[["pca"]], dims = 1:30, k.weight = 40)

ovary_2M_2[["predictions"]] <- predictions.assay
DefaultAssay(ovary_2M_2) <- "predictions"
a <- SpatialFeaturePlot(ovary_2M_2, features = c("cGC","mGC","Mesenchyma","pMesenchyma"), pt.size.factor = 6, ncol = 2, crop = TRUE)
ggsave("LPD10E-STSC-cGC.tiff",plot = a, device = "tiff", width = 10, height = 10, dpi = 300)
ggsave("LPD10E-STSC-cGC.pdf",plot = a, device = "pdf", width = 10, height = 10)












