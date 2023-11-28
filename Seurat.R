


rm(list = ls())
library(Seurat) #version=4.0.1
library(patchwork)
library(devtools)
library(harmony)
library(ggplot2)
library(dplyr)
LPD10CON <- Read10X(data.dir = "D:/mousemodel/10%LPD/analysis/control/run_count_LPD10CON/outs/filtered_feature_bc_matrix")
LPD10EXP1 <- Read10X(data.dir = "D:/mousemodel/10%LPD/analysis/LPD1/run_count_LPD10EXP1/outs/filtered_feature_bc_matrix")
LPD10EXP2 <- Read10X(data.dir = "D:/mousemodel/10%LPD/analysis/LPD2/run_count_LPD10EXP2/outs/filtered_feature_bc_matrix")

CON <- CreateSeuratObject(counts = LPD10CON, project = "CON", min.cells = 3, min.features = 200)
EXP1 <- CreateSeuratObject(counts = LPD10EXP1, project = "EXP1", min.cells = 3, min.features = 200)
EXP2 <- CreateSeuratObject(counts = LPD10EXP2, project = "EXP2", min.cells = 3, min.features = 200)

CON$stage <- "CON"
EXP1$stage <- "EXP1"
EXP2$stage <- "EXP2"


LPD10 <- merge(x=CON, y=c(EXP1,EXP2))
LPD10[["percent.mt"]] <- PercentageFeatureSet(object = LPD10, pattern = "^mt-")
head(x=LPD10@meta.data, 5)


VlnPlot(object = LPD10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size=0, ncol = 3)
plot1 <- FeatureScatter(LPD10, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(LPD10, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
LPD10_A <- subset(LPD10, subset = nFeature_RNA > 1000 & nFeature_RNA < 3000 & percent.mt <5)
LPD10_A <- subset(LPD10, subset = nFeature_RNA > 1000 & nFeature_RNA < 4500 & percent.mt <5)
plot1 <- FeatureScatter(LPD10_A, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(LPD10_A, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))


LPD10norm <- NormalizeData(object = LPD10_A, normalization.method = "LogNormalize", scale.factor = 10000)
LPD_B <- FindVariableFeatures(LPD10norm, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(LPD_B), 10)
plot3 <- VariableFeaturePlot(LPD_B)
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE)
CombinePlots(plots = list(plot3, plot4))

all.genes <- rownames(LPD_B)
LPD_B <- ScaleData(LPD_B, features = all.genes)
LPD_B <- RunPCA(LPD_B, features = VariableFeatures(object = LPD_B))
print(LPD_B[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(LPD_B, dims = 1:2, reduction = "pca")
DimPlot(LPD_B, reduction = "pca")
DimHeatmap(LPD_B, dims = 1, cells = 500, balanced = TRUE)
DimHeatmap(LPD_B, dims = 1:8, cells = 500, balanced = TRUE)
LPD_B <- JackStraw(LPD_B, num.replicate = 100)
LPD_B <- ScoreJackStraw(LPD_B, dims = 1:20)
JackStrawPlot(LPD_B, dims = 1:15)
ElbowPlot(LPD_B)

LPD_B <- FindNeighbors(LPD_B, dims = 1:18)
LPD_B <- FindClusters(LPD_B, resolution = 0.2)
head(Idents(LPD_B), 5)

pbmcsca <- RunHarmony(LPD_B, group.by.vars = "orig.ident")
pbmcsca <- RunUMAP(pbmcsca, reduction = "harmony", dims = 1:18)
pbmcsca <- FindNeighbors(pbmcsca, reduction = "harmony", dims = 1:18) 
set.seed(11)
pbmcsca <- FindClusters(pbmcsca, resolution = 0.2)
DimPlot(pbmcsca, reduction = "umap", label = T)
DimPlot(pbmcsca, reduction = "umap", group.by = "orig.ident")

FeaturePlot(object = pbmcsca, features = c("Sohlh1", "Ooep","Dppa3","Zp3"),ncol = 2)
FeaturePlot(object = pbmcsca, features = c("Amhr2", "Foxl2","Fst","Bex1"),ncol = 2)
FeaturePlot(object = pbmcsca, features = c("Wnt4", "Wt1","Ihh","Htra1"),ncol = 2)
FeaturePlot(object = pbmcsca, features = c("Pecam1", "Cdh5","Cldn5","Vwf"),ncol = 2)
FeaturePlot(object = pbmcsca, features = c("Dcn", "Dlk1","Lum","Tcf21"),ncol = 2)
FeaturePlot(object = pbmcsca, features = c("Pclaf", "Ccnb1","Mki67","Top2a"),ncol = 2)
FeaturePlot(object = pbmcsca, features = c("Star", "Cyp17a1","Hsd3b1","Dlk1"),ncol = 2)
FeaturePlot(object = pbmcsca, features = c("Epcam", "Tm4sf5","Krt19","Lgals2"),ncol = 2)


LPD.markers <- FindAllMarkers(object = pbmcsca, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(LPD.markers,file = "test2.markers.csv")

x <- pbmcsca@active.ident
table(x)
x <- factor(x,levels=c("2","1","3","5","0","4","7","6","8","9"))
pbmcsca@meta.data$cellname <-x
Idents(pbmcsca) <-"cellname"
DimPlot(pbmcsca,label = TRUE)
new.cluster.ids <-c("Oocyte_Granulosa_Mix","Granulosa","Early_Theca","Theca","Mesenchyma","pMesenchyma","Endo","Immune_cell","Periv","Epi")
names(x=new.cluster.ids) <-levels(x=pbmcsca)
LPD10Final <-RenameIdents(object = pbmcsca, new.cluster.ids)
DimPlot(LPD10Final,label = TRUE)
DimPlot(LPD10Final, group.by = "stage")
LPD10Final$stage[LPD10Final$stage =="CON"] = "LPD10_CON"
LPD10Final$stage[LPD10Final$stage =="EXP1"] = "LPD10_EXP"
LPD10Final$stage[LPD10Final$stage =="EXP2"] = "LPD10_EXP"
save(LPD10Final, file = "LPD10FINAL_oocyte.rData")

plot1 = DimPlot(LPD10Final, reduction = "umap", label = T, label.size = 2)+ NoLegend()
plot2 = DimPlot(LPD10Final, reduction = "umap", group.by = "stage")+ NoLegend()
plot2 = DimPlot(LPD10Final, reduction = "umap", group.by = "stage")
plotc <- plot1+plot2
ggsave("LPD10Final_b.tiff", plot=plotc, width=7.7, height = 3) 


table(LPD10Final@active.ident)
NPD <- subset(x=LPD10Final, stage =="LPD10_CON")
LPD <- subset(x=LPD10Final, stage =="LPD10_EXP")
table(NPD@active.ident)
table(LPD@active.ident)


LPD10.markers <- FindAllMarkers(object = LPD10Final, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(LPD10.markers,file = "LPD10_Final.markers.csv")

top5 <- LPD10.markers %>% group_by(cluster) %>% top_n(n=5, wt = avg_log2FC)
DoHeatmap(LPD10Final, features = top5$gene) + NoLegend()

Immune_cell.cells <- subset(LPD10Final, idents = "Immune_cell")
DEG_Immune_cell <- FindMarkers(Immune_cell.cells, ident.1 = "LPD10_EXP", group.by = "stage",logfc.threshold = 0, min.pct = 0)
write.csv(DEG_Immune_cell,file = "DEG_Immune_cell_ALL.csv")


VlnPlot(Theca.cells, features = "Star", group.by = "stage",cols =  c('#386CB0',"#E7298A"),pt.size = 0,combine = FALSE)






