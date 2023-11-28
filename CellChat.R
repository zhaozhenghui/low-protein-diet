
rm(list = ls())
#cellchat-LPD10
install.packages('NMF')
install.packages('rlang')
devtools::install_github("jokergoo/circlize")
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("sqjin/CellChat")

library(NMF)
library(circlize)
library(ComplexHeatmap)
library(CellChat)
install.packages('mindr')
library(mindr)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(SeuratData)

load("D:/mousemodel/LPD/LPD10/LPD10FINAL_rename_Somatic.rData")
DimPlot(SomaticCell3, reduction = "umap", label = T)
options(stringsAsFactors = FALSE)
LPD10_subset <- subset(SomaticCell3, idents = c("cGC","mGC","Early_Theca","Theca","Mesenchyma","pMesenchyma","Endo","Immune_cell"))

LPD10C <- subset(x=LPD10_subset, stage =="LPD10_CON")
LPD10E <- subset(x=LPD10_subset, stage =="LPD10_EXP")

table(LPD10C@active.ident)
save(LPD10C, file = "LPD10C_cellchat1.rData")
save(LPD10E, file = "LPD10E_cellchat1.rData")


LPD10C.input <- LPD10C@assays$RNA@data
x <- LPD10C@active.ident
table(x)
x <- factor(x,levels=c("cGC","mGC","Early_Theca","Theca","Mesenchyma","pMesenchyma","Endo","Immune_cell"))
LPD10C@meta.data$seurat_annotation <-x
table(LPD10C$seurat_annotation)
LPD10C.identity = data.frame(group = LPD10C$seurat_annotation, row.names = names(LPD10C$seurat_annotation))
unique(LPD10C.identity$group)
LPD10C.cellchat <- createCellChat(object = LPD10C.input)
LPD10C.cellchat
summary(LPD10C.cellchat)
LPD10C.cellchat <- addMeta(LPD10C.cellchat, meta = LPD10C.identity, meta.name = "labels")
LPD10C.cellchat <- setIdent(LPD10C.cellchat, ident.use = "labels")
levels(LPD10C.cellchat@idents)
groupSize <- as.numeric(table(LPD10C.cellchat@idents))


LPD10E.input <- LPD10E@assays$RNA@data
x <- LPD10E@active.ident
table(x)
x <- factor(x,levels=c("cGC","mGC","Early_Theca","Theca","Mesenchyma","pMesenchyma","Endo","Immune_cell"))
LPD10E@meta.data$seurat_annotation <-x
table(LPD10E$seurat_annotation)
LPD10E.identity = data.frame(group = LPD10E$seurat_annotation, row.names = names(LPD10E$seurat_annotation))
unique(LPD10E.identity$group)
LPD10E.cellchat <- createCellChat(object = LPD10E.input)
LPD10E.cellchat
summary(LPD10E.cellchat)
LPD10E.cellchat <- addMeta(LPD10E.cellchat, meta = LPD10E.identity, meta.name = "labels")
LPD10E.cellchat <- setIdent(LPD10E.cellchat, ident.use = "labels")
levels(LPD10E.cellchat@idents)
groupSize <- as.numeric(table(LPD10E.cellchat@idents))



object.list <- list(LPD10C = LPD10C.cellchat, LPD10E = LPD10E.cellchat)
cellchat.merge <- mergeCellChat(object.list, add.names = names(object.list))


CellChatDB <- CellChatDB.mouse
(out3 <- capture.output(str(CellChatDB)))
out4 <- paste(out3, collapse = "\n")
#mm(gsub("\\$", "# ",gsub("\\.\\. ","#", out4)), type = "text")
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
LPD10C.cellchat@DB <- CellChatDB.use
unique(CellChatDB$interaction$annotation)





LPD10C.cellchat <- subsetData(LPD10C.cellchat)
LPD10C.cellchat <- identifyOverExpressedGenes(LPD10C.cellchat)
LPD10C.cellchat <- identifyOverExpressedInteractions(LPD10C.cellchat)
LPD10C.cellchat <- projectData(LPD10C.cellchat, PPI.mouse)


LPD10C.cellchat <- computeCommunProb(LPD10C.cellchat)
LPD10C.cellchat <- computeCommunProbPathway(LPD10C.cellchat)
LPD10C.cellchat <- aggregateNet(LPD10C.cellchat)
LPD10C.cellchat@netP$pathways
head(LPD10C.cellchat@LR$LRsig)


LPD10C.cellchat@netP$pathways
levels(LPD10C.cellchat@idents)
vertex.receiver = seq(1,4)
pathways.show <- "VEGF"
netVisual_aggregate(LPD10C.cellchat, signaling = pathways.show, vertex.receiver= vertex.receiver, vertex.size = groupSize)
netVisual_aggregate(LPD10C.cellchat, signaling = c("VEGF"), layout = "circle", vertex.size = groupSize, pt.title = 20, vertex.label.cex = 1.7)
netAnalysis_contribution(LPD10C.cellchat, signaling = pathways.show)



LPD10C.cellchat <- netAnalysis_computeCentrality(LPD10C.cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(LPD10C.cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


gg1 <- netAnalysis_signalingRole_scatter(LPD10C.cellchat)
gg2 <- netAnalysis_signalingRole_scatter(LPD10C.cellchat, signaling = c("VEGF", "TGFb"))
gg1 + gg2


ht1 <- netAnalysis_signalingRole_heatmap(LPD10C.cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(LPD10C.cellchat, pattern = "incoming")
ht1 + ht2

nPatterns = 5
LPD10C.cellchat <- identifyCommunicationPatterns(LPD10C.cellchat, pattern = "outgoing", k = nPatterns, height = 10, width = 8)

p2 <- netAnalysis_river(LPD10C.cellchat, pattern = "outgoing")
ggsave("Fig3a.tiff", plot=p2, width=8, height = 8)
netAnalysis_dot(LPD10C.cellchat, pattern = "outgoing")


library("reticulate")
#reticulate::py_install(packages ='umap-learn')




nPatterns = 5
LPD10C.cellchat <- identifyCommunicationPatterns(LPD10C.cellchat, pattern = "incoming", k = nPatterns, height = 10, width = 8)
p3 <- netAnalysis_river(LPD10C.cellchat, pattern = "incoming")
ggsave("Fig3b.tiff", plot=p3, width=8, height = 8)
netAnalysis_dot(LPD10C.cellchat, pattern = "incoming")
LPD10C.cellchat <- computeNetSimilarity(LPD10C.cellchat, type = "functional")
LPD10C.cellchat <- netEmbedding(LPD10C.cellchat, type = "functional")
LPD10C.cellchat <- netClustering(LPD10C.cellchat, type = "functional")
netVisual_embedding(LPD10C.cellchat, type = "functional", label.size = 3.5)
LPD10C.cellchat <- computeNetSimilarity(LPD10C.cellchat, type = "structural")
LPD10C.cellchat <- netEmbedding(LPD10C.cellchat, type = "structural")
LPD10C.cellchat <- netClustering(LPD10C.cellchat, type = "structural")
netVisual_embedding(LPD10C.cellchat, type = "structural", label.size = 3.5)
saveRDS(LPD10C.cellchat, file = "LPD10C_cellchat.rds")





LPD10E.input <- LPD10E@assays$RNA@data
x <- LPD10E@active.ident
table(x)
x <- factor(x,levels=c("cGC","mGC","Early_Theca","Theca","Mesenchyma","pMesenchyma","Endo","Immune_cell"))
LPD10E@meta.data$seurat_annotation <-x
table(LPD10E$seurat_annotation)
LPD10E.identity = data.frame(group = LPD10E$seurat_annotation, row.names = names(LPD10E$seurat_annotation))
unique(LPD10E.identity$group)
LPD10E.cellchat <- createCellChat(object = LPD10E.input)
LPD10E.cellchat
summary(LPD10E.cellchat)
LPD10E.cellchat <- addMeta(LPD10E.cellchat, meta = LPD10E.identity, meta.name = "labels")
LPD10E.cellchat <- setIdent(LPD10E.cellchat, ident.use = "labels")
levels(LPD10E.cellchat@idents)
groupSize <- as.numeric(table(LPD10E.cellchat@idents))



CellChatDB <- CellChatDB.mouse
(out3 <- capture.output(str(CellChatDB)))
out4 <- paste(out3, collapse = "\n")
#mm(gsub("\\$", "# ",gsub("\\.\\. ","#", out4)), type = "text")
colnames(CellChatDB$interaction)
CellChatDB$interaction[1:4,1:4]
head(CellChatDB$cofactor)
head(CellChatDB$complex)
head(CellChatDB$geneInfo)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
LPD10E.cellchat@DB <- CellChatDB.use
unique(CellChatDB$interaction$annotation)





LPD10E.cellchat <- subsetData(LPD10E.cellchat)
LPD10E.cellchat <- identifyOverExpressedGenes(LPD10E.cellchat)
LPD10E.cellchat <- identifyOverExpressedInteractions(LPD10E.cellchat)
LPD10E.cellchat <- projectData(LPD10E.cellchat, PPI.mouse)


LPD10E.cellchat <- computeCommunProb(LPD10E.cellchat)
LPD10E.cellchat <- computeCommunProbPathway(LPD10E.cellchat)
LPD10E.cellchat <- aggregateNet(LPD10E.cellchat)
LPD10E.cellchat@netP$pathways
head(LPD10E.cellchat@LR$LRsig)


LPD10E.cellchat@netP$pathways
levels(LPD10E.cellchat@idents)
vertex.receiver = seq(1,4)
pathways.show <- "PDGF"
netVisual_aggregate(LPD10E.cellchat, signaling = pathways.show, vertex.receiver= vertex.receiver, vertex.size = groupSize)
netVisual_aggregate(LPD10E.cellchat, signaling = c("TGFb"), layout = "circle", vertex.size = groupSize, pt.title = 20, vertex.label.cex = 1.7)
netAnalysis_contribution(LPD10E.cellchat, signaling = pathways.show)



LPD10E.cellchat <- netAnalysis_computeCentrality(LPD10E.cellchat, slot.name = "netP")
netAnalysis_signalingRole_network(LPD10E.cellchat, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)


gg1 <- netAnalysis_signalingRole_scatter(LPD10E.cellchat)
gg2 <- netAnalysis_signalingRole_scatter(LPD10E.cellchat, signaling = c("PDGF", "TGFb"))
gg1 + gg2


ht1 <- netAnalysis_signalingRole_heatmap(LPD10E.cellchat, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(LPD10E.cellchat, pattern = "incoming")
ht1 + ht2

nPatterns = 5
LPD10E.cellchat <- identifyCommunicationPatterns(LPD10E.cellchat, pattern = "outgoing", k = nPatterns, height = 10, width = 8)

p4 <- netAnalysis_river(LPD10E.cellchat, pattern = "outgoing")
ggsave("Fig3c.tiff", plot=p4, width=8, height = 8)

netAnalysis_dot(LPD10E.cellchat, pattern = "outgoing")



nPatterns = 5
LPD10E.cellchat <- identifyCommunicationPatterns(LPD10E.cellchat, pattern = "incoming", k = nPatterns, height = 10, width = 8)
netAnalysis_river(LPD10E.cellchat, pattern = "incoming")
netAnalysis_dot(LPD10E.cellchat, pattern = "incoming")
LPD10E.cellchat <- computeNetSimilarity(LPD10E.cellchat, type = "functional")
LPD10E.cellchat <- netEmbedding(LPD10E.cellchat, type = "functional")
LPD10E.cellchat <- netClustering(LPD10E.cellchat, type = "functional")
netVisual_embedding(LPD10E.cellchat, type = "functional", label.size = 3.5)
LPD10E.cellchat <- computeNetSimilarity(LPD10E.cellchat, type = "structural")
LPD10E.cellchat <- netEmbedding(LPD10E.cellchat, type = "structural")
LPD10E.cellchat <- netClustering(LPD10E.cellchat, type = "structural")
netVisual_embedding(LPD10E.cellchat, type = "structural", label.size = 3.5)
saveRDS(LPD10E.cellchat, file = "LPD10E_cellchat.rds")


object.list <- list(LPD10C = LPD10C.cellchat, LPD10E = LPD10E.cellchat)
cellchat.merge <- mergeCellChat(object.list, add.names = names(object.list))

gg1 <- compareInteractions(cellchat.merge, show.legend = F, group = c(1,2), color.use = c('#386CB0',"#E7298A"))
gg2 <- compareInteractions(cellchat.merge, show.legend = F, group = c(1,2), color.use =  c('#386CB0',"#E7298A"), measure = "weight")
gg1 + gg2
ggsave("Fig3e.tiff", plot=gg2, width=4, height = 4)

par(mfrow = c(1,2), xpd=TRUE)
p6 <- netVisual_diffInteraction(cellchat.merge, weight.scale = T)
p7 <- netVisual_diffInteraction(cellchat.merge, weight.scale = T, measure = "weight")
ggsave("Fig3f.tiff", plot=p7, width=8, height = 4)

gg1 <- netVisual_heatmap(cellchat.merge)
gg2 <- netVisual_heatmap(cellchat.merge, measure = "weight")
gg1 + gg2

gg1 <- rankNet(cellchat.merge, mode = "comparison", stacked = T, do.stat = TRUE,color.use = c('#386CB0',"#E7298A"))
gg2 <- rankNet(cellchat.merge, mode = "comparison", stacked = F, do.stat = TRUE,color.use = c('#386CB0',"#E7298A"))
gg1 + gg2
ggsave("Fig3i.tiff", plot=gg2, width=4, height = 8)
netVisual_bubble(cellchat.merge, sources.use = 8, targets.use = c(1:8),  comparison = c(1, 2), angle.x = 45, color.text = c('#386CB0',"#E7298A"))


gg1 <- netVisual_bubble(cellchat.merge, sources.use = 8, targets.use = c(1:8),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in LPD", angle.x = 45, remove.isolate = T, color.text = c('#386CB0',"#E7298A"))
gg2 <- netVisual_bubble(cellchat.merge, sources.use = 8, targets.use = c(1:8),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LPD", angle.x = 45, remove.isolate = T, color.text = c('#386CB0',"#E7298A"))
gg1 + gg2


cellchat.merge@meta$datasets = factor(cellchat.merge@meta$datasets, levels = c("LPD10C", "LPD10E"))
plotGeneExpression(cellchat.merge, signaling = "GDF", split.by = "datasets", color.use = c('#386CB0',"#E7298A"))

saveRDS(cellchat.merge, file = "LPD10_cellchat_merge.rds")


LPD10_cellchat_merge <- readRDS("D:/mousemodel/LPD/LPD10/LPD10_cellchat_merge.rds")
netVisual_bubble(LPD10_cellchat_merge, sources.use = 5, targets.use = c(1:8),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in LPD", angle.x = 45, remove.isolate = T, color.text = c('#386CB0',"#E7298A"))









