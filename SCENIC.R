




rm(list = ls())
#LPD10 SCENIC
load("D:/mousemodel/LPD/LPD10/LPD10FINAL_rename_Somatic.rData")
library(Seurat)
DimPlot(SomaticCell3,label = TRUE,group.by = "stage")
DimPlot(SomaticCell3,label = TRUE)
LPD10C <- subset(x=SomaticCell3, stage =="LPD10_CON")
table(LPD10C@active.ident)

SCENIC_LPD10C <- subset(LPD10C, idents = c("cGC","mGC","Early_Theca","Theca","Mesenchyma","pMesenchyma","Endo","Immune_cell"))
table(SCENIC_LPD10C@active.ident)
exprMat  <-  as.matrix(SCENIC_LPD10C@assays$RNA@data)
dim(exprMat)
exprMat[1:4,1:4]
SCENIC_LPD10C@meta.data$CellType <- SCENIC_LPD10C@active.ident
cellInfo <- SCENIC_LPD10C@meta.data[,c(9,3,2)]
colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')
head(cellInfo)
table(cellInfo$CellType)


library(SCENIC)

scenicOptions <- initializeScenic(org="mgi", dbDir="cisTarget_databases", dbs= "mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather", nCores=1)
mm10_dbs <- list('500bp' = 'mm10__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather',
                 '10kb' = 'mm10__refseq-r80__10kb_up_and_down_tss.mc9nr.feather')
db_mcVersion <- 'v9'
db_path <- 'D:/mousemodel/LPD/LPD10/cisTarget_databases'
scenicOptions@settings$dbs <- mm10_dbs
scenicOptions@settings$dbDir <- db_path
scenicOptions@settings$db_mcVersion <- db_mcVersion
saveRDS(scenicOptions, file="int/scenicOptions.Rds")







### Co-expression network
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
dim(exprMat_filtered)
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1)
runGenie3(exprMat_filtered_log, scenicOptions)


### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) # Toy run settings

library(doParallel)
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings
export2loom(scenicOptions, exprMat)
saveRDS(scenicOptions, file="int/scenicOptions.Rds")



rm(list = ls()) 
library(Seurat) 
library(SCENIC)
library(doParallel)

scenicOptions=readRDS(file="int/scenicOptions.Rds")
motifEnrichment_selfMotifs_wGenes <- loadInt(scenicOptions, "motifEnrichment_selfMotifs_wGenes") 
as.data.frame(sort(table(motifEnrichment_selfMotifs_wGenes$highlightedTFs),decreasing = T))
tableSubset <- motifEnrichment_selfMotifs_wGenes[highlightedTFs=="Maff"]
viewMotifs(tableSubset)
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Runx1" & highConfAnnot==TRUE]
viewMotifs(tableSubset)



rm(list = ls()) 
library(Seurat) 
library(SCENIC)
library(doParallel)
library(SCopeLoomR)
scenicOptions=readRDS(file="int/scenicOptions.Rds")
scenicLoomPath <- getOutName(scenicOptions, "loomFile")
loom <- open_loom(scenicLoomPath)
regulons_incidMat <- get_regulons(loom)
regulons <- regulonsToGeneLists(regulons_incidMat)
regulonsAUC <- get_regulons_AUC(loom)
regulonsAucThresholds <- get_regulon_thresholds(loom)
embeddings <- get_embeddings(loom)


regulons <- loadInt(scenicOptions, "regulons")
regulons[c("Trps1", "Tagln2")]
regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))
regulonTargetsInfo <- loadInt(scenicOptions, "regulonTargetsInfo")
tableSubset <- regulonTargetsInfo[TF=="Tagln2" & highConfAnnot==TRUE]
viewMotifs(tableSubset)
write.csv(regulonTargetsInfo,file = "LPD10Ctop5_regulonTargetsInfo2.csv")


regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byCellType <- sapply(split(rownames(cellInfo), cellInfo$CellType),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byCellType_Scaled <- t(scale(t(regulonActivity_byCellType), center = T, scale=T))
library(pheatmap)
p2<- pheatmap::pheatmap(regulonActivity_byCellType_Scaled, #fontsize_row=3,
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(-2, 2, length.out = 100),
                   treeheight_row=15,treeheight_col = 15, border_color=NA, cluster_cols = T,cluster_rows = T,cellwidth = 30)

ggsave("Fig4A.tiff", plot=p2, width=6, height = 12)
topRegulators <- reshape2::melt(regulonActivity_byCellType_Scaled)
colnames(topRegulators) <- c("Regulon", "CellType", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)
write.csv(topRegulators,file = "LPD10Ctopregulators2.csv")

minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byCellType_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$CellType),
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byCellType_Binarized[which(rowSums(regulonActivity_byCellType_Binarized>minPerc)>0),]
pheatmap::pheatmap(binaryActPerc_subset, # fontsize_row=5,
                   color = colorRampPalette(c("white","pink","red"))(100), breaks=seq(0, 1, length.out = 100),
                   treeheight_row=20, cluster_cols = T, border_color=NA)


#Seurat可视化
AUCmatrix <- readRDS("int/3.4_regulonAUC.Rds")
AUCmatrix <- AUCmatrix@assays@data@listData$AUC
AUCmatrix <- data.frame(t(AUCmatrix), check.names=F)
RegulonName_AUC <- colnames(AUCmatrix)
RegulonName_AUC <- gsub(' \\(','_',RegulonName_AUC)
RegulonName_AUC <- gsub('\\)','',RegulonName_AUC)
colnames(AUCmatrix) <- RegulonName_AUC
scRNAauc <- AddMetaData(SCENIC_LPD10C, AUCmatrix)
scRNAauc@assays$integrated <- NULL
saveRDS(scRNAauc,'scRNAauc_LPD10C.rds')
##导入二进制regulonAUC矩阵
BINmatrix <- readRDS("int/4.1_binaryRegulonActivity.Rds")
BINmatrix <- data.frame(t(BINmatrix), check.names=F)
RegulonName_BIN <- colnames(BINmatrix)
RegulonName_BIN <- gsub(' \\(','_',RegulonName_BIN)
RegulonName_BIN <- gsub('\\)','',RegulonName_BIN)
colnames(BINmatrix) <- RegulonName_BIN
scRNAbin <- AddMetaData(SCENIC_LPD10C, BINmatrix)
scRNAbin@assays$integrated <- NULL
#saveRDS(scRNAbin, 'scRNAbin.rds')
#绘图
p1<-FeaturePlot(scRNAauc, features='Fos_extended_52g', reduction = 'umap',pt.size = 1)
p2<-FeaturePlot(scRNAbin, features='Nfkb2_13g', reduction = 'umap',pt.size = 1)
p1
p2
ggsave("Fig4B-Nfkb2.tiff", plot=p2, width=4.5, height = 4)











