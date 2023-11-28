



rm(list = ls())
library(DOSE)
library(GO.db)
library(org.Mm.eg.db)
library(clusterProfiler)
table(SomaticCell3@active.ident)
cGC <- read.csv("D:/mousemodel/LPD/LPD10/GO/cGC.csv",header = T)
mGC <- read.csv("D:/mousemodel/LPD/LPD10/GO/mGC.csv",header = T)
Early_Theca <- read.csv("D:/mousemodel/LPD/LPD10/GO/Early_Theca.csv",header = T)
Theca <- read.csv("D:/mousemodel/LPD/LPD10/GO/Theca.csv",header = T)
Mesenchyma <- read.csv("D:/mousemodel/LPD/LPD10/GO/Mesenchyma.csv",header = T)
pMesenchyma <- read.csv("D:/mousemodel/LPD/LPD10/GO/pMesenchyma.csv",header = T)
Endo <- read.csv("D:/mousemodel/LPD/LPD10/GO/Endo.csv",header = T)
Periv <- read.csv("D:/mousemodel/LPD/LPD10/GO/Periv.csv",header = T)
Epi <- read.csv("D:/mousemodel/LPD/LPD10/GO/Epi.csv",header = T)
Immune_cell <- read.csv("D:/mousemodel/LPD/LPD10/GO/Immune.csv",header = T)

cGC_eg = bitr(cGC$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
mGC_eg = bitr(mGC$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
Early_Theca_eg = bitr(Early_Theca$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
Theca_eg = bitr(Theca$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
Mesenchyma_eg = bitr(Mesenchyma$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
pMesenchyma_eg = bitr(pMesenchyma$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
Endo_eg = bitr(Endo$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
Periv_eg = bitr(Periv$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
Epi_eg = bitr(Epi$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")
Immume_cell_eg = bitr(Immune_cell$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Mm.eg.db")



list_data <- list(cGC_eg$ENTREZID,mGC_eg$ENTREZID,Early_Theca_eg$ENTREZID,Theca_eg$ENTREZID,Mesenchyma_eg$ENTREZID,pMesenchyma_eg$ENTREZID,Endo_eg$ENTREZID,Periv_eg$ENTREZID,Epi_eg$ENTREZID,Immume_cell_eg$ENTREZID)
names(list_data) <- c("cGC", "mGC", "Early_Theca", "Theca", "Mesenchyma", "pMesenchyma", "Endo","Periv","Epi","Immume_cell")


lapply(list_data, head)
ck <- compareCluster(geneCluster = list_data, fun = "enrichGO", OrgDb="org.Mm.eg.db")

View(as.data.frame(ck))
write.csv(ck, file = "LPD10ALL_GO.csv")
dotplot(ck)
library(ggplot2)

PGCALL_GO <- read.csv("D:/mousemodel/LPD/LPD10/GO/LPD10ALL_GOshai.csv",header = T)
ck@compareClusterResult = PGCALL_GO

head(as.data.frame(ck))
p1 <- dotplot(ck,)
ggsave("Fig2e.tiff", plot=p1, width=8, height = 8)














