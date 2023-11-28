



library(DESeq2)


#construct DEseq2 object
mydata <- read.csv("D:/20230420/低蛋白饮食课题汇报/Figure5/oocytePD21/mmLPD.csv")
class(mydata)
countdata <- as.matrix(mydata[2:7])
row.names(countdata) <- mydata$Geneid
condition <- factor(c("CON","CON","CON","KO","KO","KO"))
coldata <- data.frame(row.names=colnames(countdata), condition)
condition
coldata
dds <- DESeqDataSetFromMatrix(countData=countdata, colData = coldata, design= ~ condition )
head(dds)



dds2 <- DESeq(dds)
resultsNames(dds2)
res <- results(dds2)
summary(res)
table(res$pvalue<0.05)
table(res$padj<0.05)
res <- res[order(res$padj),]
diff_gene_deseq2 <-subset(res,padj < 0.05 & (log2FoldChange > 1 | log2FoldChange < -1))
diff_gene_deseq2 <- row.names(diff_gene_deseq2)
resdata <-  merge(as.data.frame(res),as.data.frame(counts(dds2,normalize=TRUE)),by="row.names",sort=FALSE)

up_diff <- subset(resdata,pvalue < 0.05 & log2FoldChange >1)
down_diff <- subset(resdata,pvalue < 0.05 & log2FoldChange < -1)


up_diff <- subset(resdata,padj < 0.05 & log2FoldChange >1)
down_diff <- subset(resdata,padj < 0.05 & log2FoldChange < -1)
write.csv(resdata,file= "LPD10_DEseq2_padj_resdata.csv")
write.csv(diff_gene_deseq2,file= "LPD10_DEseq2_padj_diff_gene_deseq2.csv")
write.csv(up_diff,file= "LPD10_DEseq2_pvalue_diff_gene_up.csv")
write.csv(down_diff,file= "LPD10_DEseq2_pvalue_diff_gene_down.csv")



library(ggplot2)
install.packages("ggpubr")
library(ggpubr)
install.packages("ggthemes")
library(ggthemes)
input_oocyte <- read.csv("D:/mousemodel/LPD/LPD10/LPD10_DEseq2_padj_resdata.csv",header = T)
head(input_oocyte)
input_oocyte$logp <- -log10(input_oocyte$pvalue)
ggscatter(input_oocyte, x="log2FoldChange", y="logp")+theme_base()
input_oocyte$Group = "not-significant"
input_oocyte$Group[which((input_oocyte$pvalue<0.05)&input_oocyte$log2FoldChange>0.25)]="up-regulated"
input_oocyte$Group[which((input_oocyte$pvalue<0.05)&input_oocyte$log2FoldChange< -0.25)]="down-regulated"
ggscatter(input_oocyte, x="log2FoldChange", y="logp",color = "Group")+theme_base()

input_oocyte$Label = ""
input_oocyte <- input_oocyte[order(input_oocyte$pvalue),]
up.genes <- head(input_oocyte$SYMBOL[which(input_oocyte$Group == "up-regulated")], 6)
down.genes <- head(input_oocyte$SYMBOL[which(input_oocyte$Group == "down-regulated")], 6)
oocyte.top10.genes <- c(as.character(up.genes), as.character(down.genes))
input_oocyte$Label[match(oocyte.top10.genes, input_oocyte$SYMBOL)] <- oocyte.top10.genes


ggscatter(input_oocyte, x="log2FoldChange", y="logp",color = "Group", palette = c("#3A8DCC","#A9A9A9", "#E60012"), 
          size = 1,
          label = input_oocyte$Label,
          font.label = 12,
          repel = T,
          xlab = "log2FoldChange",
          ylab = "-log10(P-value)")+theme_base()+
  geom_hline(yintercept = 1.3,linetype="dashed")+
  geom_vline(xintercept = c(-0.25, 0.25), linetype="dashed")


p <- ggscatter(input_oocyte, x="log2FoldChange", y="logp",color = "Group", palette = c("#427fc2","#A9A9A9", "#cf5d27"), 
               size = 1,
               repel = T,
               xlab = "avg_logFC",
               ylab = "-log10(P-value)")+theme_base()
ggsave("oocyte_full1.pdf", plot=p, width=8, height = 6)

library(Seurat)
load("D:/mousemodel/LPD/LPD10/LPD10FINAL_rename_Somatic.rData")
Immune.cells <- subset(SomaticCell3, idents = "Immune_cell")
DEG_Immune_Full <- FindMarkers(Immune.cells, ident.1 = "LPD10_EXP", group.by = "stage",logfc.threshold = 0, min.pct = 0)
write.csv(DEG_Immune_Full, file = "DEG_Immune_full.csv")

library(ggpubr)
library(ggthemes)
input_cGC <- read.csv("D:/mousemodel/LPD/LPD10/DEG_Theca_full.csv",header = T)
head(input_cGC)
input_cGC$logp <- -log10(input_cGC$p_val)
ggscatter(input_cGC, x="avg_log2FC", y="logp")+theme_base()
input_cGC$Group = "not-significant"
input_cGC$Group[which((input_cGC$p_val<0.05)&input_cGC$avg_log2FC>0.25)]="up-regulated"
input_cGC$Group[which((input_cGC$p_val<0.05)&input_cGC$avg_log2FC< -0.25)]="down-regulated"
ggscatter(input_cGC, x="avg_log2FC", y="logp",color = "Group")+theme_base()

input_cGC$Label = ""
input_cGC <- input_cGC[order(input_cGC$p_val),]
up.genes <- head(input_cGC$SYMBOL[which(input_cGC$Group == "up-regulated")], 6)
down.genes <- head(input_cGC$SYMBOL[which(input_cGC$Group == "down-regulated")], 6)
cGC.top10.genes <- c(as.character(up.genes), as.character(down.genes))
input_cGC$Label[match(cGC.top10.genes, input_cGC$SYMBOL)] <- cGC.top10.genes


ggscatter(input_cGC, x="avg_log2FC", y="logp",color = "Group", palette = c("#3A8DCC","#A9A9A9", "#E60012"), 
          size = 1,
          label = input_cGC$Label,
          font.label = 12,
          repel = T,
          xlab = "log2FoldChange",
          ylab = "-log10(P-value)")+theme_base()+
  geom_hline(yintercept = 2,linetype="dashed")+
  geom_vline(xintercept = c(-0.25, 0.25), linetype="dashed")


p <- ggscatter(input_cGC, x="avg_log2FC", y="logp",color = "Group", palette = c("#427fc2","#A9A9A9", "#cf5d27"), 
               size = 1,
               repel = T,
               xlab = "avg_logFC",
               ylab = "-log10(P-value)")+theme_base()
  

ggsave("Immune_full1.pdf", plot=p, width=8, height = 6)




GV <- read.csv("D:/mousemodel/T1D/T1D/vEndo/UPDEG_GO_shai.csv",header = T)
ego@result <- GV

p <- barplot(ego, orderBy = "x", showCategory=10)
p + scale_y_discrete(labels=function(x) str_wrap(x, width=40))









