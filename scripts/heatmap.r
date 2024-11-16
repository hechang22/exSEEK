library(pheatmap)
library(ggplot2)
library(magrittr)

long = read.table("./long/log.tmm.matrix.combat.txt",sep="\t") %>% t() %>% scale(center=T,scale=T)
short = read.table("./short/log.tmm.matrix.combat.txt",sep="\t") %>% t() %>% scale(center=T,scale=T)

long.feature = read.csv("./long/max20.csv",row.names=1) %>% unlist()
short.feature = read.csv("./short/predictors.csv",row.names=1) %>% unlist()
short.feature = short.feature[-1]

long.annotation = data.frame(label = read.table("./exRNA-long/metadata.txt",row.names = 1,header = T)$label)
rownames(long.annotation) = rownames(long)
short.annotation = data.frame(label = read.table("./exRNA-short/metadata.txt",row.names=1,header = T,sep="\t")$label[1:150])
rownames(short.annotation) = rownames(short)

p.l = pheatmap(long[,long.feature],show_rownames=F,annotation_row = long.annotation)
p.s = pheatmap(short[,short.feature],show_rownames=F,annotation_row = short.annotation)

png("./fig/heatmap.long.png",width=600,height=1000)
p.l
dev.off()

png("./fig/heatmap.short.png",width=600,height=1000)
p.s
dev.off()
