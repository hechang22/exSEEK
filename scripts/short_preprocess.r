library(edgeR)
library(EDASeq)
library(sva)
library(ggplot2)
library(magrittr)
library(readr)

setwd("/WORK/hechang/")

#### FILTER & NORMALIZATION

#读入count矩阵(仅限CRC&HD)
count.matrix = read_delim("./short/combined.count.csv")
label = read.table("./exRNA-short/metadata.txt",header = T,sep="\t")$label
count.matrix = count.matrix[which(label=="HD"|label=="CRC"),]

#edgeR进行低表达基因过滤 0
y = DGEList(counts=t(count.matrix))
keep = filterByExpr(y)
y = y[keep, , keep.lib.sizes=FALSE]

#edgeR进行归一化（TMM&RLE）,对比
y.tmm = calcNormFactors(y,method="TMM")
y.rle = calcNormFactors(y,method="RLE")
cpm.tmm.matrix = cpm(y.tmm)   
cpm.rle.matrix = cpm(y.rle)

png("./fig/short.plotRLE.tmm.png")
plotRLE(cpm.tmm.matrix)
dev.off()

png("./fig/short.plotRLE.rle.png")
plotRLE(cpm.rle.matrix)
dev.off()
#区别同样不大

### 去除batch effect

#对数处理
log.tmm.matrix = cpm(y.tmm,log=T)
metadata = read.table("./exRNA-short/metadata.txt",header = T,sep = "\t")
# perform batch correction with ComBat
# only one batch allowed
log.tmm.matrix.combat = ComBat(log.tmm.matrix, batch=factor(metadata$library.prepration.day[which(label=="HD"|label=="CRC")]))
write.table(log.tmm.matrix.combat,"./short/log.tmm.matrix.combat.txt",sep="\t",quote = F)

#降维进行前后对比
pca.before = prcomp(log.tmm.matrix,center=F,scale=F) 
pca.after = prcomp(log.tmm.matrix.combat,center=F,scale=F)

before = data.frame(source = metadata$library.prepration.day[which(label=="HD"|label=="CRC")],
    PC1=t(log.tmm.matrix)%*%pca.before$x[,1],
    PC2=t(log.tmm.matrix)%*%pca.before$x[,2])
after = data.frame(source=metadata$library.prepration.day[which(label=="HD"|label=="CRC")],
    PC1=t(log.tmm.matrix.combat)%*%pca.after$x[,1],
    PC2=t(log.tmm.matrix.combat)%*%pca.after$x[,2])

png("./fig/short.PCA.before.png")
ggplot(before, aes(x=PC1,y=PC2,color=source))+geom_point()+theme_bw()
dev.off()
png("./fig/short.PCA.after.png")
ggplot(after, aes(x=PC1,y=PC2,color=source))+geom_point()+theme_bw()
dev.off()

