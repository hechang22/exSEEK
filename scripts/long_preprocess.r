library(edgeR)
library(EDASeq)
library(sva)
library(ggplot2)

setwd("/WORK/hechang/")

#### FILTER & NORMALIZATION

#读入count矩阵
count.matrix = read.table("./long/count.matrix.txt",sep = "\t",header = T)
rownames(count.matrix) = count.matrix[,1]
count.matrix = count.matrix[,-1]

#edgeR进行低表达基因过滤
y = DGEList(counts=count.matrix)
keep = filterByExpr(y)
y = y[keep, , keep.lib.sizes=FALSE]

#edgeR进行归一化（TMM&RLE）,对比
y.tmm = calcNormFactors(y,method="TMM")
y.rle = calcNormFactors(y,method="RLE")
cpm.tmm.matrix = cpm(y.tmm)   
cpm.rle.matrix = cpm(y.rle)

png("./fig/plotRLE.tmm.png")
plotRLE(cpm.tmm.matrix)
dev.off()

png("./fig/plotRLE.rle.png")
plotRLE(cpm.rle.matrix)
dev.off()
#区别不大

### 去除batch effect

#对数处理
log.tmm.matrix = cpm(y.tmm,log=T)
metadata = read.table("./exRNA-long/metadata.txt",row.names = 1,header = T)
# perform batch correction with ComBat
# only one batch allowed
log.tmm.matrix.combat = ComBat(log.tmm.matrix, batch=metadata$source)
write.table(log.tmm.matrix.combat,"./long/log.tmm.matrix.combat.txt",sep="\t",quote = F)

#降维进行前后对比
pca.before = prcomp(log.tmm.matrix,center=F,scale=F) 
pca.after = prcomp(log.tmm.matrix.combat,center=F,scale=F)

before = data.frame(source = metadata$source,
    PC1=t(log.tmm.matrix)%*%pca.before$x[,1],
    PC2=t(log.tmm.matrix)%*%pca.before$x[,2])
after = data.frame(source=metadata$source,
    PC1=t(log.tmm.matrix.combat)%*%pca.after$x[,1],
    PC2=t(log.tmm.matrix.combat)%*%pca.after$x[,2])

png("./fig/PCA.before.png")
ggplot(before, aes(x=PC1,y=PC2,color=source))+geom_point()+theme_bw()
dev.off()
png("./fig/PCA.after.png")
ggplot(after, aes(x=PC1,y=PC2,color=source))+geom_point()+theme_bw()
dev.off()

png("./fig/pca.elbowplot.png")
plot(pca.after$sdev[1:100])
dev.off()

lower.matrix = t(log.tmm.matrix.combat)%*%pca.after$x[,1:20]
write.csv(lower.matrix,"./long/pca.lower.matrix.csv")
