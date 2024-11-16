## bash
## conda activate r-env
##

library(dplyr)

#获取counts时生成的id列表文件
sample_ids = read.table("/WORK/hechang/long/sample_ids.txt",sep="\n")
#声明结果矩阵
res = data.frame(temp=rep(0,11),
    row.names = c('mRNA','lncRNA','snoRNA','snRNA','srpRNA','tRNA','YRNA','intron','pseudogene','unassigned','total'))
#遍历每一个样本id，与res矩阵列合并
for(i in 1:length(sample_ids[,1])){
    t = read.table(paste0("/WORK/hechang/long/assign/",sample_ids[i,1],".txt"),header=F,sep="\t")
    res = bind_cols(res,t[,2])
}
#删除声明时的临时列，每一列为一个样本
res = select(res,-temp)
colnames(res) = sample_ids[,1]

#输出
write.csv(res,file = "/WORK/hechang/long/sort.matrix.csv")

### 画图

sort = read.csv("/WORK/hechang/long/sort.matrix.csv",row.names=1)
sum = rowSums(sort) 
piepercent = paste0(round(100*sum[1:10]/sum[11]), "%")
colors = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', 
           '#D6E7A3', '#57C3F3', '#476D87', '#D5D9E5',
           '#E59CC4', '#AB3282', '#23452F', '#BD956A')

png("./fig/pie.png")
pie(sum[1:10], labels = piepercent,main = "各类RNA比例",col=colors,family='GB1')
legend("topright", names(sum)[1:10], cex=0.8,fill=colors)
dev.off()
