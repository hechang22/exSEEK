library(dplyr)

#获取counts时生成的id列表文件
sample_ids = read.table("/WORK/hechang/short/sample_ids.txt",sep="\n")
#用第一个样本数据初始化结果矩阵
miR = read.table(paste0("/WORK/hechang/short/counts/miR/",sample_ids[1,1],".txt"),header=F,sep="\t")
piR = read.table(paste0("/WORK/hechang/short/counts/piR/",sample_ids[1,1],".txt"),header=F,sep="\t")

#遍历每一个样本id，与res矩阵列合并
for(i in 2:length(sample_ids[,1])){
    t.m = read.table(paste0("/WORK/hechang/short/counts/miR/",sample_ids[i,1],".txt"),header=F,sep="\t")
    miR = bind_cols(miR,t.m[,2])
    t.p = read.table(paste0("/WORK/hechang/short/counts/piR/",sample_ids[i,1],".txt"),header=F,sep="\t")
    piR = bind_cols(piR,t.p[,2])
}

rownames(miR) = miR[,1]
miR = miR[,-1]
rownames(piR) = piR[,1]
piR = piR[,-1]

colnames(miR)=sample_ids[,1]
colnames(piR)=sample_ids[,1]

miR = t(miR)
piR = t(piR)

#输出
write.csv(miR,file = "/WORK/hechang/short/miR.count.csv")
write.csv(piR,file = "/WORK/hechang/short/piR.count.csv")
write.csv(bind_cols(miR,piR),file = "/WORK/hechang/short/combined.count.csv")