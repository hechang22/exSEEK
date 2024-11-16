library(glmnet)
library(caret)
library(pROC)
library(magrittr)

### 读取与scale
x = read.table("./long/log.tmm.matrix.combat.txt",sep="\t") %>% t() %>% scale(center=T,scale=T)
pca = read.csv("./long/pca.lower.matrix.csv",row.names=1)
#any(is.na(data)) #FALSE
#无缺失值
y = read.table("./exRNA-long/metadata.txt",row.names = 1,header = T)$label
y[which(y!="NC")]="C"
y = factor(y)

### 划分数据集
set.seed(123) 
# 80%训练集
train.indices = createDataPartition(y,p=0.8,times = 1,list=T)$Resample1
x.train = x[train.indices,]
x.test = x[ -train.indices,]
y.train = y[train.indices]
y.test = y[ -train.indices]
pca.train = pca[train.indices,]
pca.test = pca[-train.indices,]

### 基于LR的分类
## 特征选择（基于RFE）
# 初始设置
rfFuncs$summary = twoClassSummary 

#可供选择的特征个数
rfe.sizes = c(seq(10,90,10),seq(100,900,100))

# 通过bootstrapping进行性能评估
rfectrl = rfeControl(functions=rfFuncs,
                      verbose = TRUE,
                      method="cv",number=6)
rfe.results = rfe(x.train,y.train, 
               sizes=rfe.sizes, 
               rfeControl=rfectrl,
               metric = "ROC")
png("./fig/rfe.res.png")    
plot(rfe.results)
dev.off()

saveRDS(rfe.results,"./long/rfe.res.RDS")
#rfe.results = readRDS("./long/rfe.res.RDS")

## 调参LR
params.grid = expand.grid(alpha = c(0,0.1,0.3,0.5,0.7,1),lambda = c(0,0.01,0.1,0.3,0.5,1))

tr.ctrl = trainControl(method="cv",
                        number = 5,
                        summaryFunction = twoClassSummary,
                        classProbs = TRUE)
cv.fitted.lr = train(x.train[,predictors(rfe.results)],y.train,
                   method="glmnet",
                   family="binomial",
                   metric = "ROC",
                   tuneGrid = params.grid,
                   preProcess = NULL,
                   trControl = tr.ctrl )

saveRDS(cv.fitted.lr,"./long/lr.RDS")
cv.fitted.lr$bestTune
#cv.fitted.lr = readRDS("./long/lr.RDS")
y.test.pred.lr = predict(cv.fitted.lr,newdata=x.test,type="prob")

## 调参LR(PCA)

cv.fitted.lr.pca = train(pca.train,y.train,
                   method="glmnet",
                   family="binomial",
                   metric = "ROC",
                   tuneGrid = params.grid,
                   preProcess = NULL,
                   trControl = tr.ctrl )
saveRDS(cv.fitted.lr.pca,"./long/lr.pca.RDS")
#cv.fitted.lr.pca=readRDS("./long/lr.pca.RDS")
y.test.pred.lr.pca = predict(cv.fitted.lr.pca,newdata=pca.test,type="prob")

## PCA贡献最大的20个特征进行训练 LR

#计算对20维度pca总贡献最多的20个特征
max20 = pca %>% abs() %>% rowSums() %>% order(decreasing=T) %>% head(20)
write.csv(max20,"./long/max20.csv")
#训练
cv.fitted.lr.select = train(x.train[,max20],y.train,
                   method="glmnet",
                   family="binomial",
                   metric = "ROC",
                   tuneGrid = params.grid,
                   preProcess = NULL,
                   trControl = tr.ctrl )
saveRDS(cv.fitted.lr.select,"./long/lr.pca.select.RDS")
#cv.fitted.lr.select=readRDS("./long/lr.pca.select.RDS")
y.test.pred.lr.select = predict(cv.fitted.lr.select,newdata=x.test[,max20],type="prob")


## ranger
ranger = train(x.train,y.train,method="ranger",trControl=tr.ctrl) #ranger,快速RF
saveRDS(ranger,"./long/ranger.RDS")
#ranger = readRDS("./long/ranger.RDS")
y.test.pred.ranger = predict(ranger,newdata=x.test,type="prob")


roc.lr = roc(y.test,y.test.pred.lr[,2])
roc.ranger = roc(y.test,y.test.pred.ranger[,2])
roc.lr.pca = roc(y.test,y.test.pred.lr.pca[,2])
roc.lr.select = roc(y.test,y.test.pred.lr.select[,2])
png("./fig/roc.png")
plot(roc.lr,col="black")
plot(roc.ranger,add=T,col="red")
plot(roc.lr.pca,add=T,col="blue",print.auc=T)
legend("bottomright", legend=c("LR","Ranger","LR+PCA"),col=c("black","red","blue"),lty=1)
dev.off()
png("./fig/roc.selected.png")
plot(roc.lr.select,print.auc=T)
dev.off()
