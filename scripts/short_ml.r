library(glmnet)
library(caret)
library(pROC)
library(magrittr)

### 读取与scale
x = read.table("./short/log.tmm.matrix.combat.txt",sep="\t") %>% t() %>% scale(center=T,scale=T)
#any(is.na(x)) #FALSE
#无缺失值
y = read.table("./exRNA-short/metadata.txt",header = T,sep="\t")$label[1:150]
y = factor(y)

### 划分数据集
set.seed(123) 
# 80%训练集
train.indices = createDataPartition(y,p=0.8,times = 1,list=T)$Resample1
x.train = x[train.indices,]
x.test = x[ -train.indices,]
y.train = y[train.indices]
y.test = y[ -train.indices]

### 基于LR的分类
## 特征选择（基于RFE）
# 初始设置
rfFuncs$summary = twoClassSummary 

#可供选择的特征个数
rfe.sizes = c(seq(10,100,10),200,500,700)

# 通过bootstrapping(有放回的抽样)进行性能评估
rfectrl = rfeControl(functions=rfFuncs,
                      verbose = TRUE,
                      method="cv",number=10)
rfe.results = rfe(x.train,y.train, 
               sizes=rfe.sizes, 
               rfeControl=rfectrl,
               metric = "ROC")

png("./fig/short.rfe.res.png")    
plot(rfe.results)
dev.off()

saveRDS(rfe.results,"./short/rfe.results.RDS")
write.csv(predictors(rfe.results),"./short/predictors.csv")


## 调参
params.grid.lr = expand.grid(alpha = c(0,0.5,1),lambda = c(0,0.01,0.1,1))
params.grid.svm = expand.grid(C = c(0.1,1,10),sigma = c(1,2,3,4,5))

tr.ctrl = trainControl(method="cv",
                        number = 5,
                        summaryFunction = twoClassSummary,
                        classProbs = TRUE)
cv.fitted.lr = train(x.train[,predictors(rfe.results)],y.train,
                   method="glmnet",
                   family="binomial",
                   metric = "ROC",
                   tuneGrid = params.grid.lr,
                   preProcess = NULL,
                   trControl = tr.ctrl )
#cv.fitted.lr$bestTune
# alpha lambda
#6   1   0.01

ranger = train(x.train,y.train,method="ranger",trControl=tr.ctrl) #ranger,快速RF
y.test.pred.ranger = predict(ranger,newdata=x.test,type="prob")

y.test.pred.prob.lr = predict(cv.fitted.lr,newdata=x.test,type="prob")
y.test.pred.prob.svm = predict(cv.fitted.svm,newdata=x.test,type="prob")

roc.curve.lr = roc(y.test,y.test.pred.prob.lr[,2])
roc.ranger = roc(y.test,y.test.pred.ranger[,2])

png("./fig/short.roc.png")
plot(roc.curve.lr,col="black",print.auc=T)
plot(roc.ranger,add=T,col="red")
legend("bottomright", legend=c("LR","Ranger"),col=c("black","red"),lty=1)
dev.off()
