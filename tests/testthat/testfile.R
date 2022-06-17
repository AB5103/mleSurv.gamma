library(mleSurv.gamma)
library(ucminf)
library(survival)
library(xtable)

D=veteran
D$prior <- factor(as.character(D$prior), labels = c(0, 1))
D$trt <- factor(as.character(D$trt), labels = c(0, 1))
D$celltype=as.character(D$celltype)
D=data.frame(D$karno,D$trt,D$diagtime,D$prior,D$age,D$celltype,D$time,D$status)
colnames(D)=c("karno","trt","dtime","prior","age","celltype","time","event")
dist="Exponential"
cluster="celltype"
formula=as.formula("Surv(time, event) ~ karno")
opt.model(formula,D,cluster)
R.model=bcmle(formula,D,cluster,dist)
plotSurv(R.model)
