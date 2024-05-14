
# English:
# This R code analyzes pre-microRNAs' structural features, including correlation
# analysis, visualization creation, and predictive modeling using linear and logistic 
# regression. It evaluates model performance and conducts principal component analysis 
# (PCA) for feature reduction, followed by receiver operating characteristic (ROC) 
# analysis.
# 
# Spanish:
# Este código en R analiza las características estructurales de pre-microRNAs, 
# incluyendo análisis de correlación, creación de visualizaciones y modelado 
# predictivo utilizando regresión lineal y logística. Evalúa el rendimiento del 
# modelo y realiza un análisis de componentes principales (PCA) para reducción de
# características, seguido de análisis de la curva característica del receptor (ROC).
# 
# 


rm(list=ls())
load ("X_Local_ss_sRNA.RData")
X<-as.data.frame(X)
set.seed(7)
Xval<-X[sample(row.names(X), nrow(X)/2),]
X<-X[setdiff(row.names(X),row.names(Xval)),]

pdf(file = "Matrizh.pdf")

#calcular la correlaci?n entre dos columnas
cor(X$........., X$`)))))))))`)
plot(X$........., X$`)))))))))`)

#Correlaci?n de todas las variables entre todas las varables
cormat<-cor(X)
cormat[is.na(cormat)]<-0
heatmap(cormat, col=colorRampPalette(c("blue", "white", "red"))(50))

rownames(X)
class(X)
X$class<- ifelse(grepl("pre_miRNA", rownames(X)),1,0)
table(X$class)

cor(X$class, X$.........)


xcor<-NULL
for (i in seq(1,200) ) {
  print(i)
  
  xcor[i]<-  cor(X$class, X[,i], method = "spearman")
    
}

xcor[is.na(xcor)]<-0

range(xcor)
summary(xcor)
hist(xcor)

xcor[xcor == max(xcor)]
idmax<- which.max(xcor)
ids<-xcor>= 0.35
table(ids)


X[,c(colnames(X)[ids], "class")]
head(X[,c(colnames(X)[ids], "class")])

#Generar una tabla con x y y para prueba de error de entrenamiento
tab<-data.frame(X[,ids], y=X$class)
colnames(tab)<-c(colnames(X[,ids]),"y")

head(tab)
#tab<- tab[order(tab$x),]

summary(tab)
cor(tab)

plot(((tab)))
grid(tab)

modelL<-lm(formula = y~., data = tab)
summary(modelL)
plot(tab+
       matrix(rnorm(nrow(tab)*ncol(tab), sd = 0.1),
             nrow(tab),ncol(tab)),
     pch=19,
    
     col= ifelse(modelL$fitted.values>0.5,
                 rgb(1,0,0,0.5),rgb(0.1,0.1,0.1,0.1)))
hist(modelL$fitted.values)
#abline(h=0.5)

modelS<- glm (formula=y~., data=tab, family= binomial(link = "logit"))
summary(modelS)
plot(tab+
       matrix(rnorm(nrow(tab)*ncol(tab), sd = 0.1),
              nrow(tab),ncol(tab)),
     pch=19,
     
     col= ifelse(modelS$fitted.values>0.5,
                 rgb(1,0,0,0.5),rgb(0.1,0.1,0.1,0.1)))
hist(modelS$fitted.values)


predTL<- predict(modelL,tab,type = "response")
predTS<- predict(modelS,tab,type = "response")
#predict(modelS,tab,type = "terms")

errorTL<-mean(abs(tab$y-(predTL>0.5)))
contTL<-table(tab$y, predTL>0.5)

errorTS<-mean(abs(tab$y-(predTS>0.5)))
contTS<-table(tab$y, predTS>0.5)

Xval$class<-ifelse(grepl("pre_miRNA", rownames(Xval)), 1, 0)
tabVal<-data.frame(Xval[, ids], y=Xval$class)
colnames(tabVal)<-c(colnames(Xval)[ids], "y")
head(tabVal)

predVL<- predict(modelL,tabVal,type = "response")
predVS<- predict(modelS,tabVal,type = "response")

errorVL<-mean(abs(tabVal$y-(predVL>0.5)))
contVL<-table(tabVal$y, predVL>0.5)

errorVS<-mean(abs(tabVal$y-(predVS>0.5)))
contVS<-table(tabVal$y, predVS>0.5)


errorTL
errorTS
errorVL
errorVS
contTL
contTS
contVL
contVS



pca<-prcomp(tab[,-ncol(tab)],scale=T)
biplot(pca)
plot(pca)
pca$rotation

tabpca<-data.frame(pca$x[,1:3], y=tab$y)

# tabla validacion
tabVpca<-data.frame(predict(pca, tabVal)[,1:3],y=tabVal$y)
head(tabVpca)

cor(tabVpca)
cor(tabpca)

modelpcaTL<-lm(formula=y~., dat=tabpca)
modelpcaTS<-glm(formula=y~., dat=tabpca,
                family = binomial(link="logit"))
summary(modelpcaTL)
summary(modelpcaTS)

# predecir valores

predpcaTL<-predict(modelpcaTL,tabpca,type="response")
predpcaTS<-predict(modelpcaTS,tabpca,type="response")

predpcaVL<-predict(modelpcaTL,tabVpca,type="response")
predpcaVS<-predict(modelpcaTS,tabVpca,type="response")

errorpcaTL<-mean(abs(tabpca$y-(predpcaTL>0.5)))
contpcaTL<-table(tabpca$y, predpcaTL>0.5)

errorpcaTS<-mean(abs(tabpca$y-(predpcaTS>0.5)))
contpcaTS<-table(tabpca$y, predpcaTS>0.5)

errorpcaVL<-mean(abs(tabVpca$y-(predpcaVL>0.5)))
contpcaVL<-table(tabVpca$y, predpcaVL>0.5)

errorpcaVS<-mean(abs(tabVpca$y-(predpcaVS>0.5)))
contpcaVS<-table(tabVpca$y, predpcaVS>0.5)






plot(tabpca[, 1:2], col=ifelse(tabpca$y, "red", "black"),pch=19)
points(tabVpca[,1:2], pch=ifelse(predpcaVL>0.5, 2, 1),
       col=ifelse(tabpca$y, "red", "black"))


library(pROC)
rocvl<-roc(reponse=tabVpca$y, predictor = predpcaVL)

dev.off()










