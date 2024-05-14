
# Modelo de predicción de pre miRNAs basado en sus estructuras
# Mediante la anogtación de ".","(" y ")".

# Autor= Rodrigo de Jesús Muñoz
# Fecha= 09/07/2019


rm(list=ls())

# Insertar la tabla 

# x<-read.table("X_Local_ss_sRNA.tab", sep="\t",header=T,
#                       row.names=1)

load("X_Local_ss_sRNA.RData")
x<-as.data.frame(X)
set.seed(7)
xval<-x[sample(row.names(x),nrow(x)/2),]
x<-x[setdiff(row.names(x),row.names(xval)),]

pdf(file = "Modelo_estr_sec-copia")


# calcular correlacio
cor(x$.........,x$`)))))))))`)
plot(x$.........,x$`)))))))))`)

cormat<-cor(x)
cormat[is.na(cormat)]<-0
heatmap(cormat,col=colorRampPalette(c("blue","white","red"))(50))

rownames(x)
class(x)
x$clase<- ifelse(grepl("pre_miRNA", rownames(x)),1,0)
table(x$clase)


# cor(x$clase, x$X.........)

xcor<-NULL
for (i in seq(1,200) ){
  print(i)
  
xcor[i]<-  cor(x$clase, x[,i], method = "spearman")
}
xcor[is.na(xcor)]<-0

range(xcor)
summary(xcor)
hist(xcor)

xcor[xcor == max(xcor)]
idmax<-which.max(xcor)
ids<-xcor>=0.35
table(ids)

x[,c(colnames(X)[ids], "clase")]
head(x[,c(colnames(X)[ids], "clase")])


# GEnerar tabla

tab<-data.frame(X[,ids], x$clase)
colnames(tab)<- c(colnames(X[,ids]),"y")
head(tab)
# tab<-tab[order(tab$x),] 

summary(tab)
cor(tab)

plot(((tab)))
grid()

modelL<-lm(formula = y~.,data=tab)
summary(modelL)
plot(tab+
       matrix(rnorm(nrow(tab)*ncol(tab), sd = 0.1),
             nrow(tab),ncol(tab)),
     col=ifelse(modelL$fitted.values>0.5,
                rgb(1,0,0,0.5),rgb(0.1,0.1,0.1,0.1)))
hist(modelL$fitted.values)

# abline(h=0.5)

modelS<- glm(formula=y~., data = tab, family = binomial(link = "logit"))
summary(modelS)

plot(tab+
       matrix(rnorm(nrow(tab)*ncol(tab), sd = 0.1),
              nrow(tab),ncol(tab)),
     col=ifelse(modelS$fitted.values>0.5,
                rgb(1,0,0,0.5),rgb(0.1,0.1,0.1,0.1)),pch=19)
hist(modelS$fitted.values)


predTL<-predict(modelL,tab,type = "response")
predTS<-predict(modelS,tab,type = "response")

# predict(modelS,tab,type = "terms")

errorTl<-mean(abs(tab$y-(predTS>0.5)))
contTL<-table(tab$y,predTS>0.5)

errorTS<-mean(abs(tab$y-(predTL>0.5)))
contTS<-table(tab$y,predTL>0.5)


xval$class<-ifelse(grepl("pre_miRNA",rownames(xval)), 1, 0)
tabVal<-data.frame(xval[, ids], y=xval$class)
colnames(tabVal)<-c(colnames(xval)[ids], "y")

predVL<-predict(modelL,tab,type = "response")
predVS<-predict(modelS,tab,type = "response")

errorVL<-mean(abs(tabVal$y-(predVL>0.5)))
contVL<-table(tabVal$y,predVL>0.5)

errorVS<-mean(abs(tabVal$y-(predVS>0.5)))
contVS<-table(tabVal$y,predVS>0.5)

errorTl
errorTS
errorVL
errorVLS

contTL
contTS
contVL
contVS


dev.off()

