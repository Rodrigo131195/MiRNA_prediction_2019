

# Remove list

rm(list = ls())

Tabla_row_fix<- read.table("Tabla_row_fix.txt", sep = "\t",
                           header = TRUE, row.names = 1)



# Crear vectores

z <- Tabla_row_fix$dG
x <- as.numeric(z)

w <- ifelse(Tabla_row_fix$type=="lnc_RNA",
  1,ifelse(Tabla_row_fix$type=="miRNA",
  2,ifelse(Tabla_row_fix$type=="pre_miRNA",
  3,ifelse(Tabla_row_fix$type=="tRNA",
  4,ifelse(Tabla_row_fix$type=="rRNA",
  5,ifelse(Tabla_row_fix$type=="snRNA",6,7))))))

y <- as.numeric(w)

a<- (y==1)
lnc_RNA<- x[a]
b<- (y==2)
miRNA<- x[b]
c<- (y==3)
pre_miRNA<- x[c]
d<- (y==4)
tRNA<- x[d]
e<- (y==5)
rRNA<- x[e]
f<- (y==6)
snRNA<- x[f]
g<- (y==7)
snoRNA<- x[g]



plot(x,y)
boxplot(list(snRNA,miRNA, rRNA, pre_miRNA, lnc_RNA,tRNA,snoRNA),
        names = c("snRNA","miRNA", "rRNA", "pre_miRNA", "lnc_RNA","tRNA","snoRNA"))
 
allRNA<- Tabla_row_fix$dG

miRNAd<- ifelse(Tabla_row_fix$type=="pre_miRNA", 1,0)


# Crear tabla
tab<- data.frame(x=allRNA, y = miRNAd)
tab<- tab[order(tab$x),]

# Modelo lineal

modeloL<- lm(formula=y~x, data=tab)
summary(modeloL)

# # Modelo sigmoide
# 
# sigmoide <- function(x = NULL, y = NULL, b = 1, a = 0) {
#   
#   # Cacula a partir de la definiciÃ³n de la sigmoide
#   s <- 1 / (1 + exp(-(b*(x-a))))
#   
#   return(s)
# }
# 
# sigmoide(x = 0)
# x0 <- seq(-3, 3, 0.1)
# plot(x0, sigmoide(x0, a = 0, b = 10), type = "l", lwd = 3)
# abline(h = 0.5, v = 0, col = rgb(0.1, 0.1, 0.1 ,0.5), 
#        lwd = 2, lty = 2) 
# grid()

 
modeloS <- glm(formula = y~x, data = tab, 
               family = binomial(link = "logit"))


# Gráfica

cols <- ifelse(tab$y == 1, "red", "blue")
plot(tab, xlab="Energía",ylab="miRNA o no miRNA", col=cols )
grid()
abline(h=0.5)
lines(tab$x,modeloL$fitted.values, col="black")
lines(tab$x,modeloS$fitted.values, col="orange")


predict(modeloS, data.frame(x=-70),type="response")









