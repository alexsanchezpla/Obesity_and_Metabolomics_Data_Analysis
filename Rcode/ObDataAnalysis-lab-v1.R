## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.width=6, fig.height=4, fig.path='images/', warning=FALSE)

## ----results='asis', echo=FALSE, include=FALSE, eval=FALSE---------------
## source("https://raw.githubusercontent.com/miriamMota/scripts/master/installifnot.R")
## suppressPackageStartupMessages(installBiocifnot("limma"))
## suppressPackageStartupMessages(installBiocifnot("genefilter"))
## suppressPackageStartupMessages(installifnot("FactoMineR"))
## suppressPackageStartupMessages(installifnot("mixOmics"))
## 

## ----obData--------------------------------------------------------------
load("obesityData.Rda")
dim(obData)
colnames(obData)[c(1:10, 60:70)]

## ------------------------------------------------------------------------
if (! exists("obData")) load("obesityData.Rda")
xClin <- as.matrix(obData[,5:29]); dim(xClin)
xMetab <- as.matrix(obData[,30:275]); dim(xMetab)
YDisPhe <- as.factor(obData[,3]); 
YObes <- as.factor(obData[,1]) # Y is a factor, we chose it as the
YDiab <- as.factor(obData[,2]) # Y is a factor, we chose it as the

## ------------------------------------------------------------------------
scaledClin <-scale(xClin)
round(apply(scaledClin,2,mean),2)
round(apply(scaledClin,2,sd),2)

## ------------------------------------------------------------------------
require(genefilter)
teststat <-rowttests(t(scaledClin), YObes)
head(teststat, n=25)
topDown<-order(teststat$p.value)
ranked<-teststat[topDown,]
print(top10<-ranked[1:10,])

## ------------------------------------------------------------------------
x<-ranked$dm; y<- -log(ranked$p.value)
plot(x, y, xlab="Fold Change", ylab ="-logPval", main="Volcano plot\nA vs B")
text (x[1:10], y[1:10],rownames(ranked)[1:10], cex=0.7)

## ------------------------------------------------------------------------
adjP <- p.adjust(ranked$p.value)
ranked <- cbind(ranked, adjP)
print(ranked[1:20,])


## ----clinical------------------------------------------------------------
require(FactoMineR)
clin.mfa <- MFA(base = obData[,1:29], group=c(4, 25), type=c("n","s"),
                ncp=4,name.group=c("pheno","clin"), graph=FALSE,
                num.group.sup=NULL)

## ------------------------------------------------------------------------
barplot(clin.mfa$eig[,1],main="Eigenvalues",names.arg=1:nrow(clin.mfa$eig))
plot.MFA(x= clin.mfa, choix="axes", habillage="group")
plot.MFA(x= clin.mfa, invisible="quali", choix="ind", habillage="Fen_Group")
plot.MFA(clin.mfa, invisible="ind", partial="all", habillage="Fen_Group")
plot.MFA(clin.mfa, choix="var", habillage="group")
# plotellipses(clin.mfa, habillage="Fen_Group", keepvar = "quali.sup" )

## ------------------------------------------------------------------------
res.mfa <- MFA(base = obData, group=c(4, 25, 246), type=c("n",rep("s",2)),
               ncp=5, name.group=c("pheno","clin","metab"), 
               graph=FALSE, num.group.sup=c(1))

plot.MFA(x= res.mfa, choix="axes", habillage="group")
plot.MFA(x= res.mfa, invisible="quali", choix="ind", habillage="Fen_Group")
plot.MFA(res.mfa, invisible="ind", partial="all", habillage="Fen_Group")
plot.MFA(res.mfa, choix="var", habillage="group",  lim.cos2.var=0.51)
plotellipses(res.mfa, habillage="Fen_Group", keepvar = "quali.sup" )

## ----examples1-----------------------------------------------------------
x<- 3; y<- 1:10; show(x); show(y)
z<- rnorm(100)
hist(z)

