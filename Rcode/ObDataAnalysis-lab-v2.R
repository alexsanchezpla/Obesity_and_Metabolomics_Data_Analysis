## ----global_options, include=FALSE---------------------------------------
knitr::opts_chunk$set(fig.width=6, fig.height=4, fig.path='images/', warning=FALSE)

## ----results='asis', echo=FALSE, include=FALSE, eval=FALSE---------------
## source("https://raw.githubusercontent.com/miriamMota/scripts/master/installifnot.R")
## installifnot<- function(pckgName){
##   if (!(require(pckgName, character.only = TRUE))) {
##     install.packages(pckgName, dep = TRUE)
##     require(pckgName)
##   }
## }
## suppressPackageStartupMessages(installifnot("FactoMineR"))
## suppressPackageStartupMessages(installifnot("mixOmics"))
## suppressPackageStartupMessages(installBiocifnot("limma"))

## ----obData--------------------------------------------------------------
load("obesityData.Rda")
dim(obData)
colnames(obData)[c(1:10, 60:70)]

## ----data4classif--------------------------------------------------------
if (! exists("obData")) load("obesityData.Rda")
xClin <- as.matrix(obData[,5:29]); dim(xClin)
xMetab <- as.matrix(obData[,30:275]); dim(xMetab)
YDisPhe <- as.factor(obData[,3]); 
YObes <- as.factor(obData[,1]) # Y is a factor, we chose it as the
YDiab <- as.factor(obData[,2]) # Y is a factor, we chose it as the
(tabDisPhe<- table(YObes, YDiab))
margin.table((tabDisPhe),1)
margin.table((tabDisPhe),2)

## ------------------------------------------------------------------------
require(mixOmics)
resClinDisPhe <- splsda(xClin, YDisPhe, ncomp = 3, keepX = c(10, 10, 10)) 
# where keepX is the number of variables selected for each components
# where ncomp is the number of components wanted

## ------------------------------------------------------------------------
plotIndiv(resClinDisPhe, ind.names = YDisPhe, add.legend = TRUE, plot.ellipse =TRUE)

## ------------------------------------------------------------------------
selectVar(resClinDisPhe, comp = 1)
selectVar(resClinDisPhe, comp = 2)

## ------------------------------------------------------------------------
set.seed(45)
errClinDisPhe <- perf(resClinDisPhe, validation = "Mfold", folds = 8,dist = "centroids", progressBar = FALSE)
errClinDisPhe$error.rate
errClinDisPhe$error.rate.class
plot(errClinDisPhe, type = "l")

## ------------------------------------------------------------------------
X <- as.matrix(xClin)
Y <- YDisPhe
i <- 1
samp <- sample(1:3, nrow(X), replace = TRUE) 

test <- which(samp == i) # Search which column in samp has a value of 1
train <- setdiff(1:nrow(X), test) # Keeping the column that are not in test

splsda.train <- splsda(X[train, ], Y[train], ncomp = 3, keepX = c(10, 10, 10))
newClass <- predict(splsda.train, X[test, ], method = "max.dist")
# newClass$class
(Prediction1 <- newClass$class$max.dist[, 1])
(Prediction2 <- newClass$class$max.dist[, 2])
Original <- as.character(YDisPhe[test])
cbind(Original, Prediction1, Prediction2)

## ----data4CCA------------------------------------------------------------
if (! exists("obData")) load("obesityData.Rda")
xClin <- as.matrix(obData[,5:29]); dim(xClin)
xMetab <- as.matrix(obData[,30:275]); dim(xMetab)
YDisPhe <- as.factor(obData[,3]); 
YObes <- as.factor(obData[,1]) # Y is a factor, we chose it as the
YDiab <- as.factor(obData[,2]) # Y is a factor, we chose it as the

## ------------------------------------------------------------------------
X <- as.matrix(xClin)
Y <- as.matrix(xMetab)
icSep<- imgCor(X[1:4,1:3], Y[1:4, 1:3], type="separate", interactive=FALSE)

## ------------------------------------------------------------------------
icSep<- imgCor(X, Y, type="separate", interactive=FALSE)

## ----paramTuning, eval=FALSE---------------------------------------------
## grid1 <- seq(0.01, 3, length = 30)
## grid2 <- seq(0.01, 0.1, length = 50)
## ## Depending on how fine your grid is, estim.regul may take some time
## cv.score <- tune.rcc(X, Y, grid1=grid1, grid2=grid2)
## # This call provides the following values
## # lambda1 =  0.205
## # lambda2 =  0.05573469
## # CV-score =  0.5427141
## save(cv.score, file= "cv.score.Rdata")

## ----runCCA, eval=FALSE--------------------------------------------------
## ## Run rCCA given the optimal parameters:
## lambda1 =  0.205 # cv.score$opt.lambda1
## lambda2 =  0.05573469 # cv.score$opt.lambda2
## result <- rcc(X, Y, ncomp = 3, lambda1 = lambda1, lambda2 = lambda2)

## ----runShrinkedCCA------------------------------------------------------
## Run rCCA using the Shrinkage method:
result <- rcc(X, Y, ncomp = 3, method="shrinkage")

## ------------------------------------------------------------------------
plot(result, scree.type = "barplot")

## ----samplesPlot, eval=FALSE---------------------------------------------
## indNames <- rownames(X)
## groups <- YDisPhe
## col.groups <-as.numeric(groups)
## plotIndiv(result, comp = 1:2, col = col.groups, cex=0.8)

## ----plotVar, eval=FALSE-------------------------------------------------
## plotVar(result, comp = 1:2, cutoff=0.5, cex = c(0.8, 0.8), rad.in = 0.1)
## # plotVar(result, comp = 1:2, cutoff=0.90, X.label=TRUE, Y.label=TRUE, cex = c(0.8, 0.8))
## # plot3dVar(result, cutoff = 0.6, cex = c(1.5, 1.5), label.axes.box = "axes")

## ----relNetwork, eval=FALSE----------------------------------------------
## # The user can set a correlation threshold to only represent variables above that threshold
## network(result, comp = 1:3, interactive = FALSE, cutoff = 0.5)

## ----examples1-----------------------------------------------------------
x<- 3; y<- 1:10; show(x); show(y)
z<- rnorm(100)
hist(z)

