## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----installation from CRAN, eval = FALSE-------------------------------------
#  install.packages("DiffCorr")

## ----installation from GitHub, eval = FALSE-----------------------------------
#  install.packages("devtools")
#  install.packages(c("igraph", "fdrtool"))
#  
#  if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
#  BiocManager::install(c("pcaMethods", "multtest"))
#  
#  library(devtools)
#  install_github("afukushima/DiffCorr")

## ----setup, message = FALSE---------------------------------------------------
library(DiffCorr)

## ----Golub dataset------------------------------------------------------------
golub.df <- read.table("https://coxpress.sourceforge.net/golub.txt", 
                       sep = "\t", header = TRUE, row.names = 1)
dim(golub.df)

## ----clustering---------------------------------------------------------------
hc.mol1 <- cluster.molecule(golub.df[, 1:27], "pearson", "average")  ## ALL (27 samples)
hc.mol2 <- cluster.molecule(golub.df[, 28:38], "pearson", "average") ## AML (11 samples)

## ----cutting tree-------------------------------------------------------------
g1 <- cutree(hc.mol1, h = 0.4)
g2 <- cutree(hc.mol2, h = 0.4)
##
res1 <- get.eigen.molecule(data = golub.df, groups = g1)
res2 <- get.eigen.molecule(data = golub.df, groups = g2)

## ----visualizing modules------------------------------------------------------
gg1 <- get.eigen.molecule.graph(res1)
plot(gg1, layout = layout.fruchterman.reingold(gg1))

gg2 <- get.eigen.molecule.graph(res2)
plot(gg2, layout = layout.fruchterman.reingold(gg2))

## ----writing, eval = FALSE----------------------------------------------------
#  write.modules(g1, res1, outfile = "module1_list.txt")
#  write.modules(g2, res2, outfile = "module2_list.txt")

## ----examination--------------------------------------------------------------
for (i in 1:length(res1$eigen.molecules)) {
  for (j in 1: length(res2$eigen.molecules)) {
    r <- cor(res1$eigen.molecules[[i]],res2$eigen.molecules[[j]], method = "spearman")
    if (abs(r) > 0.8) {
      print(paste("(i, j): ", i, " ", j, sep = ""))
      print(r)
    }
  }
}

cor(res1$eigen.molecules[[2]], res2$eigen.molecules[[8]], method = "spearman")
plot(res1$eigen.molecules[[2]], res2$eigen.molecules[[8]])
plot(res1$eigen.molecules[[21]], res2$eigen.molecules[[24]])

## ----examination of groups of interest graphically----------------------------
plotDiffCorrGroup(golub.df, g1, g2, 21, 24, 1:27, 28:38,
                    scale.center = TRUE, scale.scale = TRUE,
                    ylim=c(-5,5))

## ----export, eval = FALSE-----------------------------------------------------
#  comp.2.cc.fdr(output.file = "res.txt", golub.df[, 1:27], golub.df[, 28:38], threshold = 0.05, save = TRUE)

## ----data---------------------------------------------------------------------
data(AraMetLeaves)
dim(AraMetLeaves)

## ----AraMetLeaves-------------------------------------------------------------
colnames(AraMetLeaves)
?AraMetLeaves

## ----DiffCorr for AraMetLeaves, eval = FALSE----------------------------------
#  comp.2.cc.fdr(output.file = "Met_DiffCorr_res.txt",
#                log10(AraMetLeaves[, 1:17]),   ## Col-0 (17 samples)
#                log10(AraMetLeaves[, 18:37]),  ## tt4 (20 samples)
#                method = "pearson",
#                threshold = 1.0, save = TRUE)
#  

