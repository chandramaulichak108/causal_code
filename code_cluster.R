## Unzip the zip file
unzip(zipfile = 'cupc-20230605T042950Z-001.zip', exdir = '.')

## Shell function
shell_call <- function(command, ...) {
  result <- system(command, intern = TRUE, ...)
  cat(paste0(result, collapse = "\n"))
}

shell_call("ls cupc") ## Just a check
cat("\n")
args <- commandArgs(trailingOnly = TRUE)
#print(as.numeric(args))
param_d = as.numeric(args[1])
alpha_d = as.numeric(args[2])
class_meta = args[-c(1,2)]
print(args)
print(param_d)
print(alpha_d)
print(class_meta)
setwd('cupc')

shell_call("nvcc -O3 --shared -Xcompiler -fPIC -o Skeleton.so cuPC-S.cu") ## Required

print("required done!")

# install.packages('pak')

# library(pak)
# pak::pkg_install(c("graph", "RBGL", "Rgraphviz","pcalg","Seurat","tictoc"))

library(pcalg)
library(graph)
library(MASS)
library(tictoc)
library(igraph)
library(Seurat)

meta = read.csv("~/Causal/Code for cluster/CP/meta.tsv",sep="\t")
# meta = read.csv("~/Causal/Code for cluster/Healthy/meta.tsv",sep="\t")
print(levels(factor(meta$Cluster)))

library(Matrix)

gene_filter = function(x,pct,cluster){
  counts <- GetAssayData(x)[,meta$Cluster %in% cluster]
  #   counts <- GetAssayData(x)
  print(dim(counts))

  genes.percent.expression <- rowMeans(counts>0 )*100   
  genes.filter <- names(genes.percent.expression[genes.percent.expression>pct])  #select genes expressed in at least 1% of cells
  counts.sub <- counts[genes.filter,]
  return(counts.sub)
}

#gene_filter = function(x,pct){
#  # counts <- GetAssayData(x)[,meta$Cluster %in% cluster]
#  counts <- GetAssayData(x)
#  print(dim(counts))

#  genes.percent.expression <- rowMeans(counts>0 )*100
#  genes.filter <- names(genes.percent.expression[genes.percent.expression>pct])  #select genes expressed in at least 1% of cells
#  counts.sub <- counts[genes.filter,]
#  return(counts.sub)
#}

data_d = readRDS("~/Causal/Code for cluster/CP/chronic_pancreatitis.rds")
#data_d = readRDS("~/Causal/Code for cluster/Healthy/adult_pancreas.rds")

data = gene_filter(data_d,param_d,class_meta)

#data = gene_filter(data_d,15)

print(dim(data))

source("cuPC.R")

############# faltu
n_btp = 50
pval = 0.05

# data_new = t(data)[sample(1:dim(data)[2],size = 2000),]
data_new = t(data)
# data_new = round(data_new)

library(pcalg)

V <- colnames(data_new)
V = 1:dim(data_new)[2]

V = as.character(V)

## define sufficient statistics

# suffStat <- list(dm = as.matrix(data_new), nlev = rep(10,dim(data_new)[2]), adaptDF = FALSE) #nlev has number of levels in each variable in datamatrix
suffStat <- list(dm = round(as.matrix(data_new)), nlev = rep(10,dim(data_new)[2]), adaptDF = FALSE)


## estimate CPDAG
start_time <- Sys.time()
#pc.D <- pc(suffStat, indepTest = disCItest, alpha = 0.05, labels = V, verbose = TRUE)

tic()
#cuPC_fit <- cu_pc(suffStat, indepTest = disCItest, alpha = alpha_d, labels = V, verbose = TRUE)
cuPC_fit <- pc(suffStat, indepTest = disCItest, alpha = alpha_d, labels = V, verbose = TRUE,numCores = 20)
print("The total time consumed by cuPC is:")
toc()
print(cuPC_fit)
# pc.D.cp <- pc(suffStat, indepTest = disCItest, alpha = 0.2, labels = V, verbose = TRUE,numCores = 5)
end_time <- Sys.time()
###################


# save.image(file = paste("~/Causal/Code for cluster/CP/diseased_Activated_Stellate", param,".RData"))
# save.image(file = "~/Causal/Code for cluster/CP/diseased_Activated_Stellate20.RData")
# save.image(file = "~/Causal/Code for cluster/Healthy/healthy_Activated_Stellate_Quiescent_Stellate15.RData")
#save.image(file = "~/Causal/Code for cluster/Healthy/healthy15.RData")
#save.image(file = "~/Causal/Code for cluster/CP/diseased_Acinar-Reg+20.RData")
#save.image(file = "~/Causal/Code for cluster/CP/diseased_Acinar-REG+40.RData")
save.image(file = paste("~/Causal/Code for cluster/CP/diseased_",class_meta,param_d,"_",alpha_d,".RData",sep=""))
print("done!")



