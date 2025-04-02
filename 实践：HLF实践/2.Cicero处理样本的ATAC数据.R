if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("cicero") #这个感觉不太行啊
BiocManager::install("monocle")
install.packages("devtools")
devtools::install_github('cole-trapnell-lab/monocle3')

install.packages("sf")
install.packages("units")
BiocManager::install("renv")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor'))
devtools::install_github("r-quantities/units")

config <- c(units="--with-udunits2-include=/disk1/cai026/anaconda3/envs/celloracle/bin/ --with-udunits2-lib=/disk1/cai026/anaconda3/envs/celloracle/include/")
install.packages("units",
                 configure.args = config)
options(configure.args=config)
renv::install("units")
library(units)

BiocManager::install(c("Gviz", "GenomicRanges", "rtracklayer"))
install.packages("remotes")
remotes::install_github("cole-trapnell-lab/monocle3")
devtools::install_local("/disk1/cai026/CellOracle/HLF-Celloracle-GSM4829413/monocle3-master")

#Installing 3 packages: units, sf, spdep
library(sf)
#总结一下失败的原因：因为安装monocle3(devtools本地安装)需要用到3个包 units sf spdeq
#但是安装units包报错的是需要 sudo install udunits2-devel 来配置libudunits2.so

library(cicero)
library(monocle3)
library(SingleCellExperiment)
library(VGAM)
library(igraph)

# Create a folder to save results
getwd()
setwd("/disk1/cai026/CellOracle/HLF-Celloracle-GSM4829413/")
output_folder <- "/disk1/cai026/CellOracle/HLF-Celloracle-GSM4829413/cicero_output"
dir.create(output_folder)

####1.加载数据####
#参考教程：https://www.jianshu.com/p/cd0c3c4c11b9
#加载cell-by-peak的count矩阵
# read in matrix data using the Matrix package
indata <- Matrix::readMM("/disk1/cai026/CellOracle/HLF-Celloracle-GSM4829413/data/GSM4829412_cll_atac_filtered_matrix.mtx")
# binarize the matrix
indata@x[indata@x > 0] <- 1

# 加载cell的metadata
# format cell info
cellinfo <- read.table("/disk1/cai026/CellOracle/HLF-Celloracle-GSM4829413/data/GSM4829412_cll_atac_filtered_barcodes.tsv")
row.names(cellinfo) <- cellinfo$V1
names(cellinfo) <- "cells"

# 加载peak的metadata
# format peak info
peakinfo <- read.table("/disk1/cai026/CellOracle/HLF-Celloracle-GSM4829413/data/GSM4829412_cll_atac_peaks.bed")
names(peakinfo) <- c("chr", "bp1", "bp2")
peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
row.names(peakinfo) <- peakinfo$site_name

row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)
head(indata)

# 使用newCellDataSet函数构建CDS对象
# make CDS
fd <- methods::new("AnnotatedDataFrame", data = peakinfo)
pd <- methods::new("AnnotatedDataFrame", data = cellinfo)
input_cds <-  suppressWarnings(newCellDataSet(indata,
                                              phenoData = pd,
                                              featureData = fd,
                                              expressionFamily=VGAM::binomialff(),
                                              lowerDetectionLimit=0))

input_cds@expressionFamily@vfamily <- "binomialff"
input_cds <- monocle::detectGenes(input_cds)

# 数据初步过滤
#Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 

#Cicero将数据存储在CellDataSet（CDS）类的对象中，该类继承自Bioconductor的ExpressionSet类。我们可以使用fData: 获取feature的元信息
# 查看feature的metadata
head(fData(input_cds))
# 查看cell的metadata pData: 获取cell/sample的元信息
head(pData(input_cds))
# 查看cell-by-peak的count矩阵
head(exprs(input_cds))
#对比cicero提供的数据集  并没有问题（cai099）

###接下来完全可以按照celloracle的教程来

####2.Qauality check and Filtering####

# Visualize peak_count_per_cell
hist(Matrix::colSums(exprs(input_cds)))

# Filter cells by peak_count
# Please set an appropriate threshold values according to your data 
max_count <-  15000
min_count <- 2000
input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) >= min_count] 
input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) <= max_count] 

####3.Process Cicero-CDS object####
# Data preprocessing
set.seed(2017)

#参考：https://blog.csdn.net/u012110870/article/details/115511986
#发现下面的detect_genes是monocle3里面的函数
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

###由于monocles包的下载需要units包 再需要sudo yum install udunits2-devel 所以不用教程里面的函数了 直接照cai099的来
#Cicero使用k最近邻方法来创建一个重叠的细胞集，并基于细胞相似度的降维坐标图（如来自tSNE或DDRTree图）来构建这些集合。
set.seed(2017)
input_cds <- detectGenes(input_cds)
input_cds <- estimateSizeFactors(input_cds)
# 使用tSNE方法进行数据降维
input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                             reduction_method = 'tSNE', norm_method = "none")
#到这里报错：NA subscripts in x[i,j] not supported for 'x' inheriting from sparseMatrix
# 检查并处理 NA 值
if (any(is.na(exprs(input_cds)))) {
  exprs(input_cds)[is.na(exprs(input_cds))] <- 0
}

# 如果数据量非常大，先进行 PCA 降维
if (ncol(exprs(input_cds)) > 1000) {
  pca_result <- prcomp(as.matrix(exprs(input_cds)), center = TRUE, scale. = TRUE)
  reduced_data <- predict(pca_result, newdata = as.matrix(exprs(input_cds)))[, 1:6]  # 取前6个主成分
  exprs(input_cds) <- reduced_data
}

# 确保数据不是稀疏矩阵
if (inherits(exprs(input_cds), "sparseMatrix")) {
  exprs(input_cds) <- as.matrix(exprs(input_cds))
}
# 增加内存限制
#memory.limit(size = 8000)  # 设置内存限制为 8 GB

# 使用tSNE方法进行数据降维
input_cds <- reduceDimension(input_cds, max_components = 2, num_dim=6,
                             reduction_method = 'tSNE', norm_method = "none")


# 提取tSNE降维后的坐标
tsne_coords <- t(reducedDimA(input_cds))
row.names(tsne_coords) <- row.names(pData(input_cds))
# 使用make_cicero_cds函数构建cicero CDS对象=CDS对象+降维后的坐标，降维后的坐标reduce_coordinates应该采用data.frame或矩阵的形式，其中行名称与CDS的pData表中的细胞ID相匹配，列为降维对象的坐标
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = tsne_coords)


####4.运行Cicero包####
#Cicero包的主要功能是通过计算评估全基因组中位点之间的共可及性（成对的），从而预测顺式调控相互作用。
#有两种方法：1.直接函数一步出结果（推荐） 2.分步使用多个函数，较灵活  下面用第一种
# 加载内置的人hg19参考基因组的坐标信息
data("human.hg19.genome")
chromosome_length <- mhuman.hg19.genome

# Run the main function
conns <- run_cicero(cicero_cds, chromosome_length) # Takes a few minutes to run
# 查看运行的结果
head(conns)

####5. Save results for the next step####
all_peaks <- row.names(exprs(input_cds))
write.csv(x = all_peaks, file = paste0(output_folder, "/all_peaks.csv"))
write.csv(x = conns, file = paste0(output_folder, "/cicero_connections.csv"))


