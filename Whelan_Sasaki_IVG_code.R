
library(dplyr) #version 1.1.4
library(Seurat) #version 5.1.0
library(SeuratWrappers) #version 0.3.1
library(patchwork) #version 1.2.0
library(gplots) #version 3.1.3.1
library(ggplot2) #version 3.5.1
library(SeuratObject) #version 5.0.2
library(monocle3) #version 1.3.1
library(cowplot) #version 1.1.3
library(viridis) #version 0.6.5
library(tidyverse) #version 2.0.0
library(forcats) #version 1.0.0
library(scCustomize) #version 2.1.2
library(stringr) #version 1.5.1
library(ComplexHeatmap) #version 2.14.0
library(circlize) #version 0.4.16
library(lattice) #version 0.22-6
library(scales) #version 1.3.0
library(ggridges) #version 0.5.6
library(reshape2) #version 1.4.4
library(ggrastr) #version 1.0.2
library(tibble) #version 3.2.1
library(scico) #version 1.5.0
library(RColorBrewer) #version 1.1-3
library(DoubletFinder) #version 2.0.3

load("HS_reclustered_CC_integrated.Robj") #all human with cell cycle removal
load("NCG_germ4_2023.07.23.Robj") #this is the final in vitro xrTestis object
load("HS_NCG_2023.07.19.Robj") #this is the combined in vivo and in vitro object
load("NCG_germ_macaque4_2024.01.03.Robj") #this is the macaque xrTestis object
load("spermatogonia_comparison2.Robj") #this is the subset of spermatogonia


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### PART 1          ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### Human xrTestis  ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


setwd("/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/")


#Read in the 10X data. Here samples have been aligned against both human and mouse separately in order to remove mouse cells later.

#d123
NCG2_Hu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_2_Hu/filtered_feature_bc_matrix")
NCG2_Mu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_2_Mu/raw_feature_bc_matrix")

#d127
NCG3_Hu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_3_Hu/filtered_feature_bc_matrix")
NCG3_Mu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_3_Mu/raw_feature_bc_matrix")

#d180
NCG4_Hu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_4_Hu/filtered_feature_bc_matrix")
NCG4_Mu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_4-Mu/raw_feature_bc_matrix")

#d181
NCG5_Hu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_5_Hu/filtered_feature_bc_matrix")
NCG5_Mu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_5_Mu/raw_feature_bc_matrix")

#d251
NCG8_Hu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_8_Hu/outs/filtered_feature_bc_matrix")
NCG8_Mu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_8_Mu/NCG_8_raw_feature_bc_matrix")

#d245
NCG9_Hu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_9_Hu/filtered_feature_bc_matrix")
NCG9_Mu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_9_Mu/raw_feature_bc_matrix")

#c91+180
NCG19_Hu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_19_Hu/filtered_feature_bc_matrix")
NCG19_Mu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_19-Mu/raw_feature_bc_matrix")

#c41+d120
NCG_c41_Hu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_c41_Hu/filtered_feature_bc_matrix")
NCG_c41_Mu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_c41_Mu/raw_feature_bc_matrix")

#d120
NCG_182_2_Hu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_182_2-HU/filtered_feature_bc_matrix")
NCG_182_2_Mu <- Read10X(data.dir = "~/Documents/PROJECTS/ksasaki/scRNAseq/NCG_182_2-Mu/raw_feature_bc_matrix")

#ips
C1_191121_Hu <- Read10X_h5("~/Documents/PROJECTS/ksasaki/scRNAseq/191121-C1-Hu/filtered_feature_bc_matrix.h5")
C1_191121_Mu <- Read10X_h5("~/Documents/PROJECTS/ksasaki/scRNAseq/191121-C1-Mu/raw_feature_bc_matrix.h5")

#human testis
HS12H_Hu <- Read10X(data.dir = "/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/Hs12H/filtered_feature_bc_matrix")
HS12H_Mu <- Read10X(data.dir = "/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/HS12_gonad_mouse/raw_feature_bc_matrix")



human.list <- c(NCG2_Hu, NCG3_Hu, NCG4_Hu,
                NCG5_Hu, NCG8_Hu,
                NCG9_Hu, NCG19_Hu, NCG_c41_Hu, NCG_182_2_Hu,
                C1_191121_Hu, HS12H_Hu)
mouse.list <- c(NCG2_Mu, NCG3_Mu, NCG4_Mu, 
                NCG5_Mu, NCG8_Mu, 
                NCG9_Mu, NCG19_Mu, NCG_c41_Mu, NCG_182_2_Mu,
                C1_191121_Mu, HS12H_Mu)
names.data <- c("NCG2", "NCG3", "NCG4", 
                "NCG5", "NCG8", 
                "NCG9", "NCG19", "NCG_c41", "NCG_182_2",
                "C1_191121", "HS12H")
human.list.original <- human.list

mouse.or.human.list <- human.list


#This script iterates through every sample and assigns a "human", "mouse" or "unknown" identity to each cell based on whether human or mouse counts were higher.
for(current.dataset in 1:length(human.list)){
  # current.dataset <- 6
  current.dataset.hu <- human.list[[current.dataset]]
  current.dataset.mu <- mouse.list[[current.dataset]]
  colnames.hu <- colnames(current.dataset.hu)
  # colnames.mu <- colnames(current.dataset.mu)
  # intersect(colnames.hu, colnames.mu)
  length(colnames.hu)
  current.dataset.mu.filtered <- current.dataset.mu[,colnames.hu]
  colnames.mu <- colnames(current.dataset.mu.filtered)
  length(colnames.mu)
  number.of.cells <- length(intersect(colnames.hu, colnames.mu))
  if(number.of.cells != length(colnames.hu)){print("ERROR!!!")}
  counts.hu <- colSums(current.dataset.hu)
  counts.mu <- colSums(current.dataset.mu.filtered)
  human.or.mouse <- counts.hu
  for(current.cell in 1:(length(human.or.mouse))){
      # current.cell <- 1
    # if((counts.hu[current.cell])>6000){
      if((counts.hu[current.cell]/counts.mu[current.cell])>2){
        human.or.mouse[current.cell] <- "human"
      } else {
        if((counts.hu[current.cell]/counts.mu[current.cell])<0.5){
          human.or.mouse[current.cell] <- "mouse"
        } else {
          human.or.mouse[current.cell] <- "unknown"
        }
      }
    # }else{human.or.mouse[current.cell] <- "lowQ"}
  }
  print(names.data[current.dataset])
  print(table(human.or.mouse))
  mouse.or.human.list[[current.dataset]] <- human.or.mouse

}

#imlc
C3_191203_hu <- Read10X(data.dir = "/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/191203-C3/filtered_feature_bc_matrix")
C1_191213_hu <- Read10X(data.dir = "/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/191213-C1/filtered_feature_bc_matrix")



human.list2 <- c(NCG2_Hu, NCG3_Hu, NCG4_Hu,
                NCG5_Hu, NCG8_Hu,
                NCG9_Hu, NCG19_Hu, NCG_c41_Hu, NCG_182_2_Hu,
                C1_191121_Hu, C3_191203_hu, C1_191213_hu)
human.list <- human.list2

for(a in 1:length(human.list)){
  print(length(Cells(human.list[[a]])))
}

inflection.point1 <- c(1:length(human.list))
log_lib_size_at_inflection1 <- c(1:length(human.list))

par(mfrow = c(3,3))

#this code is meant to calculate an inflexion point for each sample. That doesn't work very well but the graphs are informative to manually set minimum UMI values.

for(y in 1:length(human.list)){
  #y=1
  print(names.data[y])
  expression.df <- as.data.frame(human.list[y])

  umi_per_barcode <- colSums(expression.df)
  barcode_rank <- rank(-umi_per_barcode)
  # plot(barcode_rank, umi_per_barcode,
  #      #ylim=c(1,2000),
  #      xlim=c(1,10000))

  log_lib_size <- log10(umi_per_barcode)
  #plot(barcode_rank, log_lib_size, xlim=c(1,10000))

  o <- order(barcode_rank)
  log_lib_size <- log_lib_size[o]
  barcode_rank <- barcode_rank[o]

  rawdiff <- diff(log_lib_size)/diff(barcode_rank)
  inflection <- which(rawdiff == min(rawdiff[500:4000], na.rm=TRUE))



  plot(barcode_rank, log_lib_size, xlim=c(1,8000),
       pch = 20,
       main = paste(y, names.data[y], "cells", inflection, sep="_", "UMIs", round(10^log_lib_size[inflection]))
  )

  abline(v=inflection, col="blue", lwd=2)
  abline(h=log_lib_size[inflection], col="green", lwd=2)


  inflection.point1[y] <- inflection
  log_lib_size_at_inflection1[y] <- log_lib_size[inflection]
}

min.UMIs1 <- 10^log_lib_size_at_inflection1 #something is wrong with this, just make it manually

min.UMIs1 <- vector("numeric" , length(human.list)+2)


#these are the minimum UMI values by sample
min.UMIs1[1] <- 10^4 #NCG2
min.UMIs1[2] <- 10^4 #NCG3
min.UMIs1[3] <- 10^3.8 #NCG4
min.UMIs1[4] <- 10^3.7 #NCG5
min.UMIs1[5] <- 10^3.5 #NCG8
min.UMIs1[6] <- 10^3.8 #NCG9
min.UMIs1[7] <- 10^3.7 #NCG19
min.UMIs1[8] <- 10^3.7 #NCG_c41
min.UMIs1[9] <- 10^3.8 #NCG_182_2
min.UMIs1[10] <- 10^4 #C1191121
min.UMIs1[11] <- 10^4 #HS12H


# These are global minimum and maximum gene and % mitochondria values:
min.features <- 100
max.features <- 3000
max.mito <- 20

# #RENAME THE CELL NAMES TO MATCH VELOCITY 
# 
# cell.names.NCG2_Hu <- gsub("-1", "", colnames(NCG2_Hu))
# cell.names.NCG2_Hu <- paste("NCG_2_Hu:", cell.names.NCG2_Hu, "x", sep = "")
# cell.names.NCG2_Hu -> colnames(NCG2_Hu)
# 
# cell.names.NCG3_Hu <- gsub("-1", "", colnames(NCG3_Hu))
# cell.names.NCG3_Hu <- paste("NCG_3_Hu:", cell.names.NCG3_Hu, "x", sep = "")
# cell.names.NCG3_Hu -> colnames(NCG3_Hu)
# 
# cell.names.NCG4_Hu <- gsub("-1", "", colnames(NCG4_Hu))
# cell.names.NCG4_Hu <- paste("NCG_4_Hu:", cell.names.NCG4_Hu, "x", sep = "")
# cell.names.NCG4_Hu -> colnames(NCG4_Hu)
# 
# cell.names.NCG5_Hu <- gsub("-1", "", colnames(NCG5_Hu))
# cell.names.NCG5_Hu <- paste("NCG_5_Hu:", cell.names.NCG5_Hu, "x", sep = "")
# cell.names.NCG5_Hu -> colnames(NCG5_Hu)
# 
# cell.names.NCG8_Hu <- gsub("-1", "", colnames(NCG8_Hu))
# cell.names.NCG8_Hu <- paste("NCG_8_Hu:", cell.names.NCG8_Hu, "x", sep = "")
# cell.names.NCG8_Hu -> colnames(NCG8_Hu)
# 
# cell.names.NCG9_Hu <- gsub("-1", "", colnames(NCG9_Hu))
# cell.names.NCG9_Hu <- paste("NCG_9_hu:", cell.names.NCG9_Hu, "x", sep = "")
# cell.names.NCG9_Hu -> colnames(NCG9_Hu)
# 
# cell.names.NCG19_Hu <- gsub("-1", "", colnames(NCG19_Hu))
# cell.names.NCG19_Hu <- paste("NCG_19_Hu:", cell.names.NCG19_Hu, "x", sep = "")
# cell.names.NCG19_Hu -> colnames(NCG19_Hu)
# 
# cell.names.NCG_c41_Hu <- gsub("-1", "", colnames(NCG_c41_Hu))
# cell.names.NCG_c41_Hu <- paste("NCG_c41_1_Hu:", cell.names.NCG_c41_Hu, "x", sep = "")
# cell.names.NCG_c41_Hu -> colnames(NCG_c41_Hu)
# 
# cell.names.NCG_182_2_Hu <- gsub("-1", "", colnames(NCG_182_2_Hu))
# cell.names.NCG_182_2_Hu <- paste("NCG_182_2-HU:", cell.names.NCG_182_2_Hu, "x", sep = "")
# cell.names.NCG_182_2_Hu -> colnames(NCG_182_2_Hu)
# 
# cell.names.C1_191121_Hu <- gsub("-1", "", colnames(C1_191121_Hu))
# cell.names.C1_191121_Hu <- paste("191121-C1-Hu:", cell.names.C1_191121_Hu, "x", sep = "")
# cell.names.C1_191121_Hu -> colnames(C1_191121_Hu)
# 
# cell.names.C3_191203_hu <- gsub("-1", "", colnames(C3_191203_hu))
# cell.names.C3_191203_hu <- paste("191203-C3_trimmed:", cell.names.C3_191203_hu, "x", sep = "")
# cell.names.C3_191203_hu -> colnames(C3_191203_hu)
# 
# cell.names.C1_191213_hu <- gsub("-1", "", colnames(C1_191213_hu))
# cell.names.C1_191213_hu <- paste("191213-C1_trimmed:", cell.names.C1_191213_hu, "x", sep = "")
# cell.names.C1_191213_hu -> colnames(C1_191213_hu)

#make seurats

NCG2 <- CreateSeuratObject(NCG2_Hu) # project = "NCG2")
NCG2[["orig.ident"]] <- "NCG2"
NCG3 <- CreateSeuratObject(NCG3_Hu, project = "NCG3")
NCG3[["orig.ident"]] <- "NCG3"
NCG4 <- CreateSeuratObject(NCG4_Hu, project = "NCG4")
NCG4[["orig.ident"]] <- "NCG4"
NCG5 <- CreateSeuratObject(NCG5_Hu, project = "NCG5")
NCG5[["orig.ident"]] <- "NCG5"
NCG8 <- CreateSeuratObject(NCG8_Hu, project = "NCG8")
NCG8[["orig.ident"]] <- "NCG8"
NCG9 <- CreateSeuratObject(NCG9_Hu, project = "NCG9")
NCG9[["orig.ident"]] <- "NCG9"
NCG19 <- CreateSeuratObject(NCG19_Hu, project = "NCG19")
NCG19[["orig.ident"]] <- "NCG19"
NCG_c41 <- CreateSeuratObject(NCG_c41_Hu, project = "NCG_c41")
NCG_c41[["orig.ident"]] <- "NCG_c41"
NCG_182_2 <- CreateSeuratObject(NCG_182_2_Hu, project = "NCG_182_2")
NCG_182_2[["orig.ident"]] <- "NCG_182_2"
# NCG_18_2_6M_1130 <- CreateSeuratObject(NCG_18_2_6M_1130_Hu, project = "NCG_18_2_6M_1130")
# NCG_18_2_6M_1123 <- CreateSeuratObject(NCG_18_2_6M_1123_Hu, project = "NCG_18_2_6M_1123")
C1_191121 <- CreateSeuratObject(C1_191121_Hu, project = "C1_191121")
C1_191121[["orig.ident"]] <- "C1_191121"
C3_191203 <- CreateSeuratObject(C3_191203_hu, project = "C3_191203")
C3_191203[["orig.ident"]] <- "C3_191203"
C1_191213 <- CreateSeuratObject(C1_191213_hu, project = "C1_191213")
C1_191213[["orig.ident"]] <- "C1_191213"


list.of.seurats <- c(NCG2, NCG3, NCG4, NCG5, NCG8, NCG9, NCG19, NCG_c41, NCG_182_2,
                     C1_191121, C3_191203, C1_191213)

experiment <- c("transplant", "transplant", "transplant", "transplant", "transplant", "transplant", "transplant", "transplant", "transplant",
                "hiPSCs", "hPGCLC", "iMELCs")

stage.data <- c("d123", "d127", "d180",
                "d181", "d251",
                "d245", "c91+d180", "c41+120", "d120",
                "hiPSCs", "hPGCLCs", "iMELCs")

names.data <- c("NCG2", "NCG3", "NCG4", 
                "NCG5", "NCG8", 
                "NCG9", "NCG19", "NCG_c41", "NCG_182_2",
                "C1_191121", "C3_191203", "C1_191213")

for(current.sample in 1:length(list.of.seurats)){
  # current.sample <- 1
  print(current.sample)
  current.seurat <- list.of.seurats[[current.sample]]
  if(current.sample <= (length(human.list.original))-1){
    current.mouse.or.human <- mouse.or.human.list[[current.sample]]
    current.seurat$species <- current.mouse.or.human
  }else{
    current.seurat$species <- "human"
  }
  current.seurat$replicate <- names.data[current.sample]
  current.seurat$experiment <- experiment[current.sample]
  # current.seurat$experiment <- "human in vivo"
  current.seurat$stage <- stage.data[current.sample]
  current.seurat[['percent.mito']] <- PercentageFeatureSet(current.seurat, pattern = "^MT-")
  list.of.seurats[[current.sample]] <- current.seurat
}

for(current.sample in 1:(length(list.of.seurats))){
  # current.sample <- 1
  current.seurat <- list.of.seurats[[current.sample]]
  current.seurat <- subset(x = current.seurat,
                           # subset = nFeature_RNA > min.features & nFeature_RNA < max.features & percent.mito < max.mito & nCount_RNA > min.UMIs[current.sample])
                           subset = percent.mito < max.mito & nCount_RNA > min.UMIs1[current.sample]) #
  # current.seurat <- subset(x = current.seurat, subset = species == "human")
  list.of.seurats[[current.sample]] <- current.seurat
}



number.of.samples <- length(list.of.seurats)
number.of.cells.per.sample <-  data.frame(
  SampleName = c(1:number.of.samples),
  CellNumber = c(1:number.of.samples),
  GeneNumber = c(1:number.of.samples),
  SampleNumber = c(1:number.of.samples),
  SampleOrigin = c(1:number.of.samples),
  medianUMI = c(1:number.of.samples),
  medianGenes = c(1:number.of.samples),
  medianMito =c(1:number.of.samples),
  cutoff = c(1:number.of.samples)
)

for(a in 1:length(list.of.seurats)){
  number.of.cells.per.sample[a,4]<-a
  number.of.cells.per.sample[a,1]<-names.data[a]
  number.of.cells.per.sample[a,2]<-length(Cells(list.of.seurats[[a]]))
  number.of.cells.per.sample[a,3]<-length(rownames(list.of.seurats[[a]]))
  number.of.cells.per.sample[a,5]<-stage.data[a]
  number.of.cells.per.sample[a,6]<-median(list.of.seurats[[a]]$nCount_RNA)
  number.of.cells.per.sample[a,7]<-median(list.of.seurats[[a]]$nFeature_RNA)
  number.of.cells.per.sample[a,8]<-median(list.of.seurats[[a]]$percent.mito)
  number.of.cells.per.sample[a,9]<-min.UMIs1[a]
  #print(number.of.cells.per.sample[a,])
}

print(number.of.cells.per.sample)

Cells.merged <- merge(
  x = list.of.seurats[[1]],
  y = list.of.seurats[2:length(list.of.seurats)])

VlnPlot(object = Cells.merged,
        #log = FALSE,
        pt.size = 0,
        group.by = "orig.ident",
        #y.max = 2000,
        #ncol = 3,
        #cols = c("blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue"),
        features = c("nFeature_RNA", "nCount_RNA", "percent.mito")
)

#Integration

# current.list <- c(HS_SCT_subset, list.of.seurats)
current.list <- list.of.seurats[c(1:9, 11)]

for (i in 1:length(x = current.list)) {
  current.list[[i]] <- NormalizeData(object = current.list[[i]], verbose = FALSE)
  current.list[[i]] <- FindVariableFeatures(object = current.list[[i]], nfeatures = 2000, verbose = FALSE)
}

current.anchors <- FindIntegrationAnchors(object.list = current.list, dims = 1:30)
# save(current.anchors, file="current.anchors.Robj")
current.integrated <- IntegrateData(anchorset = current.anchors, dims = 1:30)
# save(current.integrated, file="current.integrated.Robj")

DefaultAssay(object = current.integrated) <- "integrated"
current.integrated <- ScaleData(object = current.integrated, vars.to.regress = c("nCount_RNA", "percent.mito"), verbose = FALSE)
current.integrated <- RunPCA(object = current.integrated,
                             npcs = 30,
                             verbose = TRUE)
current.integrated <- FindNeighbors(object = current.integrated)
current.integrated <- FindClusters(object = current.integrated,
                                   resolution = 1.1,
                                   #algorithm = 1,
                                   verbose = TRUE
)
# current.integrated <- RunTSNE(object = current.integrated, dims = 1:20)
current.integrated <- RunUMAP(object = current.integrated, reduction = "pca",
                              dims = 1:30)

DimPlot(object = current.integrated, reduction = "umap", group.by = "stage", label = T)
DimPlot(object = current.integrated, reduction = "umap", group.by = "species", split.by = "experiment")

Idents(current.integrated) <- "species"
DimPlot(current.integrated)
current.integrated <- subset(current.integrated, idents = c("human"), invert = F)
DimPlot(current.integrated)

#recluster
DefaultAssay(current.integrated) <- "integrated"
current.integrated <- FindNeighbors(object = current.integrated)
current.integrated <- FindClusters(object = current.integrated,
                                   resolution = 1.1,
                                   #algorithm = 1,
                                   verbose = TRUE
)
current.integrated <- RunUMAP(object = current.integrated, reduction = "pca",
                              dims = 1:30)

DimPlot(object = current.integrated)

NCG_germ <- current.integrated

#CELL CYCLE REGRESSION


#These are Seurat's default list of cell cycle genes

s.genes <- (cc.genes.updated.2019$s.genes)
g2m.genes <- (cc.genes.updated.2019$g2m.genes)

#Cell cycle scoring:
DefaultAssay(NCG_germ) <- "RNA"
NCG_germ <- CellCycleScoring(NCG_germ, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)



fp1 <- FeaturePlot(NCG_germ, features = c("S.Score"), min.cutoff = 0) + labs(title = "S score")
fp2 <- FeaturePlot(NCG_germ, features = c("G2M.Score"), min.cutoff = 0) + labs(title = expression(bold(paste(G[2],"-M score",sep=""))))
wrap_plots(fp1, fp2)

DimPlot(NCG_germ, group.by = "Phase")

setwd("/Users/ewhelan/Documents/scRNAseq/Rscripts/savefiles")
save(NCG_germ, file = "NCG_germ_newCellNames.2023.07.27.Robj")
# load("NCG_germ_newCellNames.2023.07.27.Robj")

#Regress out the cell cycle genes
#This did not adequately remove the cell cycle effects, so it was achieved by splitting and re-integrating the samples, below, as described on the Satija lab website
# 
# DefaultAssay(NCG_germ) <- "integrated"
# NCG_germ2 <- NCG_germ
# 
# NCG_germ2 <- ScaleData(NCG_germ2, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(NCG_germ))
# 
# 
# 
# #Germ.culture2 is the variation with cell cycle regression
# DefaultAssay(NCG_germ2) <- "integrated"
# # NCG_germ2 <- FindVariableFeatures(object = NCG_germ2, 
# #                                       verbose = FALSE)
# n.pcs <- 30
# NCG_germ2 <- RunPCA(NCG_germ2, npcs = n.pcs)#, features = VariableFeatures(Germ.culture),npcs = 30, nfeatures.print = 10)
# NCG_germ2 <- FindNeighbors(object = NCG_germ2,
#                            k.param = n.pcs)
# NCG_germ2 <- FindClusters(object = NCG_germ2, 
#                           # resolution = 0.34,
#                           algorithm = 1,
#                           verbose = TRUE
# )
# # Germ.culture2<- RunTSNE(object = Germ.culture2, dims = 1:n.pcs, reduction = "pca")
# NCG_germ2 <- RunUMAP(object = NCG_germ2, 
#                      dims = 1:n.pcs,
#                      reduction = "pca")
# 
# 
# DefaultAssay(NCG_germ2) <- "RNA"
# DimPlot(NCG_germ2, #split.by = "experiment", 
#         group.by = "Phase")
# 



#Instead split and re-integrate samples
DefaultAssay(NCG_germ) <- "RNA"
Idents(NCG_germ) <- "Phase"
NCG_S <- subset(NCG_germ, idents = "S")
NCG_G1 <- subset(NCG_germ, idents = "G1")
NCG_G2M <- subset(NCG_germ, idents = "G2M")

n.dims <- 25

current.list <- c(NCG_S, NCG_G1, NCG_G2M)
for (i in 1:length(x = current.list)) {
  # current.list[[i]] <- NormalizeData(object = current.list[[i]], verbose = FALSE)
  current.list[[i]] <- FindVariableFeatures(object = current.list[[i]], nfeatures = 3000, verbose = FALSE)
}

current.anchors <- FindIntegrationAnchors(object.list = current.list, dims = 1:n.dims)
# save(current.anchors, file="current.anchors.Robj")
current.integrated <- IntegrateData(anchorset = current.anchors, dims = 1:n.dims)
# save(current.integrated, file="current.integrated.Robj")


DefaultAssay(object = current.integrated) <- "integrated"
current.integrated <- ScaleData(object = current.integrated, verbose = FALSE)
current.integrated <- RunPCA(object = current.integrated,
                             npcs = n.dims,
                             verbose = TRUE)
current.integrated <- FindNeighbors(object = current.integrated)
current.integrated <- FindClusters(object = current.integrated,
                                   resolution = 1.1,
                                   #algorithm = 1,
                                   verbose = TRUE
)
# current.integrated <- RunTSNE(object = current.integrated, dims = 1:20)
current.integrated <- RunUMAP(object = current.integrated, reduction = "pca",
                              dims = 1:n.dims)

DimPlot(object = current.integrated, reduction = "umap", group.by = "Phase", label = T)
DimPlot(object = current.integrated, reduction = "umap", group.by = "experiment")

NCG_germ3 <- current.integrated

setwd("/Users/ewhelan/Documents/scRNAseq/Rscripts/savefiles")
save(NCG_germ3, file="NCG_germ3_newCellNames_2023.07.27")
load("NCG_germ3_newCellNames_2023.07.27")

p1 <- DimPlot(NCG_germ3, group.by = "seurat_clusters", split.by = "stage")
p2 <- FeaturePlot(NCG_germ3, features = c("STRA8"), split.by = "stage") & DarkTheme() & scale_color_viridis(option = "C")
plot_grid(p1, p2, ncol = 1)
# "POU5F1", "TFAP2C", "DDX4", "STRA8", "PIWIL4", "SOHLH2", "SYCP3", "SYCP1", "MEIOB"

#examine gene expression:


Keren.gene.list <- c("POU5F1", "NANOG", "TCL1A", "PDPN", "PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1",
                     "UTF1", "KIT", "ASB9", "NANOS2", "PIWIL4", "MORC1", "MAGEB2", "TDRD1", "RHOXF1", "MAGEC2", "ZBTB16", "EGR4", "RET", "MAGEA4",
                     "SOHLH1", "DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8", "DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9", "ZCWPW1", "RAD51AP2", "SYCP3",
                     "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1", "MLH3", "PIWIL1", "MLH1", "SPATA8", "CCDC112", "OVOL1", "OVOL2",
                     "CCDC42", "ACR", "CDRT15", "TEKT1", "H1FOO", "TJP3", "FAM24A", "CATSPER3", "SPACA3", "SPACA1", "ACRV1", "SAXO1", "CAPZA3", "HEMGN", "CA2", "PRM3",
                     "TNP2", "TSSK6", "SMCP", "SPEM2", "LELP1", "MKI67", "UBE2C", "TOP2A")
#from Cell Reports paper "Defining the cellular origin of seminoma by transcriptional and epigenetic mapping to the normal human germline"

gene.list.PGC <- c("POU5F1", "NANOG", "TCL1A", "PDPN", "PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1",
                   "UTF1", "KIT", "ASB9", "NANOS2")
gene.list.spermatogonia <- c(                 
  "PIWIL4", "MORC1", "MAGEB2", "TDRD1", "RHOXF1", "MAGEC2", "ZBTB16", "EGR4", "RET")
gene.list.diff.spermatgonia <- c(
  
  "MAGEA4",
  "SOHLH1", "DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8")

gene.list.early.meiosis <- c(
  "DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9", "ZCWPW1", "RAD51AP2", "SYCP3",
  "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1", "MLH3")

gene.list.late.meiosis <- c(
  "PIWIL1", "MLH1", "SPATA8", "CCDC112", "OVOL1", "OVOL2",
  "CCDC42", "ACR", "CDRT15", "TEKT1", "H1FOO", "TJP3", "FAM24A", "CATSPER3", "SPACA3", "SPACA1", "ACRV1", "SAXO1")
gene.list.spermatids <- c("CAPZA3", "HEMGN", "CA2", "PRM3",
                          "TNP2", "TSSK6", "SMCP", "SPEM2", "LELP1", "MKI67", "UBE2C", "TOP2A")


gene.list.sertoli <- c("CLU", "AMHR2", "SOX9", "CTSL")
genes.macrophage <- c("APOE", "DAB2", "CD74", "ADGRE1")
genes.endothelial <- c("VWF", "TIE1", "TEK")
genes.myoid <- c("ACTA2", "MYH11", "MYL6", "PDGFRB")
genes.leydig <- c("CYP17A1", "CYP11A1", "STAR", "HSD3B1")
genes.T <- c("CD8A", "CD4", "CD3E", "CD28")
genes.telocytes <- c("DCN", "GSN", "TCF21")

germ.cells <- c("POU5F1", "NANOG", "ZBTB16", "SOHLH1", "PRDM9", "ACRV1", "SPACA1", "PRM1", "TNP2")
somatic.cells <- c("SOX9", "APOE", "VWF", "ACTA2", "STAR", "DCN", "TCF21", "ACTA2", "DAB2")

apoptotic.genes <- c("ADD1","AIFM3","ANKH","ANXA1","APP","ATF3","AVPR1A","BAX","BCAP31","BCL10","BCL2L1",
                     "BCL2L10","BCL2L11","BCL2L2","BGN","BID","BIK","BIRC3","BMF","BMP2","BNIP3L","BRCA1",
                     "BTG2","BTG3","CASP1","CASP2","CASP3","CASP4","CASP6","CASP7","CASP8","CASP9","CAV1",
                     "CCNA1","CCND1","CCND2","CD14","CD2","CD38","CD44","CD69","CDC25B","CDK2","CDKN1A","CDKN1B",
                     "CFLAR","CLU","CREBBP","CTH","CTNNB1","CYLD","DAP","DAP3","DCN","DDIT3","DFFA","DIABLO","DNAJA1","DNAJC3","DNM1L","DPYD","EBP",
                     "EGR3","EMP1","ENO2","ERBB2","ERBB3","EREG","ETF1","F2","F2R","FAS","FASLG","FDXR","FEZ1","GADD45A","GADD45B","GCH1","GNA15",#"GPX1",
                     "GPX3","GPX4","GSN","GSR","GSTM1","GUCY2D",#"H1-0",
                     "HGF","HMGB2","HMOX1","HSPB1","IER3","IFITM3","IFNB1","IFNGR1","IGF2R","IGFBP6",
                     "IL18","IL1A","IL1B","IL6","IRF1","ISG20","JUN","KRT18","LEF1","LGALS3","LMNA","PLPPR4","LUM","MADD","MCL1","MGMT","MMP2","NEDD9",
                     "NEFH","PAK1","PDCD4","PDGFRB","PEA15","PLAT","PLCB2","PMAIP1","PPP2R5B","PPP3R1","PPT1","PRF1","PSEN1","PSEN2","PTK2","RARA",
                     "RELA","RETSAT","RHOB","RHOT2","RNASEL","ROCK1","SAT1","SATB1","SC5D","SLC20A1","SMAD7","SOD1","SOD2","SPTAN1","SQSTM1","TAP1",
                     "TGFB2","TGFBR3","TIMP1","TIMP2","TIMP3","TNF","TNFRSF12A","TNFSF10","TOP2A","TSPO","TXNIP","VDAC2","WEE1","XIAP")


gene.list.PGC.E <- c("POU5F1", "NANOG", "TCL1A", "PDPN" )
gene.list.PGC.L <- c("PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1", "UTF1")
gene.list.prospermatogonia <- c("KIT", 
                                # "ASB9", 
                                "NANOS2")
gene.list.spermatogonia <- c(                 
  "PIWIL4", "MORC1", "MAGEB2", "TDRD1", "RHOXF1", "MAGEC2", "ZBTB16", "EGR4", "RET")

gene.list.diff.spermatgonia.E <- c(
  "MAGEA4",
  "SOHLH1")

gene.list.diff.spermatgonia.L   <- c("DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8")

gene.list.preleptotene <- c("DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9")

gene.list.early.meiosis <- c("ZCWPW1", "RAD51AP2", "SYCP3",
                             "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1", "MLH3")


NCG_germ3_temp <- NCG_germ3


#detect cells with the matching barcode found in fastq files
#I searched for cells with a Meiob sequence manually using BLAST and the fastq files and selected 8 cells at random and selected 4 of those to project.

cells.NCG <- Cells(NCG_germ3)
library(stringr)

cell.1 <- "CCACGTTGTTATGTCG"
# cell.2 <- "TAACACGCACCATTCC"
# cell.3 <- "TGCGGGTAGTCTACCA"
cell.4 <- "GGGAGATTCGCAGTTA"
# cell.5 <- "GTCGCGAGTCTGTGAT"
# cell.6 <- "CGTAGTATCAACACCA"
cell.7 <- "AAGCGTTGTGCAATAA"
cell.8 <- "CAACCAATCCCTTGGT" #cell8
table(str_detect(cells.NCG, cell.8))

cell.of.interest1 <- cells.NCG[str_detect(cells.NCG, cell.1)]
cell.of.interest4 <- cells.NCG[str_detect(cells.NCG, cell.4)]
cell.of.interest7 <- cells.NCG[str_detect(cells.NCG, cell.7)]
cell.of.interest8 <- cells.NCG[str_detect(cells.NCG, cell.8)]

cells.of.interest.list <- c(cell.of.interest1,
                            cell.of.interest4, cell.of.interest7, cell.of.interest8)


### FIGURE 3 F & G ####

DimPlot(NCG_germ4, cells.highlight = cells.of.interest.list)
FeaturePlot(NCG_germ4, features = "MEIOB")

### ### ### ### ### ###


#Next, remove doublet-heavy clusters and add module scores

DefaultAssay(NCG_germ3_temp) <- "RNA"

table(Idents(NCG_germ3_temp))

DefaultAssay(NCG_germ3_temp) <- "RNA"
suppressMessages(require(DoubletFinder))
library(DoubletFinder)
NCG_germ3_temp <- NormalizeData(NCG_germ3_temp)
NCG_germ3_temp = FindVariableFeatures(NCG_germ3_temp, verbose = F)
NCG_germ3_temp = ScaleData(NCG_germ3_temp, vars.to.regress = c("nFeature_RNA", "percent.mito"),
                           verbose = F)
NCG_germ3_temp = RunPCA(NCG_germ3_temp, verbose = F, npcs = 25)
NCG_germ3_temp = RunUMAP(NCG_germ3_temp, dims = 1:25, verbose = F)

nExp <- round(ncol(NCG_germ3_temp) * 0.04)  # expect 4% doublets
# nExp <- 1000 #increase this to account for closer to what clusters 8 2 and 4
NCG_germ3_temp <- doubletFinder_v3(NCG_germ3_temp, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

# name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(NCG_germ3_temp@meta.data)[grepl("DF.classification", colnames(NCG_germ3_temp@meta.data))]


cowplot::plot_grid(ncol = 2, DimPlot(NCG_germ3_temp, group.by = "orig.ident") + NoAxes(),
                   DimPlot(NCG_germ3_temp, group.by = DF.name) + NoAxes())

FeaturePlot(NCG_germ3_temp, features = gene.list, min.cutoff = 0, order = T)

NCG_germ3_temp2 <- NCG_germ3

NCG_germ3_temp2$doublets <- NCG_germ3_temp$DF.classifications_0.25_0.09_461

DimPlot(NCG_germ3_temp2, group.by = "doublets")

Idents(NCG_germ3_singlets) <- "seurat_clusters"
DimPlot(NCG_germ3_singlets, label = T)
VlnPlot(NCG_germ3_singlets, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"))


NCG_germ3_singlets <- subset(NCG_germ3_singlets,
                             subset = nCount_RNA < 70000)


#SUBSET to remove doublet-heavy clusters##

DimPlot(NCG_germ3, label = T)
NCG_germ4 <- subset(NCG_germ3, idents = c(14, 15, 3), invert = T) #7
DimPlot(NCG_germ4, label = T)

n.dims = 30
DefaultAssay(object = NCG_germ4) <- "RNA"
NCG_germ4 = NormalizeData(NCG_germ4)
NCG_germ4 = FindVariableFeatures(NCG_germ4, verbose = F)
NCG_germ4 <- ScaleData(object = NCG_germ4, verbose = T) #, vars.to.regress = c("nFeature_RNA", "percent.mito"))
NCG_germ4 <- RunPCA(object = NCG_germ4,
                    npcs = n.dims,
                    verbose = TRUE)
NCG_germ4 <- FindNeighbors(object = NCG_germ4)
NCG_germ4 <- FindClusters(object = NCG_germ4,
                          resolution = 1.1,
                          #algorithm = 1,
                          verbose = TRUE
)

# NCG_germ4 <- RunTSNE(object = NCG_germ4, reduction = "pca",
#                      dims = 1:n.dims)

NCG_germ4 <- RunUMAP(object = NCG_germ4, reduction = "pca",
                     dims = 1:n.dims)


DimPlot(NCG_germ4, reduction = "umap", label = T)

DefaultAssay(NCG_germ4) <- "RNA"
FeaturePlot(NCG_germ4, features = c("MAGEC2", "DDX4", "UTF1"),reduction = "umap", min.cutoff = 0, ncol = 3, order = T) & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(NCG_germ4, features = c("TFAP2C", "DDX4", "KIT"),reduction = "umap", min.cutoff = 0, ncol = 3, order = T) & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(NCG_germ4, features = c("MAGEC2", "DDX4", "UTF1"),reduction = "umap", min.cutoff = 0, ncol = 3, order = T) & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(NCG_germ4, features = c("TFAP2C", "DDX4", "NANOG"),reduction = "umap", min.cutoff = 0, ncol = 3, order = T) & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(NCG_germ4, features = c("TFAM", "DDX4", "MKI67"),reduction = "umap", min.cutoff = 0, ncol = 3, order = T) & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(NCG_germ4, features = c("TFAP2C", "DDX4", "MKI67"),reduction = "umap", min.cutoff = 0, ncol = 3, order = T) & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(NCG_germ4, features = c("MAGEC2", "DDX4", "GFRA1"),reduction = "umap", min.cutoff = 0, ncol = 3, order = T) & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(NCG_germ4, features = c("MKI67", "DDX4", "KIT"),reduction = "umap", min.cutoff = 0, ncol = 3, order = T) & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(NCG_germ4, features = c("REC8", "DDX4", "MAGEA6"),reduction = "umap", min.cutoff = 0, ncol = 3, order = T) & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(NCG_germ4, features = c("TFAP2C", "DDX4", "KIT"),reduction = "umap", min.cutoff = 0, ncol = 3, order = T) & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(NCG_germ4, features = c("SYCP3", "DDX4", "MAGEA6"),reduction = "umap", min.cutoff = 0, ncol = 3, order = T) & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(NCG_germ4, features = c("KIT", "DDX4", "PIWIL4"),reduction = "umap", min.cutoff = 0, ncol = 3, order = T) & DarkTheme() & scale_color_viridis(option = "C")

VlnPlot(NCG_germ4, features = c("TFAP2C", "UTF1"))

vector.of.cell.types <- c("PGCLCs", 
  "PGCs (early)", 
  "PGCs (late)", 
  "M prospermatogonia",
  "M/T1 prospermaogonia", 
  "T1 prospermatogonia", 
  "Spermatogonia",
  "Diff. spermatogonia (early)", 
  "Diff. spermatogonia (late)", 
  "Preleptotene")


Idents(NCG_germ4) <- "experiment"
NCG_germ4[["UMAP"]] <- NCG_germ4[["umap"]]
NCG_germ4subset <- subset(NCG_germ4, idents = "transplant")
Idents(NCG_germ4subset) <- "seurat_clusters"




###PSEUDOTIME STARTS HERE###



cds <- as.cell_data_set(NCG_germ4subset)
# cds <- as.cell_data_set(NCG_germ4)


cds <- cluster_cells(cds, reduction_method = c("UMAP"), #cluster_method = "louvain",
                     resolution = 0.2e-4,
                     k = 8,
                     # k = 6,
                     partition_qval = 0.05
)

plot_cells(cds, group_label_size = 5,
           show_trajectory_graph = FALSE)


# for(a in 25:150){
cds <- learn_graph(cds,
                   close_loop = F,
                   learn_graph_control=list(ncenter=50, minimal_branch_len=1), #100, 3 is pretty good, 30, 2 is also.
                   use_partition = F)
print(plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)) #+ ggtitle(a))
# }
# Figure 2C #######################################
cds <- order_cells(cds)#, root_cells = max.avp)
p1 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,
           label_branch_points = FALSE)+ NoLegend()

p2 <- DimPlot(NCG_germ4, group.by = "cell.type") + NoLegend()

plot_grid(p2, p1)



pseudotime.values <- pseudotime(cds, reduction_method = "UMAP")


NCG_germ4[["pseudotime"]] <- pseudotime.values




#add module scores


gene.list.PGC.E <- c("POU5F1", "NANOS3", "MIF", "ACTG2")
gene.list.PGC.L <- c("SOX15", "L1TD1", "PDPN", "TUBA1C")
gene.list.M <- c("KIT", "CCND2", "EEF1G", "CCPG1" )
gene.list.M2T1 <- c("MAP4K4", "ASB9", "PRXL2A", "TLE5")
gene.list.T1 <- c("DCAF4L1", "ID1", "ID3", "MORC1")
gene.list.spermatogonia <- c("UTF1", "PAGE2", "RHOXF1", "MAGEB2")
gene.list.diff.spermatgonia.E <- c("SOHLH1", "KCNQ2", "MAGEC2", "MAGEA4")
gene.list.diff.spermatgonia.L   <- c("DMRT1", "SOHLH2", "ESX1", "CTCFL")
gene.list.preleptotene <- c("RAD51AP2", "SCML1", "TEX19", "PRDM9")
gene.list.leptotene <- c("TEX101", "SYCP1", "DMRTC2", "TDRG1")
gene.list.zygotene <- c("CTAGE1", "EQTN", "SPACA1", "GOLGA6L2")


NCG_germ4 <- AddModuleScore(NCG_germ4,
                                 features = list(apoptotic.genes),
                                 name="Apoptotic")
NCG_germ4 <- AddModuleScore(NCG_germ4,
                                 features = list(gene.list.PGC.E),
                                 name="PGCE")
NCG_germ4 <- AddModuleScore(NCG_germ4,
                                 features = list(gene.list.PGC.L),
                                 name="PGCL")
NCG_germ4 <- AddModuleScore(NCG_germ4,
                            features = list(gene.list.M),
                            name="M")
NCG_germ4 <- AddModuleScore(NCG_germ4,
                            features = list(gene.list.M2T1),
                            name="M2T")
NCG_germ4 <- AddModuleScore(NCG_germ4,
                            features = list(gene.list.T1),
                            name="T")
NCG_germ4 <- AddModuleScore(NCG_germ4,
                                 features = list(gene.list.spermatogonia),
                                 name="spg")
NCG_germ4 <- AddModuleScore(NCG_germ4,
                                 features = list(gene.list.diff.spermatgonia.E),
                                 name="diffspgE")
NCG_germ4 <- AddModuleScore(NCG_germ4,
                                 features = list(gene.list.diff.spermatgonia.L),
                                 name="diffspgL")
NCG_germ4 <- AddModuleScore(NCG_germ4,
                                 features = list(gene.list.preleptotene),
                                 name="preleptotene")
NCG_germ4 <- AddModuleScore(NCG_germ4,
                                 features = list(gene.list.leptotene),
                                 name="leptotene")
NCG_germ4 <- AddModuleScore(NCG_germ4,
                            features = list(gene.list.zygotene),
                            name="zygotene")


# scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

# Plot scores

plot_grid(
  
  FeaturePlot(NCG_germ4,
              features = "Apoptotic1",min.cutoff = 0,max.cutoff = 2, order = T, label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ4,
              features = "PGCE1",min.cutoff = 0, max.cutoff = 2,order = T, label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ4,
              features = "PGCL1",min.cutoff = 0, max.cutoff = 2,order = T, label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ4,
              features = "M1", min.cutoff = 0, max.cutoff = 2,order = T,label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  FeaturePlot(NCG_germ4,
              features = "M2T1", min.cutoff = 0, max.cutoff = 2,order = T,label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  FeaturePlot(NCG_germ4,
              features = "T1", min.cutoff = 0, max.cutoff = 2,order = T,label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ4,
              features = "spg1",min.cutoff = 0, max.cutoff = 2,order = T, label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ4,
              features = "diffspgE1",min.cutoff = 0, max.cutoff = 2,order = T, label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ4,
              features = "diffspgL1",min.cutoff = 0, max.cutoff = 2,order = T, label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ4,
              features = "preleptotene1", min.cutoff = 0,max.cutoff = 2, order = T,label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ4,
              features = "leptotene1", min.cutoff = 0,max.cutoff = 2,order = T,label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ4,
              features = "zygotene1", min.cutoff = 0, max.cutoff = 2,order = T,label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme()
  
)

# Assign idents 

Idents(NCG_germ4) <- "seurat_clusters"
DimPlot(NCG_germ4, label = T)
new.cluster.ids <- character(length(levels(Idents(NCG_germ4))))

PGCLC.ids <- c(10, 15, 13)+1
PGC.E.ids <- c(8)+1
PGC.L.ids <- c(4, 11, 5)+1
M.ids <- c(3, 1)+1
M2T1.ids <- c(12, 6)+1
T1.ids <- c(9)+1
spermatogonia.ids <- c(0, 2)+1
earlydiff.ids <- c(7)+1
latediff.ids <- c(14, 16)+1
preleplep.ids <- c(17)+1

new.cluster.ids[PGCLC.ids] <- "PGCLCs"
new.cluster.ids[PGC.E.ids] <- "PGCs (early)"
new.cluster.ids[PGC.L.ids] <- "PGCs (late)"
new.cluster.ids[M.ids] <- "M prospermatogonia"
new.cluster.ids[M2T1.ids] <- "M/T1 prospermaogonia"
new.cluster.ids[T1.ids] <- "T1 prospermatogonia"
new.cluster.ids[spermatogonia.ids] <- "Spermatogonia"
new.cluster.ids[earlydiff.ids] <- "Diff. spermatogonia (early)"
new.cluster.ids[latediff.ids] <- "Diff. spermatogonia (late)"
new.cluster.ids[preleplep.ids] <- "Preleptotene spermatocytes"



new.cluster.ids

names(x = new.cluster.ids) <- levels(x = NCG_germ4)
NCG_germ4 <- RenameIdents(object = NCG_germ4, new.cluster.ids)
NCG_germ4$cell.type <- Idents(NCG_germ4)
DimPlot(NCG_germ4)
NCG_germ4$cell.type <- factor(x = NCG_germ4$cell.type, levels = c("PGCLCs", 
                                                                  "PGCs (early)", 
                                                                  "PGCs (late)",
                                                                  "M prospermatogonia",
                                                                  "M/T1 prospermaogonia",
                                                                  "T1 prospermatogonia",
                                                                  "Spermatogonia",
                                                                  "Diff. spermatogonia (early)",
                                                                  "Diff. spermatogonia (late)",
                                                                  "Preleptotene spermatocytes"
))
Idents(NCG_germ4) <- NCG_germ4$cell.type
#Plot out designations with a feature of interest.
DimPlot(NCG_germ4, label = F, reduction = "umap", group.by = "cell.type")

setwd("~/Documents/scRNAseq/Rscripts/savefiles")
# save(NCG_germ4, file="NCG_germ4_newCellNames_2023.07.27.Robj")
load("NCG_germ4_newCellNames_2023.07.27.Robj")


NCG_germ4_newcells <- NCG_germ4
new.cell.names <- (Cells(NCG_germ4_newcells))
old.cell.names <- (Cells(NCG_germ4))


#HEATMAP

Keren.gene.list <- c("POU5F1", "NANOG", "TCL1A", "PDPN", "PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1","DDX4",
                     "UTF1", "KIT", "ASB9", "NANOS2", "PIWIL4", "MORC1", "MAGEB2","MAGEA3", "TDRD1", "RHOXF1", "MAGEC2", "GFRA1", "ZBTB16", "EGR4", "RET", "MAGEA4",
                     "SOHLH1", "DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8", "DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9", "ZCWPW1", "RAD51AP2", "SYCP3",
                     "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1", "MLH3", "PIWIL1", "MLH1", "SPATA8", "CCDC112", "OVOL1", "OVOL2",
                     "CCDC42", "ACR", "CDRT15", "TEKT1", "H1FOO", "TJP3", "FAM24A", "CATSPER3", "SPACA3", "SPACA1", "ACRV1", "SAXO1", "CAPZA3", "HEMGN", "CA2", "PRM3",
                     "TNP2", "TSSK6", "SMCP", "SPEM2", "LELP1", "MKI67", "UBE2C", "TOP2A")

Keren.scatter.PGCE <- c("HMGA2", "ABHD12B", "SPRR2F", "DDT", "PRDM1", "KHDC3L", "TFAP2C", "NANOG", "YBX1", "MIF")
Keren.scatter.PGCL1 <- c("TEX14", "GAGE12H", "GAGE2A", "XIST", "DCAF4L1", "MAEL", "DAZL", "SOX17", "UTF1")
Keren.scatter.PGCL2 <- c("CCND1", "ID3", "MAGEB2", "FHL2", "ID1", "MAGED2", "TFAP2C", "SOX17", "NANOG", "UTF1", "POU5F1")
Keren.scatter.M1 <- c("CDKN2C", "GCNA", "E2F8", "HIST1H3B", "SYCP2", "FOXM1", "CDC25A", "XIST", "MKI67")
Keren.scatter.M2 <- c("PRDM14", "ACTG2", "KHDC3L", "UBE2C", "DPPA5", "GPHA2", "MKI67", "RASD1", "CENPA", "NANOG", "SFRP2", "POU5F1", "APOC1", "PRDM1")
Keren.scatter.trans1 <- c("PAGE2", "ASB9", "RHOXF1", "MORC1", "TEX15", "KRBOX1", "PAGE2B", "TDRD1", "ID1", "STK31")
Keren.scatter.trans2 <- c("DGKK", "NANOS3", "ASB9", "MDK", "MAP4K4", "STMN1", "IGFBP2")
Keren.scatter.T1 <- c("FOXD1", "DPPA5", "EGR4", "SIX1", "NKX6-2", "PIWIL4", "VCX3B", "MAGEC2", "DPPA2", "MORC1", "MAGEB2")
Keren.scatter.T2 <- c("GAGE12H", "NANOS2", "MEG3", "LY86-AS1", "XIST", "HELLPAR", "DCAF4L1")
Keren.scatter.spg1 <- c("RET", "HMGB4", "POU3F1", "TSPAN33", "FOXD1", "C19orf84", "SCT", "KCNQ2", "MAGEA4", "VCX", "TCF3", "UTF1", "PIWIL4")
Keren.scatter.spg2 <- c("RET", "HMGB4", "FIGLA", "PHOX2A", "HES5", "MAPK3", "EGR4", "PIWIL4", "TSPAN33", "C19orf84", "RGS14", "POU3F1", "FOXD1", "SCT", "MORC1", "UTF1")
Keren.scatter.diffspgE1 <- c("HIST1H3B", "HMGB4", "EGR1", "EIF3K", "ASB9", "HIST1H1A", "L1TD1", "MAGEA4")
Keren.scatter.diffspgE2 <- c("RESF1", "DUSP6", "ID3", "EGR4", "SIX1", "DPPA2", "C17orf49", "FGFR3", "ID2", "PCSK1N", "ASB9", "MAGEA4")
Keren.scatter.diffspgL1 <- c("STRA8", "SYCE2", "TEX19", "MKI67", "KIT", "HORMAD1", "MEIOC", "SCML1", "CENPU", "TOP2A", "SYCP3", "SMC1B", "HSPA5", "HIST1H1A")
Keren.scatter.diffspgL2 <- c("CCND1", "SOHLH1", "HIST1H2AB", "UBE2C", "KIT", "RHOXF1", "ESX1", "ITGA6", "TDRD1", "CENPF", "ELACVL2", "UCHL1")
Keren.scatter.prelep1 <- c("HOTAIR", "PRDM9", "FBXO47", "SPO11", "HORMAD2", "RNF212B", "CCDC172", "TEX12", "MEIOB", "SYCP1", "TEX101", "SYCP2", "RAD51AP2", "ZCWPW1")

Keren.gene.list2 <- unique(c(Keren.scatter.PGCE,Keren.scatter.PGCL1,Keren.scatter.PGCL2,Keren.scatter.M1, Keren.scatter.M2, Keren.scatter.trans1,
                           Keren.scatter.trans2, Keren.scatter.T1, Keren.scatter.T2, Keren.scatter.spg1, Keren.scatter.spg2, Keren.scatter.diffspgE1, 
                           Keren.scatter.diffspgE2, Keren.scatter.diffspgL1, Keren.scatter.diffspgL2, Keren.scatter.prelep1))


PGCLC.markers <- c("UBE2C", "SOX17", "NANOG", "OVOL2",
"PDPN", "L1TD1", "TFAP2C", "POU5F1",
"PEG10", "CA2")
PGCE.markers <- c("TCL1A", "NANOS3", "KHDC3L", "GPHA2", 
                  "APOC1", "KHDC3L", "ACTG2")
PGCL.markers <- c("TJP3", "TOP2A", "MKI67", 
                  "DPPA5", "CENPA"  )
M.markers <- c("GAGE2A", "GAGE12H", "SDC4",
               "SPATA16", "MAPK11", "ETV5")
intermediate.markers<- c("ASB9", "DGKK", "FOXA2", 
                         "TSPAN8", "CITED2", "FATE1")
T1.markers <-c("TDRD1", "MORC1", "SOHLH1",
               "DDX4", "EGR4", "ZCWPW1", 
               "RHOXF1", "NANOS2", "ID1", 
               "TEX15", "RHOXF1" 
)
Spg.markers <-c("MAGEB2", "MAGEC2",
                "CDRT15", "SPACA3", "SPO11", 
                "PIWIL4", "TEKT1")
DiffE.markers <- c("MAGEA4", "SCML2", "ESX1", "CCDC42",
                   "DAZL", "DMRTB1", "SPATA8", 
                   "SOHLH2", "PIWIL1","PAGE2" )
DiffL.markers <- c("STRA8", "MYBL1", "DMC1",
                   "SYCP3", "GAL3ST1", "HORMAD1",
                   "TEX19", "KRBOX1" )
Prelep.markers <- c("RAD51AP2", "MEIOB", "PRDM9",
                    "CAPZA3", "SYCP1")

FeaturePlot(NCG_germ4, features = Spg.markers, order = T) & DarkTheme() & scale_color_viridis(option = "C")


DimPlot(NCG_germ4)
Idents(NCG_germ4) <- (NCG_germ4$experiment)
# NCG_germ4_subset <- subset(NCG_germ4, idents = c("hPGCLC"), invert = T)
NCG_germ4_subset <- NCG_germ4
Idents(NCG_germ4_subset) <- NCG_germ4_subset$cell.type
DimPlot(NCG_germ4_subset)
# all.markers.subset <- filter(all.markers_NCG, avg_log2FC > 1)
PGCLCs.markers.all <- FindMarkers(NCG_germ4_subset, ident.1 = "PGCLCs")
PGCs.E.markers.all <- FindMarkers(NCG_germ4_subset, ident.1 = "PGCs (early)")
PGCs.L.markers.all <- FindMarkers(NCG_germ4_subset, ident.1 = "PGCs (late)")
M.ProSpermatogonia.markers.all <- FindMarkers(NCG_germ4_subset, ident.1 = "M prospermatogonia")
M2T1.markers.all <- FindMarkers(NCG_germ4_subset, ident.1 = "M/T1 prospermatogonia")
T.ProSpermatogonia.markers.all <- FindMarkers(NCG_germ4_subset, ident.1 = "T1 prospermatogonia")
Spg.markers.all <- FindMarkers(NCG_germ4_subset, ident.1 = "Spermatogonia")
Diff.Spg.E.markers.all <- FindMarkers(NCG_germ4_subset, ident.1 = "Diff. spermatogonia (early)")
Diff.Spg.L.markers.all <- FindMarkers(NCG_germ4_subset, ident.1 = "Diff. spermatogonia (late)")
Pre.Leptotene.markers.all <- FindMarkers(NCG_germ4_subset, ident.1 = "Preleptotene spermatocytes")

PGCLCs.markers.all$gene <- rownames(PGCLCs.markers.all) 
PGCs.E.markers.all$gene <- rownames(PGCs.E.markers.all) 
PGCs.L.markers.all$gene <- rownames(PGCs.L.markers.all) 
M.ProSpermatogonia.markers.all$gene <- rownames(M.ProSpermatogonia.markers.all) 
M2T1.markers.all$gene <- rownames(M2T1.markers.all)
T.ProSpermatogonia.markers.all$gene <- rownames(T.ProSpermatogonia.markers.all) 
Spg.markers.all$gene <- rownames(Spg.markers.all) 
Diff.Spg.E.markers.all$gene <- rownames(Diff.Spg.E.markers.all) 
Diff.Spg.L.markers.all$gene <- rownames(Diff.Spg.L.markers.all) 
Pre.Leptotene.markers.all$gene <- rownames(Pre.Leptotene.markers.all) 

# setwd("/Users/ewhelan/Desktop/xrTestisPaper/NCG/NCG_marker_genes")
# setwd("/Users/ewhelan/Desktop")
write.csv(PGCs.E.markers.all, "PGCs.E.markers.all.csv")
write.csv(PGCs.L.markers.all, "PGCs.L.markers.all.csv")
write.csv(M.ProSpermatogonia.markers.all, "M.ProSpermatogonia.markers.all.csv")
write.csv(M2T1.markers.all, "M2T1.markers.all.csv")
write.csv(T.ProSpermatogonia.markers.all, "T.ProSpermatogonia.markers.all.csv")
write.csv(Spg.markers.all, "Spg.markers.all.csv")
write.csv(Diff.Spg.E.markers.all, "Diff.Spg.E.markers.all.csv")
write.csv(Diff.Spg.L.markers.all, "Diff.Spg.L.markers.all.csv")
write.csv(Pre.Leptotene.markers.all, "Pre.Leptotene.markers.all.csv")

#instead, pick the top passing filter genes.

log2FC_cutoff <- 0.585

PGCLCs.markers.ordered <- PGCLCs.markers.all %>% filter(p_val_adj < 0.05 & avg_log2FC > log2FC_cutoff)
PGCs.E.markers.ordered <- PGCs.E.markers.all %>% filter(p_val_adj < 0.05 & avg_log2FC > log2FC_cutoff )
PGCs.L.markers.ordered <- PGCs.L.markers.all %>% filter(p_val_adj < 0.05 & avg_log2FC > log2FC_cutoff )
M.ProSpermatogonia.markers.ordered <- M.ProSpermatogonia.markers.all %>% filter(p_val_adj < 0.05 & avg_log2FC > log2FC_cutoff )
M2T1.markers.ordered <- M2T1.markers.all %>% filter(p_val_adj < 0.05 & avg_log2FC > log2FC_cutoff )
T.ProSpermatogonia.markers.ordered <- T.ProSpermatogonia.markers.all %>% filter(p_val_adj < 0.05 & avg_log2FC > log2FC_cutoff )
Spg.markers.ordered <- Spg.markers.all %>% filter(p_val_adj < 0.05 & avg_log2FC > log2FC_cutoff)
Diff.Spg.E.markers.ordered <- Diff.Spg.E.markers.all %>% filter(p_val_adj < 0.05 & avg_log2FC > log2FC_cutoff)
Diff.Spg.L.markers.ordered <- Diff.Spg.L.markers.all %>% filter(p_val_adj < 0.05 & avg_log2FC > log2FC_cutoff)
Pre.Leptotene.markers.ordered <- Pre.Leptotene.markers.all %>% filter(p_val_adj < 0.05 & avg_log2FC > log2FC_cutoff)

#pick the top X genes per cluster
genes.per.cluster <- 3000
PGCLCs.markers.ordered <- PGCLCs.markers.ordered %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
PGCs.E.markers.ordered <- PGCs.E.markers.ordered %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
PGCs.L.markers.ordered <- PGCs.L.markers.ordered %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
M.ProSpermatogonia.markers.ordered <- M.ProSpermatogonia.markers.ordered %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
M2T1.markers.ordered <- M2T1.markers.ordered %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
T.ProSpermatogonia.markers.ordered <- T.ProSpermatogonia.markers.ordered %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
Spg.markers.ordered <- Spg.markers.ordered %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
Diff.Spg.E.markers.ordered <- Diff.Spg.E.markers.ordered %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
Diff.Spg.L.markers.ordered <- Diff.Spg.L.markers.ordered %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
Pre.Leptotene.markers.ordered <- Pre.Leptotene.markers.ordered %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)

PGCLCs.markers.ordered$cell.type <- "PGCLC" 
PGCs.E.markers.ordered$cell.type <- "PGCs.E" 
PGCs.L.markers.ordered$cell.type <- "PGCs.L" 
M.ProSpermatogonia.markers.ordered$cell.type <- "M" 
M2T1.markers.ordered$cell.type <- "M2T1" 
T.ProSpermatogonia.markers.ordered$cell.type <- "T"
Spg.markers.ordered$cell.type <- "Spg"
Diff.Spg.E.markers.ordered$cell.type <- "Diff.Spg.E"
Diff.Spg.L.markers.ordered$cell.type <- "Diff.Spg.L"
Pre.Leptotene.markers.ordered$cell.type <- "Prelep.lep"

dim(PGCLCs.markers.ordered)
dim(PGCs.E.markers.ordered)
dim(PGCs.L.markers.ordered)
dim(M.ProSpermatogonia.markers.ordered)
dim(M2T1.markers.ordered)
dim(T.ProSpermatogonia.markers.ordered)
dim(Spg.markers.ordered)
dim(Diff.Spg.E.markers.ordered)
dim(Diff.Spg.L.markers.ordered)
dim(Pre.Leptotene.markers.ordered)


topNmarkers <- rbind(PGCLCs.markers.ordered, PGCs.E.markers.ordered, PGCs.L.markers.ordered, M.ProSpermatogonia.markers.ordered, M2T1.markers.ordered,
                     T.ProSpermatogonia.markers.ordered, Spg.markers.ordered, Diff.Spg.E.markers.ordered, Diff.Spg.L.markers.ordered, Pre.Leptotene.markers.ordered)

table(topNmarkers$cell.type)

intersect(rownames(topNmarkers), Keren.gene.list)
intersect(rownames(topNmarkers), Keren.gene.list2)
intersect(rownames(topNmarkers), "KCND1")
intersect("KCND1", Keren.gene.list2)

Keren.scatter.PGCE <- c("HMGA2", "ABHD12B", "SPRR2F", "DDT", "PRDM1", "KHDC3L", "TFAP2C", "NANOG", "YBX1", "MIF")
Keren.scatter.PGCL1 <- c("TEX14", "GAGE12H", "GAGE2A", "XIST", "DCAF4L1", "MAEL", "DAZL", "SOX17", "UTF1")
Keren.scatter.PGCL2 <- c("CCND1", "ID3", "MAGEB2", "FHL2", "ID1", "MAGED2", "TFAP2C", "SOX17", "NANOG", "UTF1", "POU5F1")
Keren.scatter.M1 <- c("CDKN2C", "GCNA", "E2F8", "HIST1H3B", "SYCP2", "FOXM1", "CDC25A", "XIST", "MKI67")
Keren.scatter.M2 <- c("PRDM14", "ACTG2", "KHDC3L", "UBE2C", "DPPA5", "GPHA2", "MKI67", "RASD1", "CENPA", "NANOG", "SFRP2", "POU5F1", "APOC1", "PRDM1")
Keren.scatter.trans1 <- c("PAGE2", "ASB9", "RHOXF1", "MORC1", "TEX15", "KRBOX1", "PAGE2B", "TDRD1", "ID1", "STK31")
Keren.scatter.trans2 <- c("DGKK", "NANOS3", "ASB9", "MDK", "MAP4K4", "STMN1", "IGFBP2")
Keren.scatter.T1 <- c("FOXD1", "DPPA5", "EGR4", "SIX1", "NKX6-2", "PIWIL4", "VCX3B", "MAGEC2", "DPPA2", "MORC1", "MAGEB2")
Keren.scatter.T2 <- c("GAGE12H", "NANOS2", "MEG3", "LY86-AS1", "XIST", "HELLPAR", "DCAF4L1")
Keren.scatter.spg1 <- c("RET", "HMGB4", "POU3F1", "TSPAN33", "FOXD1", "C19orf84", "SCT", "KCNQ2", "MAGEA4", "VCX", "TCF3", "UTF1", "PIWIL4")
Keren.scatter.spg2 <- c("RET", "HMGB4", "FIGLA", "PHOX2A", "HES5", "MAPK3", "EGR4", "PIWIL4", "TSPAN33", "C19orf84", "RGS14", "POU3F1", "FOXD1", "SCT", "MORC1", "UTF1")
Keren.scatter.diffspgE1 <- c("HIST1H3B", "HMGB4", "EGR1", "EIF3K", "ASB9", "HIST1H1A", "L1TD1", "MAGEA4")
Keren.scatter.diffspgE2 <- c("RESF1", "DUSP6", "ID3", "EGR4", "SIX1", "DPPA2", "C17orf49", "FGFR3", "ID2", "PCSK1N", "ASB9", "MAGEA4")
Keren.scatter.diffspgL1 <- c("STRA8", "SYCE2", "TEX19", "MKI67", "KIT", "HORMAD1", "MEIOC", "SCML1", "CENPU", "TOP2A", "SYCP3", "SMC1B", "HSPA5", "HIST1H1A")
Keren.scatter.diffspgL2 <- c("CCND1", "SOHLH1", "HIST1H2AB", "UBE2C", "KIT", "RHOXF1", "ESX1", "ITGA6", "TDRD1", "CENPF", "ELACVL2", "UCHL1")
Keren.scatter.prelep1 <- c("HOTAIR", "PRDM9", "FBXO47", "SPO11", "HORMAD2", "RNF212B", "CCDC172", "TEX12", "MEIOB", "SYCP1", "TEX101", "SYCP2", "RAD51AP2", "ZCWPW1")
# intersect(Keren.gene.list, "ACR")

# write.csv(topNmarkers, file=paste0("top", genes.per.cluster, "markers.csv"))
write.csv(topNmarkers, file="cutoff.markers.log2fc0.585.cap.3000.csv")

# Heatmap of select genes by cluster ----

Idents(NCG_germ4_subset) <- NCG_germ4_subset$cell.type
DimPlot(NCG_germ4_subset)

#fix the typo in M/T1 propermaogonia
levels(NCG_germ4_subset$cell.type)
levels(NCG_germ4_subset$cell.type)[levels(NCG_germ4_subset$cell.type) == "M/T1 prospermaogonia"] <- "M/T1 prospermatogonia"
levels(NCG_germ4_subset$cell.type)
Idents(NCG_germ4_subset) <- NCG_germ4_subset$cell.type

data.PGCLC <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = "PGCLCs")])
data.pgcE <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = "PGCs (early)")])
data.pgcL <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = "PGCs (late)")])
data.m <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = "M prospermatogonia")])
data.mt1 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = "M/T1 prospermatogonia")])
data.t1 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = "T1 prospermatogonia")])
data.undiff <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = "Spermatogonia")])
data.diffE <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = "Diff. spermatogonia (early)")])
data.diffL <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = "Diff. spermatogonia (late)")])
data.Spermatocytes <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = "Preleptotene spermatocytes")])

av.counts.PGCLC <- apply(data.PGCLC, 1, mean)
av.counts.pgcE <- apply(data.pgcE, 1, mean)
av.counts.pgcL <- apply(data.pgcL, 1, mean)
av.counts.m <- apply(data.m, 1, mean)
av.counts.mt1 <- apply(data.mt1, 1, mean)
av.counts.t1 <- apply(data.t1, 1, mean)
av.counts.undiff <- apply(data.undiff, 1, mean)
av.counts.diffE <- apply(data.diffE, 1, mean)
av.counts.diffL <- apply(data.diffL, 1, mean)
av.counts.Spermatocytes <- apply(data.Spermatocytes, 1, mean)

av.counts.HS_NCG_only_NCG <- cbind(av.counts.PGCLC, av.counts.pgcE, av.counts.pgcL, av.counts.m, av.counts.mt1,av.counts.t1,
                                   av.counts.undiff, av.counts.diffE, av.counts.diffL, av.counts.Spermatocytes)


av.counts.HS_NCG_only_NCG.df <- as.data.frame(av.counts.HS_NCG_only_NCG)

gene.name.list <- rownames(av.counts.HS_NCG_only_NCG.df)

av.counts.HS_NCG_only_NCG.df$gene.name <- gene.name.list

# genes.of.interest <- Keren.gene.list[c(1:48)] #this is just Keren's list of genes
# genes.of.interest <- (all.markers.subset$gene) #this is a log-fold cutoff
# genes.of.interest <- rownames(topNmarkers) #this is just taking the top N markers
# genes.of.interest <- topNmarkers$gene #this is just taking the top N markers
genes.of.interest <- "DNMT3A"

av.counts.genes.of.interest <- av.counts.HS_NCG_only_NCG.df %>%
  dplyr::filter(gene.name %in% genes.of.interest)


av.counts.genes.of.interest.ordered <- av.counts.genes.of.interest[match(genes.of.interest, av.counts.genes.of.interest$gene.name),]

integer.col <- ncol(av.counts.genes.of.interest.ordered)

av.counts.genes.of.interest.ordered <- av.counts.genes.of.interest.ordered[,-integer.col]

av.counts.genes.of.interest.matrix1 <- as.matrix(av.counts.genes.of.interest.ordered)
colnames(av.counts.genes.of.interest.matrix1) <- c("PGCLCs", "PGCs (early)", "PGCs (late)", "M prospermatogonia", 
                                                   "M/T1 prospermaogonia", "T1 prospermatogonia", "Spermatogonia", 
                                                   "Diff. spermatogonia (early)", "Diff. spermatogonia (late)", "Preleptotene spermatocytes")

av.counts.genes.of.interest.matrix.no.na <- na.omit(av.counts.genes.of.interest.matrix1)

head(av.counts.genes.of.interest.matrix.no.na)




dev.off()
l.mat <- rbind(c(0,4), c(0, 3), c(2,1), c(0,0))
l.mat

# heatmap.2(av.counts.genes.of.interest.matrix.filtered, Colv = FALSE,
#           Rowv = FALSE,
#           trace = "none",
#           # col =  colorRampPalette(colors=c("blue","white","red"))(100),
#           col = inferno(100),
#           density.info="none",
#           srtCol = -60, #angle of labels
#           adjCol = c(0,0), #offset of labels, should be centered now
#           cexRow=1.0, cexCol=1.3,
#           lmat= l.mat,
#           lwid = c(0.5, 5),
#           lhei=c(1, 1, 5, 1),
#           scale = "row"
#           
# )



av.counts.genes.of.interest.matrix.no.na.t <- as.data.frame(t(av.counts.genes.of.interest.matrix.no.na))
av.counts.genes.of.interest.matrix.no.na.t$cell.type <- rownames(av.counts.genes.of.interest.matrix.no.na.t)

av.counts.genes.of.interest.matrix.no.na.t$cell.type <- factor(av.counts.genes.of.interest.matrix.no.na.t$cell.type, levels = c(
  "PGCLCs", "PGCs (early)", "PGCs (late)", "M prospermatogonia", 
  "M/T1 prospermaogonia", "T1 prospermatogonia", "Spermatogonia", 
  "Diff. spermatogonia (early)", "Diff. spermatogonia (late)", 
  "Preleptotene spermatocytes"
))

av.counts <- av.counts.genes.of.interest.matrix.no.na.t

# Reorder the levels of cell.type according to the specified order
av.counts$cell.type <- factor(av.counts$cell.type, levels = c(
  "PGCLCs", "PGCs (early)", "PGCs (late)", "M prospermatogonia", 
  "M/T1 prospermaogonia", "T1 prospermatogonia", "Spermatogonia", 
  "Diff. spermatogonia (early)", "Diff. spermatogonia (late)", 
  "Preleptotene spermatocytes"
))

# Create the line plot
line_plot <- ggplot(av.counts, aes(x = cell.type, y = DNMT3A, group = 1)) +  # Set group = 1
  geom_line() +
  geom_point() +  # Add points to show individual values
  labs(title = "Expression of TET1 by Cell Type",
       x = "Cell Type",
       y = "TET1 Expression") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Print the plot
print(line_plot)

FeaturePlot(NCG_germ4, features = "APOBEC1") & scale_color_viridis(option = "C") & DarkTheme()


#repeat for cluster heatmap
Idents(NCG_germ4) <- "cell.type"
DimPlot(NCG_germ4)
DimPlot(NCG_germ4, group.by = "seurat_clusters", label = T)
NCG_germ4_subset <- subset(NCG_germ4, idents = c("PGCLCs"), invert = T)
NCG_germ4_subset <- NCG_germ4
Idents(NCG_germ4_subset) <- "seurat_clusters"
DimPlot(NCG_germ4_subset, label = T)

order.of.clusters <- c(8, 4, 11, 5, 3, 1, 12, 6, 9, 0, 2, 7, 14, 16, 17)


# all.markers.subset <- filter(all.markers_NCG, avg_log2FC > 1)
cluster0.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 0)
cluster1.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 1)
cluster2.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 2)
cluster3.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 3)
cluster4.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 4)
cluster5.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 5)
cluster6.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 6)
cluster7.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 7)
cluster8.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 8)
cluster9.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 9)
cluster10.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 10) 
cluster11.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 11)
cluster12.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 12)
cluster13.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 13)
cluster14.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 14)
cluster15.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 15)
cluster16.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 16)
cluster17.markers <-FindMarkers(NCG_germ4_subset, ident.1 = 17)


# setwd("/Users/ewhelan/Desktop/NCG/NCG_marker_genes")
# write.csv(PGCs.E.markers.all, "PGCs.E.markers.all.csv")



genes.per.cluster <- 100


cluster0.markers.ordered <- cluster0.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster1.markers.ordered <- cluster1.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster2.markers.ordered <- cluster2.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster3.markers.ordered <- cluster3.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster4.markers.ordered <- cluster4.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster5.markers.ordered <- cluster5.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster6.markers.ordered <- cluster6.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster7.markers.ordered <- cluster7.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster8.markers.ordered <- cluster8.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster9.markers.ordered <- cluster9.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster10.markers.ordered <- cluster10.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster11.markers.ordered <- cluster11.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster12.markers.ordered <- cluster12.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster13.markers.ordered <- cluster13.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster14.markers.ordered <- cluster14.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster15.markers.ordered <- cluster15.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster16.markers.ordered <- cluster16.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster17.markers.ordered <- cluster17.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)

cluster0.markers.ordered$cluster <- 0
cluster1.markers.ordered$cluster <- 1
cluster2.markers.ordered$cluster <- 2
cluster3.markers.ordered$cluster <- 3
cluster4.markers.ordered$cluster <- 4
cluster5.markers.ordered$cluster <- 5
cluster6.markers.ordered$cluster <- 6
cluster7.markers.ordered$cluster <- 7
cluster8.markers.ordered$cluster <- 8
cluster9.markers.ordered$cluster <- 9
cluster10.markers.ordered$cluster <- 10
cluster11.markers.ordered$cluster <- 11
cluster12.markers.ordered$cluster <- 12
cluster13.markers.ordered$cluster <- 13
cluster14.markers.ordered$cluster <- 14
cluster15.markers.ordered$cluster <- 15
cluster16.markers.ordered$cluster <- 16
cluster17.markers.ordered$cluster <- 17


topNmarkers <- rbind(cluster10.markers.ordered, cluster15.markers.ordered, cluster13.markers.ordered,
                     cluster8.markers.ordered, cluster4.markers.ordered, cluster11.markers.ordered, cluster5.markers.ordered, 
                     cluster3.markers.ordered, cluster1.markers.ordered,cluster12.markers.ordered, cluster6.markers.ordered, 
                     cluster9.markers.ordered, cluster0.markers.ordered,  cluster2.markers.ordered, cluster7.markers.ordered, 
                     cluster14.markers.ordered, cluster16.markers.ordered, cluster17.markers.ordered)
                     
                     
                     
                     
               


# Heatmap of select genes by cluster ----

Idents(NCG_germ4_subset) <- NCG_germ4_subset$seurat_clusters
DimPlot(NCG_germ4_subset)

data.0 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 0)])
data.1 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 1)])
data.2 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 2)])
data.3 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 3)])
data.4 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 4)])
data.5 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 5)])
data.6 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 6)])
data.7 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 7)])
data.8 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 8)])
data.9 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 9)])
data.10 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 10)])
data.11<- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 11)])
data.12 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 12)])
data.13 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 13)])
data.14 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 14)])
data.15 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 15)])
data.16 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 16)])
data.17 <- as.matrix(GetAssayData(NCG_germ4_subset, slot = "data")[, WhichCells(NCG_germ4_subset, ident = 17)])

av.counts.0 <- apply(data.0, 1, mean)
av.counts.1 <- apply(data.1, 1, mean)
av.counts.2 <- apply(data.2, 1, mean)
av.counts.3 <- apply(data.3, 1, mean)
av.counts.4 <- apply(data.4, 1, mean)
av.counts.5 <- apply(data.5, 1, mean)
av.counts.6 <- apply(data.6, 1, mean)
av.counts.7 <- apply(data.7, 1, mean)
av.counts.8 <- apply(data.8, 1, mean)
av.counts.9 <- apply(data.9, 1, mean)
av.counts.10 <- apply(data.10, 1, mean)
av.counts.11 <- apply(data.11, 1, mean)
av.counts.12 <- apply(data.12, 1, mean)
av.counts.13 <- apply(data.13, 1, mean)
av.counts.14 <- apply(data.14, 1, mean)
av.counts.15 <- apply(data.15, 1, mean)
av.counts.16 <- apply(data.16, 1, mean)
av.counts.17 <- apply(data.17, 1, mean)

av.counts.HS_NCG_only_NCG <- cbind(av.counts.10, av.counts.15, av.counts.13, 
                                   av.counts.8, av.counts.4,av.counts.11, av.counts.5, av.counts.3,
                                   av.counts.1, av.counts.12, av.counts.6,av.counts.9, 
                                   av.counts.0, av.counts.2, 
                                    av.counts.7, 
                                     av.counts.14, 
                                   av.counts.16, av.counts.17)


av.counts.HS_NCG_only_NCG.df <- as.data.frame(av.counts.HS_NCG_only_NCG)

gene.name.list <- rownames(av.counts.HS_NCG_only_NCG.df)

av.counts.HS_NCG_only_NCG.df$gene.name <- gene.name.list

# genes.of.interest <- Keren.gene.list[c(1:48)] #this is just Keren's list of genes
# genes.of.interest <- (all.markers.subset$gene) #this is a log-fold cutoff
genes.of.interest <- rownames(topNmarkers) #this is just taking the top N markers

av.counts.genes.of.interest <- av.counts.HS_NCG_only_NCG.df %>%
  dplyr::filter(gene.name %in% genes.of.interest)


av.counts.genes.of.interest.ordered <- av.counts.genes.of.interest[match(genes.of.interest, av.counts.genes.of.interest$gene.name),]

integer.col <- ncol(av.counts.genes.of.interest.ordered)

av.counts.genes.of.interest.ordered <- av.counts.genes.of.interest.ordered[,-integer.col]

av.counts.genes.of.interest.matrix1 <- as.matrix(av.counts.genes.of.interest.ordered)
colnames(av.counts.genes.of.interest.matrix1) <- c("cluster10", "cluster15", "cluster13", 
  
  "cluster8", "cluster4", "cluster11", "cluster5", "cluster3", "cluster1", "cluster12", 
                                                   "cluster6", "cluster9", "cluster0", "cluster2", "cluster7", "cluster14", "cluster16", "cluster17")
  
av.counts.genes.of.interest.matrix.no.na <- na.omit(av.counts.genes.of.interest.matrix1)


l.mat <- rbind(c(0,4), c(0, 3), c(2,1), c(0,0))
l.mat

dev.off()
heatmap.2(av.counts.genes.of.interest.matrix.no.na, Colv = FALSE,
          Rowv = FALSE,
          trace = "none",
          # col =  colorRampPalette(colors=c("blue","white","red"))(100),
          
          # col = rev(brewer.pal(n = 11, name = "Spectral")),
          # col = rev(brewer.pal(n = 11, name = "RdYlGn")),
          # col = rev(brewer.pal(n = 11, name = "RdYlBu")),
          # col = rev(brewer.pal(n = 11, name = "RdGy")),
          col = rev(brewer.pal(n = 11, name = "RdBu")),
          # col = rev(brewer.pal(n = 11, name = "PuOr")),
          # col = rev(brewer.pal(n = 11, name = "PRGn")),
          # col = rev(brewer.pal(n = 11, name = "PiYG")),
          # col = rev(brewer.pal(n = 11, name = "BrBG")),
          # col = inferno(100),
          # col = viridis(100),
          density.info="none",
          srtCol = -90, #angle of labels
          adjCol = c(0,0), #offset of labels, should be centered now
          key = T,
          cexRow=1.0, cexCol=1.3,
          lmat= l.mat,
          lwid = c(0.5, 5),
          lhei=c(1, 1, 5, 1),
          scale = "row"
          
)



scaled_data <- scale(av.counts.genes.of.interest.matrix.no.na, center = TRUE, scale = TRUE)

# Create a heatmap with row scaling
heatmap(scaled_data, Rowv = NA, Colv = NA, col = cm.colors(256), scale = "row",
        margins = c(5,10), main = "Heatmap with Row Scaling")







# Reverse the order of rows in the scaled data
scaled_data <- scale(av.counts.genes.of.interest.matrix.no.na, center = TRUE, scale = TRUE)
scaled_data_reversed <- scaled_data[nrow(scaled_data):1, ]

heatmap(scaled_data_reversed, Rowv = NA, Colv = NA, col = viridis(256, option = "A"), scale = "row",
        margins = c(5, 30))

scaled_data_reversed <- as.data.frame(t(scale(t(av.counts.genes.of.interest.matrix.no.na), center = TRUE, scale = TRUE)))
scaled_data_reversed <- scaled_data[nrow(scaled_data):1, ]

# Convert row names into a column
scaled_data_reversed <- scaled_data_reversed %>%
  rownames_to_column(var = "Gene")

# Melt the data to long format
scaled_data_long <- scaled_data_reversed %>%
  pivot_longer(cols = -Gene, names_to = "Cluster", values_to = "Value")

gene_order <- rownames(av.counts.genes.of.interest.matrix.no.na)
cluster_order <- c("cluster10", "cluster15", "cluster13", 
                   
                   "cluster8", "cluster4", "cluster11", "cluster5", "cluster3", "cluster1", "cluster12", 
                   "cluster6", "cluster9", "cluster0", "cluster2", "cluster7", "cluster14", "cluster16", "cluster17")


scaled_data_long$Gene <- factor(scaled_data_long$Gene, levels =gene_order)
scaled_data_long$Cluster <- factor(scaled_data_long$Cluster, levels = cluster_order)

r_heatmap <- ggplot(scaled_data_long, aes(x=Cluster, y=Gene))+
  geom_tile(aes(fill=Value))+
  scale_fill_gradient2(low = "blue", mid = "white", high = "red")+
  # scale_fill_viridis(option = "magma")+
  # geom_point(aes(size=RSS_size), alpha=0.6)+
  ggtitle("RAT")+
  theme(axis.text.x = element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        legend.position = "none")
# theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
r_heatmap



genes <- c("POU5F1", "NANOG", "TCL1A", "PDPN", "PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1",
           "UTF1", "KIT", "ASB9", "NANOS2", "PIWIL4", "MORC1", "MAGEB2", "TDRD1", "RHOXF1", "MAGEC2", "ZBTB16", "EGR4", "RET", "MAGEA4",
           "SOHLH1", "DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8", "DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9", "ZCWPW1", "RAD51AP2", "SYCP3",
           "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1")

DimPlot(HS_NCG_only_NCG)
# NCG_germ4_subset <- subset(HS_NCG_only_NCG, idents = c("Pachytene spermatocytes", "Diplotene/2ndary spermatocytes", "Spermatids"), invert = T)
DimPlot(NCG_germ4_subset)
all.markers <- FindAllMarkers(NCG_germ4_subset, min.pct = 0, logfc.threshold = 0, features = genes)
PGCLC.markers <- FindMarkers(NCG_germ4_subset, ident.1 = "PGCLCs", min.pct = 0, logfc.threshold = 0, features = genes)
PGCs.E.markers <- FindMarkers(NCG_germ4_subset, ident.1 = "PGCs (early)", min.pct = 0, logfc.threshold = 0, features = genes)
PGCs.L.markers <- FindMarkers(NCG_germ4_subset, ident.1 = "PGCs (late)", min.pct = 0, logfc.threshold = 0, features = genes)
M.ProSpermatogonia.markers <- FindMarkers(NCG_germ4_subset, ident.1 = "M prospermatogonia", min.pct = 0, logfc.threshold = 0, features = genes)
M2T1.markers <- FindMarkers(NCG_germ4_subset, ident.1 = "M/T1 prospermatogonia", min.pct = 0, logfc.threshold = 0, features = genes)
T.ProSpermatogonia.markers <- FindMarkers(NCG_germ4_subset, ident.1 = "T1 prospermatogonia", min.pct = 0, logfc.threshold = 0, features = genes)
Spg.markers <- FindMarkers(NCG_germ4_subset, ident.1 = "Spermatogonia", min.pct = 0, logfc.threshold = 0, features = genes)
Diff.Spg.E.markers <- FindMarkers(NCG_germ4_subset, ident.1 = "Diff. spermatogonia (early)", min.pct = 0, logfc.threshold = 0, features = genes)
Diff.Spg.L.markers <- FindMarkers(NCG_germ4_subset, ident.1 = "Diff. spermatogonia (late)", min.pct = 0, logfc.threshold = 0, features = genes)
Pre.Leptotene.markers <- FindMarkers(NCG_germ4_subset, ident.1 = "Prelep/Lep/Zygotene", min.pct = 0, logfc.threshold = 0, features = genes)


PGCs.E.markers.ordered <- PGCs.E.markers %>% 
  arrange(avg_log2FC) %>% 
  slice_head(n = 50)



graph.data <- graph.data %>%
  arrange(Celltype)

order.of.genes <- as.vector(genes)


# graph.data$Pathway2 <- factor(graph.data$Pathway, levels = rev(unique(graph.data$Pathway)))
all.markers$gene <- factor(all.markers$gene, levels = rev(genes))
all.markers$cluster <- factor(all.markers$cluster, levels = c("PGCs (early)", 
                                                              "PGCs (late)", 
                                                              "M prospermatogonia",
                                                              "M/T1 prospermatogonia", 
                                                              "T1 prospermatogonia", 
                                                              "Spermatogonia",
                                                              "Diff. spermatogonia (early)", 
                                                              "Diff. spermatogonia (late)", 
                                                              "Prelep/Lep/Zygotene"))


ggplot(data = all.markers, na.rm = T,
       aes(x = cluster, y = (gene), 
           color = avg_log2FC, size = pct.1))+ 
  geom_point() +
  # scale_color_gradient(high = "red" , low = "blue") +
  scale_color_viridis(option = "inferno")+
  scale_size(range = c(-7, 7)) +
  theme_bw() + 
  ylab("Gene") + 
  xlab("Cell Type") + 
  ggtitle("Genes expressed by cell type") +
  scale_y_discrete(position = "right")+
  theme(axis.text.x = element_text( #face = "bold",
    color = "black", size = 14, angle = -60, hjust = 0
  ),
  axis.text.y = element_text( #face = "bold",
    color = "black",  hjust = 1)
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "left")




#pairwise comparison



all.markers$cluster <- factor(all.markers$cluster, levels = c("PGCs (early)", 
                                                              "PGCs (late)", 
                                                              "M prospermatogonia",
                                                              "M/T1 prospermatogonia", 
                                                              "T1 prospermatogonia", 
                                                              "Spermatogonia",
                                                              "Diff. spermatogonia (early)", 
                                                              "Diff. spermatogonia (late)", 
                                                              "Prelep/Lep/Zygotene"))

###
NCG_germ4_subset <- NCG_germ4
p.val.cutoff <- 0.05
log2fc.cutoff <- 0.25

###


DimPlot(NCG_germ4_subset)

PGCLC.vs.PGCE <- FindMarkers(NCG_germ4_subset, ident.1 = "PGCLCs", ident.2 = "PGCs (early)", logfc.threshold = log2fc.cutoff, min.pct = 0)
PGCE.vs.PGCL <- FindMarkers(NCG_germ4_subset, ident.1 = "PGCs (early)", ident.2 = "PGCs (late)", logfc.threshold = log2fc.cutoff, min.pct = 0)
PGCL.vs.M <- FindMarkers(NCG_germ4_subset, ident.1 = "PGCs (late)", ident.2 = "M prospermatogonia", logfc.threshold = log2fc.cutoff, min.pct = 0)
M.vs.MT1 <- FindMarkers(NCG_germ4_subset, ident.1 = "M prospermatogonia", ident.2 = "M/T1 prospermatogonia", logfc.threshold = log2fc.cutoff, min.pct = 0)
MT1.vs.T1 <- FindMarkers(NCG_germ4_subset, ident.1 = "M/T1 prospermatogonia", ident.2 = "T1 prospermatogonia", logfc.threshold = log2fc.cutoff, min.pct = 0)
T1.vs.Spg <- FindMarkers(NCG_germ4_subset, ident.1 = "T1 prospermatogonia", ident.2 = "Spermatogonia", logfc.threshold = log2fc.cutoff, min.pct = 0)
Spg.vs.DiffSpg <- FindMarkers(NCG_germ4_subset, ident.1 = "Spermatogonia", ident.2 = "Diff. spermatogonia (early)", logfc.threshold = log2fc.cutoff, min.pct = 0)
DiffSpg.vs.DiffSpg <- FindMarkers(NCG_germ4_subset, ident.1 = "Diff. spermatogonia (early)", ident.2 = "Diff. spermatogonia (late)", logfc.threshold = log2fc.cutoff, min.pct = 0)
DiffSpg.vs.Prelep <- FindMarkers(NCG_germ4_subset, ident.1 = "Diff. spermatogonia (late)", ident.2 = "Preleptotene spermatocytes", logfc.threshold = log2fc.cutoff, min.pct = 0)
# DiffSpg.vs.Prelep.all <- FindMarkers(NCG_germ4_subset, ident.1 = "Diff. spermatogonia (late)", ident.2 = "Prelep/Lep/Zygotene", logfc.threshold = 0, min.pct = 0)

setwd("~/Documents/scRNAseq/Rscripts/savefiles/xrTestis.cell.type.comparisons")
write.csv(PGCE.vs.PGCL, file = "PGCE.vs.PGCL.csv")
write.csv(PGCL.vs.M, file = "PGCL.vs.M.csv")
write.csv(M.vs.MT1, file = "M.vs.MT1.csv")
write.csv(MT1.vs.T1, file = "MT1.vs.T1.csv")
write.csv(T1.vs.Spg, file = "T1.vs.Spg.csv")
write.csv(Spg.vs.DiffSpg, file = "Spg.vs.DiffSpg.csv")
write.csv(DiffSpg.vs.DiffSpg, file = "DiffSpg.vs.DiffSpg.csv")
write.csv(DiffSpg.vs.Prelep, file = "DiffSpg.vs.Prelep.csv")


list.of.comparisons <- list(PGCLC.vs.PGCE, PGCE.vs.PGCL, PGCL.vs.M, M.vs.MT1, MT1.vs.T1, T1.vs.Spg, Spg.vs.DiffSpg, DiffSpg.vs.DiffSpg, DiffSpg.vs.Prelep)
list.of.counts <- list(av.counts.PGCLC, av.counts.pgcE, av.counts.pgcL, av.counts.m, av.counts.mt1, av.counts.t1, av.counts.undiff, av.counts.diffE, av.counts.diffL, av.counts.Spermatocytes)

list.of.cell.types <- c("PGCLCs", 
                        "PGCs (early)", 
                        "PGCs (late)", 
                        "M prospermatogonia",
                        "M/T1 prospermatogonia", 
                        "T1 prospermatogonia", 
                        "Spermatogonia",
                        "Diff. spermatogonia (early)", 
                        "Diff. spermatogonia (late)", 
                        "Preleptotene spermatocytes")

list.of.cell.types_abbr <- c("PGCLCs", 
                             "PGCs.E",
                             "PGCs.L",
                             "M",
                             "M.T1", 
                             "T1",
                             "Spg",
                             "DiffSpg.E",
                             "DiffSpg.L",
                             "Prelep")

# setwd("~/Desktop/list.of.genes")


for(number in 1:(length(list.of.cell.types)-1)){
  # number <- 1
  
  current.comparison <- list.of.comparisons[[number]]
  
  current.counts <- list.of.counts[[number]]
  # current.counts <- (log2(exp((current.counts))))
  current.counts2 <- list.of.counts[[number+1]]
  # current.counts2 <- (log2(exp((current.counts2))))
  current.comparison$gene <- rownames(current.comparison)
  # current.comparison.filtered <- filter(current.comparison, p_val_adj < p.val.cutoff)
  
  
  current.comparison.filtered.up <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC > log2fc.cutoff)
  current.comparison.filtered.down <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC < -log2fc.cutoff)
  current.comparison.filtered.significant <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff)
  
  current.comparison.dataframe <- data.frame(current.counts, current.counts2)
  # current.counts.df <- data.frame(current.counts, current.counts2)
  current.comparison.dataframe$log2fc <- current.counts - current.counts2
  current.comparison.dataframe$regulation <- "not"
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc > log2fc.cutoff] <- "up"
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc < -log2fc.cutoff] <- "down"
  current.comparison.dataframe$regulation[!(current.comparison.dataframe$gene.name %in% current.comparison.filtered.significant$gene)] <- "not"
  

  # #simply deduct count values
  # current.comparison.dataframe$difference <- current.comparison.dataframe$current.counts - current.comparison.dataframe$current.counts2
  # 
  
  
  
  current.comparison.subset <- dplyr::filter(current.comparison.dataframe, current.counts > 0 | current.counts2 >0)
  # current.comparison.sig.up <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.up$gene)
  # current.comparison.sig.down <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.down$gene)
  current.comparison.subset <- na.omit(current.comparison.subset)
  current.comparison.subset$DEGs <- "Not significant"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene] <- "Up"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene] <- "Down"
  print(paste0(list.of.cell.types[number], " vs ", list.of.cell.types[number+1], 
               " Up: ", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]), 
               " Down:", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])))
  
  print(table(current.comparison.subset$regulation))
  
  
  # write.csv(current.comparison.subset, file=paste0(list.of.cell.types_abbr[number], ".csv"))
  
  # assign(paste0("p", number),
  #        (ggplot() + geom_point(data = current.comparison.subset, aes(x = current.counts, y = current.counts2,  color = DEGs)) +
  #                                scale_color_manual(values=c('#E69F00','#999999', '#56B4E9')) +
  #                                theme_classic() + labs(title=paste(list.of.cell.types[number], " vs. ", list.of.cell.types[number+1], ", q", p.val.cutoff, ", log2fc", log2fc.cutoff),
  #                                                       x=list.of.cell.types[number], y = list.of.cell.types[number+1]) + theme(plot.title = element_text(size=9)))
  # )
  assign(paste0("p", number),
         (ggplot() + geom_point(data = current.comparison.subset, aes(x = current.counts, y = current.counts2,  color = regulation)) +
            scale_color_manual(values=c('#E69F00','#999999', '#56B4E9')) +
            theme_classic() + labs(title=paste(list.of.cell.types[number], " vs. ", list.of.cell.types[number+1], ", q", p.val.cutoff, ", log2fc", log2fc.cutoff),
                                   x=list.of.cell.types[number], y = list.of.cell.types[number+1]) + theme(plot.title = element_text(size=9)))
  )

}


plot_grid(ncol = 3, p1, p2, p3, p4, 
          p5, p6, p7, p8, p9)





# gene.of.interest <- "STRA8"

# gene.of.interest.PGCE <- c("SPRR2F", "ABHD12B", "HMGA2", "DDT", "KHDC3L", "NANOG", "SOX17", "YBX1", "MIF", "TFAP2C")
# # gene.of.interest.PGCL <- c("DDX4", "DAZL", "DCAF4L1","GAGE1", "TEX14", "GAGE12H", "GAGE2A", "XIST", "DCAF4L1", "ELAVL2", "TFRAP2C", "UTF1",
# #                           "CCND1", "ID3", "MAGEB2", "FHL2", "ID1", "MAGED2", "TFRAP2C", "SOX17", "NANOG", "UTF1", "POU5F1", "SOX15")
# # gene.of.interest.PGCL <- c("TFAP2C", "NANOG", "SOZ17", "POU5F1")
# gene.of.interest.M <- c("GCNA", "HIST1H3B", "FOXM1", "SYCP2", "E2F8", "CDC25A", "XIST", "MKI67","GCNA","TEX15", "SPC25","CDC25A",
#                         "ACTG2", "DPPA5", "UBE2C", "KHDC3L", "GPHA2", "CENPA", "NANOG", "APOC1", "POU5F1")
# # gene.of.interest.M <- c("UTF1")
# gene.of.interest.transitional <- c("NANOS2", "PAGE2", "RHOXF1", "ASB9", "TEX15", "TDRD1", "PAGE2B",
#                                    "DGKK", "NANOS3", "MDK", "MAP4K4", "STMN1", "IGFBP2")
# gene.of.interest.T1 <- c("FOXF1", "DPPA5", "EGR4", "SIX1", "PIWIL4", "VCX3B", "MAGEC2", "DPPA2", "MORC1", "MAGEB2",
#                           "GAGE12H", "NANOS2", "MEG3", "EGR4", "XIST", "HELLPAR", "DCAF4L1")
# gene.of.interest.spg <- c("RET", "POU3F1", "HMGB4", "TSPAN33", "FOXD1", "SCT", "KCNQ2", "MAGEA4", "VCX", "TCF3", "UTF1", "PIWIL4", 
#                           "UTF1", "MORC1", "SCT", "FOXD1", "POU3F1", "RGS14", "TSPAN33", "PIWIL4", "EGR4", "MAPK3", "HES5", "PHOX2A", "FIGLA")
# gene.of.interest.earlydiff <- c("HIST1H1A", "L1TD1", "ASB9", "EIF3K", "EGR1", "HIST1H3B", "HMGB4",
#                                 "RESF1", "ID3", "DUSP6", "EGR4", "SIX1", "DPPA2", "FGFR3", "ID4", "ID2", "PCSK1N", "ASB9", "MAGEA4")
# gene.of.interest.latediff <- c("TEX19", "MEIOC", "SYCP3", "SYCE2", "RAD51AP2", "HORMAD1", "SCML1", "CENPU", "TOP2A", "SMC1B", "HSPA5", "HIST1H1A",
#                                "ELAVL2", "CENPF", "TDRD1", "ESX1", "ITGA6", "RHOXF1", "KIT", "UBE2C", "HIST1H2AB", "SOHLH1", "CCND1")
# gene.of.interest.prelep <- c("TEX101", "MEIOB", "SPO11", "ZCWPW1", "RAD51AP2", "SYCP1", "SYCP2", "TEX12", "CCNB3", "RNF212B", "FBXO47", "HOTAIR")


# gene.of.interest.PGCE <- c("SPRR2F", "ABHD12B", "HMGA2", "DDT", "KHDC3L", "NANOG", "SOX17", "YBX1", "MIF", "TFAP2C")
# # gene.of.interest.PGCL <- c("DDX4", "DAZL", "DCAF4L1","GAGE1", "TEX14", "GAGE12H", "GAGE2A", "XIST", "DCAF4L1", "ELAVL2", "TFRAP2C", "UTF1",
# #                           "CCND1", "ID3", "MAGEB2", "FHL2", "ID1", "MAGED2", "TFRAP2C", "SOX17", "NANOG", "UTF1", "POU5F1", "SOX15")
# # gene.of.interest.PGCL <- c("TFAP2C", "NANOG", "SOZ17", "POU5F1")
# gene.of.interest.M <- c("GCNA", "HIST1H3B", "FOXM1", "SYCP2", "E2F8", "CDC25A", "XIST", "MKI67","GCNA","TEX15", "SPC25","CDC25A",
#                         "ACTG2", "DPPA5", "UBE2C", "KHDC3L", "GPHA2", "CENPA", "NANOG", "APOC1", "POU5F1")
# # gene.of.interest.M <- c("UTF1")
# gene.of.interest.transitional <- c("NANOS2", "PAGE2", "RHOXF1", "ASB9", "TEX15", "TDRD1", "PAGE2B",
#                                    "DGKK", "NANOS3", "MDK", "MAP4K4", "STMN1", "IGFBP2")
# gene.of.interest.T1 <- c("FOXF1", "DPPA5", "EGR4", "SIX1", "PIWIL4", "VCX3B", "MAGEC2", "DPPA2", "MORC1", "MAGEB2",
#                          "GAGE12H", "NANOS2", "MEG3", "EGR4", "XIST", "HELLPAR", "DCAF4L1")
# gene.of.interest.spg <- c("RET", "POU3F1", "HMGB4", "TSPAN33", "FOXD1", "SCT", "KCNQ2", "MAGEA4", "VCX", "TCF3", "UTF1", "PIWIL4", 
#                           "UTF1", "MORC1", "SCT", "FOXD1", "POU3F1", "RGS14", "TSPAN33", "PIWIL4", "EGR4", "MAPK3", "HES5", "PHOX2A", "FIGLA")
# gene.of.interest.earlydiff <- c("HIST1H1A", "L1TD1", "ASB9", "EIF3K", "EGR1", "HIST1H3B", "HMGB4",
#                                 "RESF1", "ID3", "DUSP6", "EGR4", "SIX1", "DPPA2", "FGFR3", "ID4", "ID2", "PCSK1N", "ASB9", "MAGEA4")
gene.of.interest.PGCLC <- c("TFAP2C", #"SOX17", 
                            # "YBX1", #"MIF", 
                            "POU5F1", "NANOS3", "DND1", #"PDPN", 
                            "ID3", #"LIN28A", 
                            # "ID2", 
                            "EPCAM", #"CXCR4",# "DPPA4", 
                            "MSX2", "ID1", "NANOG", #"KRT8", 
                            "KRT18" , #"CDH1", "SOX15", 
                            "KRT19",
                            "DNMT1")
gene.of.interest.PGCE1 <- c(#"HORMAD1", 
  "ACTG2", #"ETV5", 
  "KIT", #"PIWIL2", 
  "DDIT4", #"UTF1", #"THY1",
                            "TCL1B", "TCL1A", "TDRD12" #"PRAME"
                            )
gene.of.interest.PGCE2 <- c("TFAP2C", "SOX17", "YBX1", "MIF", "POU5F1")
gene.of.interest.PGCL1 <- c("DDX4", "DAZL", "MKI67",  "TOP2A", "ID4")
gene.of.interest.PGCL2 <- c("ACTG2",  "POU5F1", "S100A13", "CCNB1", "PCSK1N" )
gene.of.interest.M1 <- c("DAZL", "ID4", "DDX43", "DDX4", "E2F8")
gene.of.interest.M2 <- c("NANOG", "ACTG2", "MKI67", "NANOS3", "POU5F1", "STMN1", "CENPA")
gene.of.interest.transitional1 <- c("MDK", "ASB9", "CITED2", "MAP4K4", "TDRD1", "ID3")
gene.of.interest.transitional2 <- c("MDK", "ASB9", "NANOS3", "POU5F1", "TSPAN6", "APOE", "STMN1", "MEG3")
gene.of.interest.T1 <- c("NANOS2", "RHOXF1", "TDRD1", "MAGEB2", "TEX15", "DCAF4L1", "MAGEC2", "PAGE2B",
                         "SIX1", "DDX4", "FOXP1", "TCF3", "VCX", 'EGR4', "VCX3B", "DPPA2")
gene.of.interest.T1_2 <- c("MEG3", "NANOS2", "TPM2")
gene.of.interest.spg1 <- c("PIWIL4", "SIX1", "HES5", "MAGEB2",  "TCF3", "UTF1", "EGR4")
gene.of.interest.spg2 <- c("DPPA2", "SIX1", "FGFR3", "PIWIL4", "EGR4", "ID4", "SOHLH1", "UTF1")
gene.of.interest.earlydiff1 <- c("ASB9", "STRA8","PRSS21", "SYCP3", "MLEC", "PARP4", "REC8", "LARP1", "DNMT1", "PIWIL1", "TOP2B")
gene.of.interest.earlydiff2 <- c("KIT", "PIWIL1", "ASB9", "ETV7", "L1TD1", "SCML1", "MKI67", "TEX15", "HORMAD1", "MEIOC", "SYCP3")
gene.of.interest.latediff1 <- c("SYCP3", "SMC1B", "MEIOC", "MAGEA4", "UBE2C", "ESX1", "MKI67")
gene.of.interest.latediff2 <- c("CENPF", "KIT", "ITGA6", "UBE2C", "ELAVL2",  "RHOXF1")
gene.of.interest.prelep1 <- c("TEX101", "MEIOB", "HORMAD1", "SYCP2", "SYCP1", "ZCWPW1","SCML1", "TOP2A","SYCP3", "SPO11", "PRDM9")

gene.list.r0 <- c(gene.of.interest.PGCLC, gene.of.interest.PGCE1)
gene.list.r1 <- c(gene.of.interest.PGCE2, gene.of.interest.PGCL1)
gene.list.r2 <- c(gene.of.interest.PGCL2, gene.of.interest.M1)
gene.list.r3 <- c(gene.of.interest.M2, gene.of.interest.transitional1)
gene.list.r4 <- c(gene.of.interest.transitional2, gene.of.interest.T1)
gene.list.r5 <- c(gene.of.interest.T1_2, gene.of.interest.spg1)
gene.list.r6 <- c(gene.of.interest.spg2, gene.of.interest.earlydiff1)
gene.list.r7 <- c(gene.of.interest.earlydiff2, gene.of.interest.latediff1)
gene.list.r8 <- c(gene.of.interest.latediff2, gene.of.interest.prelep1)




# list.of.lists.of.genes <- list(gene.of.interest.PGCE,
#                             gene.of.interest.PGCL,
#                             gene.of.interest.M,
#                             gene.of.interest.transitional,
#                             gene.of.interest.T1,
#                             gene.of.interest.spg,
#                             gene.of.interest.earlydiff,
#                             gene.of.interest.latediff,
#                             gene.of.interest.prelep)

list.of.lists.of.genes <- list(gene.list.r0,
                               gene.list.r1,
                               gene.list.r2,
                               gene.list.r3,
                               gene.list.r4,
                               gene.list.r5,
                               gene.list.r6,
                               gene.list.r7,
                               gene.list.r8)



#CELL TYPE COMPARISON CHARTS
#starting with PGC-E
cell.type.colors <- c("#a6cee2", "#2179b4", "#b4d88b", "#36a047", "#f6999a", "#f15a29", "#fdbf6f", "#f7941d", "#cab3d6", "#6b3f98")


log2fc.cutoff <- 0.25

for(number in 1:(length(list.of.cell.types)-1)){
  # number <- 1
  print(number)
  # gene.of.interest <- c(list.of.lists.of.genes[[number]], list.of.lists.of.genes[[number+1]])
  gene.of.interest <- list.of.lists.of.genes[[number]]
  # gene.of.interest <- NULL
  current.comparison <- list.of.comparisons[[number]]
  current.counts <- list.of.counts[[number]]
  current.counts2 <- list.of.counts[[number+1]]
  
  current.comparison$gene <- rownames(current.comparison)
  # current.comparison.filtered <- filter(current.comparison, p_val_adj < p.val.cutoff)
  current.comparison.filtered.up <- filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC > 0)
  current.comparison.filtered.down <- filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC < 0)
  current.comparison.filtered.significant <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff)
  

  
  
  current.comparison.dataframe <- data.frame(current.counts, current.counts2)
  current.comparison.dataframe$log2fc <- current.counts - current.counts2
  current.comparison.dataframe$regulation <- "Nonsig"
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc > log2fc.cutoff] <- "Up"
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc < -log2fc.cutoff] <- "Down"
  
  current.comparison.dataframe$regulation[!(current.comparison.dataframe$gene.name %in% current.comparison.filtered.significant$gene)] <- "Nonsig"
  
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  # current.comparison.subset <- filter(current.comparison.dataframe, av.counts.diffL > 0 | av.counts.Spermatocytes>0)
  current.comparison.subset <- dplyr::filter(current.comparison.dataframe, current.counts > 0 | current.counts2 >0)
  # current.comparison.sig.up <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.up$gene)
  # current.comparison.sig.down <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.down$gene)
  current.comparison.subset <- na.omit(current.comparison.subset)
  current.comparison.subset$DEGs <- "Not significant"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene] <- "Up"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene] <- "Down"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% gene.of.interest] <- "GeneOfInterest"
  current.comparison.subset$regulation[current.comparison.subset$gene.name %in% gene.of.interest] <- "GeneOfInterest"
  print(paste0(list.of.cell.types[number], " vs ", list.of.cell.types[number+1], 
               " Up: ", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]), 
               " Down:", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])))
  
  colnames(current.comparison.subset) <- c("HS.counts", "NCG.counts", "log2FC_manual", "regulation", "gene.name", "DEGs")
  print(table(current.comparison.subset$regulation))

  degs_order <- c( "Down", "GeneOfInterest", "Not significant", "Up")
  
  
  current.comparison.subset_ordered <- current.comparison.subset %>%
    mutate(DEGs = factor(DEGs, levels = degs_order)) 
  
  setwd("~/Documents/scRNAseq/Rscripts/savefiles/xrTestis.cell.type.comparisons")
  write.csv(current.comparison.subset_ordered, file=paste0(list.of.cell.types_abbr[number], "_vs_", list.of.cell.types_abbr[number+1], "_manual_comparison.csv"))

  assign(paste0("r", number), 

         ggplot(data = current.comparison.subset_ordered, aes(x = HS.counts, y = NCG.counts,
                                                                                   label = gene.name,
                                                                                   color = regulation)) +
           geom_point(data = filter(current.comparison.subset_ordered, regulation == "Nonsig"), alpha = 0.5) + # Plot "Not Significant" points first
           geom_point(data = filter(current.comparison.subset_ordered, regulation != "Nonsig")) + # Plot other points
           geom_text(data = filter(current.comparison.subset_ordered, regulation != "Nonsig")) +
           geom_point(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest")) +
           geom_text(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest"), hjust = -0.1, vjust = -0.1) +


           scale_color_manual(values = c(cell.type.colors[number+1], "black", '#999999', cell.type.colors[number])) +
           theme_classic() +
           xlim(0, 3.5) + ylim(0, 3.5) +
           labs(title = paste(list.of.cell.types[number], " vs. ", list.of.cell.types[number + 1], ", q", p.val.cutoff, ", log2fc", log2fc.cutoff),
                x = list.of.cell.types[number], y = list.of.cell.types[number + 1]) +
           theme(plot.title = element_text(size = 9),
                 panel.grid.major = element_line(color = "gray", linetype = "dotted"),
                 panel.grid.minor = element_line(color = "lightgray", linetype = "dotted")) +
           geom_abline(intercept = 0, slope = 1, linetype = "dotted")
         )
}

plot_grid(ncol = 3, r1, r2, r3, r4, 
          r5, r6, r7, r8, r9) 

r6
r2
r3
r4
r5
r6
r7
r8




#detect cells with the matching barcode found in fastq files
#repeat using the in-vivo aligned data
cells.NCG <- Cells(HS_NCG_only_NCG)


#TGCGGGTAGTCTACCA

cell.1 <- "CCACGTTGTTATGTCG"
# cell.2 <- "TAACACGCACCATTCC"
# cell.3 <- "TGCGGGTAGTCTACCA"
cell.4 <- "GGGAGATTCGCAGTTA"
# cell.5 <- "GTCGCGAGTCTGTGAT"
# cell.6 <- "CGTAGTATCAACACCA"
cell.7 <- "AAGCGTTGTGCAATAA"
cell.8 <- "CAACCAATCCCTTGGT" #cell8
table(str_detect(cells.NCG, cell.8))

cell.of.interest1 <- cells.NCG[str_detect(cells.NCG, cell.1)]
# cell.of.interest3 <- cells.NCG[str_detect(cells.NCG, cell.3)]
cell.of.interest4 <- cells.NCG[str_detect(cells.NCG, cell.4)]
cell.of.interest7 <- cells.NCG[str_detect(cells.NCG, cell.7)]
cell.of.interest8 <- cells.NCG[str_detect(cells.NCG, cell.8)]

cells.of.interest.list <- c(cell.of.interest1, 
                            # cell.of.interest3, 
                            cell.of.interest4, cell.of.interest7, cell.of.interest8)

DimPlot(HS_NCG_only_NCG, cells.highlight = cells.of.interest.list) + NoLegend()


FeaturePlot(HS_NCG_only_NCG,
            features = "MEIOB", label = F, repel = FALSE) +
  scale_color_viridis(option="plasma") + DarkTheme()



###PSEUDOTIME STARTS HERE###
#this time using the in-vivo-aligned data

HS_NCG_only_NCG2 <- HS_NCG_only_NCG


cds <- as.cell_data_set(HS_NCG_only_NCG2)
cds <- as.cell_data_set(HS_NCG_only_NCG2)


cds <- cluster_cells(cds, reduction_method = c("UMAP"), #cluster_method = "louvain",
                     resolution = 5e-5, 
                     k = 8,
                     partition_qval = 0.05
)

plot_cells(cds, show_trajectory_graph = FALSE)
cds <- learn_graph(cds, 
                   close_loop = F,
                   learn_graph_control=list(ncenter=300, minimal_branch_len=7), #250, 18 is pretty good
                   use_partition = F)
print(plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)) #+ ggtitle(a))

# Figure 2C #######################################
cds <- order_cells(cds)#, root_cells = max.avp)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)+ NoLegend()


cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(HS_NCG_only_NCG2[["RNA"]])
cds_subset <- cds["SYCP3",]
plot_genes_in_pseudotime(cds_subset)

pseudotime.values <- pseudotime(cds, reduction_method = "UMAP")


HS_NCG_only_NCG2[["pseudotime"]] <- pseudotime.values


# c("SYCP1", "MEIOB", "MEIOC", "SPO11")
#"STRA8", "EGR4", "ESX1", "MAGEC2"

FeaturePlot(HS_NCG_only_NCG2, features = c("SYCP1", "MEIOB", "MEIOC", "SPO11"))  & scale_color_viridis(option="magma") & DarkTheme()
DimPlot(HS_NCG_only_NCG2, group.by = "experiment")
DimPlot(HS_NCG_only_NCG2, group.by = "orig.ident")
DimPlot(HS_NCG_only_NCG2, group.by = "seurat_clusters", label = T)
Idents(HS_NCG_only_NCG2) <- "seurat_clusters"
cluster.16 <- subset(HS_NCG_only_NCG2, idents = 16)
Idents(cluster.16) <- "orig.ident"
table(Idents(cluster.16))

setwd("~/Documents/scRNAseq/Rscripts/savefiles")
save(HS_NCG_only_NCG2, file = "HS_NCG_only_NCG2_2023_05_12.Robj")


### heatmap ###


genes <- c("POU5F1", "NANOG", "TCL1A", "PDPN", "PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1",
           "UTF1", "KIT", "ASB9", "NANOS2", "PIWIL4", "MORC1", "MAGEB2", "TDRD1", "RHOXF1", "MAGEC2", "ZBTB16", "EGR4", "RET", "MAGEA4",
           "SOHLH1", "DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8", "DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9", "ZCWPW1", "RAD51AP2", "SYCP3",
           "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1")

genes <- rownames(topNmarkers)

variable.features <- VariableFeatures(NCG_germ4_subset)

# "MLH3", "PIWIL1", "MLH1", "SPATA8", "CCDC112", "OVOL1", "OVOL2",
# "CCDC42", "ACR", "CDRT15", "TEKT1", "H1FOO", "TJP3", "FAM24A", "CATSPER3", "SPACA3", "SPACA1", "ACRV1", "SAXO1", "CAPZA3", "HEMGN", "CA2", "PRM3",
# "TNP2", "TSSK6", "SMCP", "SPEM2", "LELP1", "MKI67", "UBE2C", "TOP2A")

HS_NCG_only_NCG2 <- NCG_germ4

DefaultAssay(HS_NCG_only_NCG2) <- "RNA"


#can I do this with seurat objects? Need to pull out the normalized counts

pt.matrix.all <- HS_NCG_only_NCG2@assays$RNA@data
pseudotime <- as.data.frame(HS_NCG_only_NCG2$pseudotime)

dim(pt.matrix.all)
# gene.names.all <- rownames(pt.matrix.all)

# length(genes)
# length(intersect(genes, gene.names.all))
# setdiff(genes, gene.names.all)

pt.matrix.all[1:4, 1:4]




intersect.genes <- (intersect(genes, rownames(pt.matrix.all)))
genes <- intersect.genes
# genes <- variable.features

pt.matrix.subset <- pt.matrix.all[genes,]
t.pt.matrix.all <- t(pt.matrix.subset)
t.pt.matrix.all[1:4, 1:4]


pseudotime[1:4,1]
t.pt.matrix.all <- as.data.frame(t.pt.matrix.all)
dim(t.pt.matrix.all)
t.pt.matrix.all$pseudotime <- pseudotime[,1]
dim(t.pt.matrix.all)
t.pt.matrix.all[1:5, 1:5]
t.pt.matrix.sorted <- t.pt.matrix.all[order(t.pt.matrix.all$pseudotime),]
t.pt.matrix.sorted[1:5, ]


# pt.matrix <- t(as.matrix(t.pt.matrix.sorted[,1:(length(genes))]))
pt.matrix <- t(as.matrix(t.pt.matrix.sorted))
pt.matrix[, 1:5]

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=5)$y}))
# pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=2)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
dim(pt.matrix)
pt.matrix <- pt.matrix[c(1:length(genes)),]
dim(pt.matrix)


#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 16),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  use_raster                   = FALSE,
  cluster_rows                 = TRUE,
    cluster_row_slices           = TRUE,
  cluster_columns              = FALSE)

print(hthc)

VlnPlot(NCG_germ4, features = c("NANOS2", "ID1", "MEIOB"))



#plot pseudotime values by sample

pseudotime.values <- NCG_germ4$pseudotime
range((pseudotime.values), na.rm = TRUE)

orig.ident <- NCG_germ4$orig.ident
stage.values <- NCG_germ4$stage

pseudotime.df <- data.frame(stage.values, pseudotime.values)

stripplot(pseudotime.values)


table(stage.values)


ggplot(pseudotime.df, aes(x=pseudotime.values, y=0, color = stage.values)) +
  geom_point(size = 4)  +
  # annotate("segment",x=1,xend=20, y=0, yend=0, size=2) +
  # annotate("segment",x=1,xend=1, y=-0.1,yend=0.1, size=2) +
  # annotate("segment",x=20,xend=20, y=-0.1,yend=0.1, size=2) +
  # geom_text(aes(label = x), col="white") +
  scale_x_continuous(limits = c(1,max(pseudotime.df$pseudotime.values))) +
  scale_y_continuous(limits = c(-1,1)) +
  # scale_color_manual(values = unname(colours)) #+ 
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())



q_colors =  12
v_colors =  viridis(q_colors, option = "B")
show_col(v_colors)

stripchart(pseudotime.values ~ stage.values, data = pseudotime.df, col = v_colors[1:9], 
           pch = 16, frame = FALSE, method = "jitter")


hist(pseudotime.df$pseudotime.values, breaks = 100)

#Try with density lines

table(pseudotime.df$stage.values)

# 
# c41_120.values <-subset(pseudotime.df,stage.values=='c41+120',select=('pseudotime.values'))
# c91_d180.values <-subset(pseudotime.df,stage.values=='c91+d180',select=('pseudotime.values'))
# d120.values <-subset(pseudotime.df,stage.values=='d120',select=('pseudotime.values'))
# d123.values <-subset(pseudotime.df,stage.values=='d123',select=('pseudotime.values'))
# d127.values <-subset(pseudotime.df,stage.values=='d127',select=('pseudotime.values'))
# d180.values <-subset(pseudotime.df,stage.values=='d180',select=('pseudotime.values'))
# d181.values <-subset(pseudotime.df,stage.values=='d181',select=('pseudotime.values'))
# d245.values <-subset(pseudotime.df,stage.values=='d245',select=('pseudotime.values'))
# d251.values <-subset(pseudotime.df,stage.values=='d251',select=('pseudotime.values'))
# 
# plot(density(c41_120.values$pseudotime.values),col='red')
# lines(density(c91_d180.values$pseudotime.values),col='blue')
# 
# 
# plot(density(d120.values$pseudotime.values), col = "orange1")
# lines(density(d123.values$pseudotime.values),col='orange2')
# lines(density(d127.values$pseudotime.values),col='orange3')
# lines(density(d180.values$pseudotime.values),col='red1')
# lines(density(d181.values$pseudotime.values),col='red2')
# lines(density(d245.values$pseudotime.values),col='red3')
# lines(density(d251.values$pseudotime.values),col='red4')
# 
# #Lets try in ggplot2
# 
colour.values <- c("khaki1", "khaki1", "red1", "red1", "red1", "green", "green", "blue", "blue" )
# colour.values <- c("orange", "orange1", "red1", "red2", "red3", "red4", "darkred", "chocolate4", "brown" )



a <- ggplot(pseudotime.df, aes(x = pseudotime.values))
a + geom_density(aes(fill = stage.values, alpha = 0.4))+
  # scale_color_manual(values = ) +
  scale_fill_manual(values = colour.values)


pseudotime.df_transplant <- filter(pseudotime.df, stage.values %in% c("d120", "d123", "d127", 
                                                                   "d180", "d181", "d245", "d251"))


ggplot(pseudotime.df_transplant, aes(x = pseudotime.values, y = stage.values, fill = ..x..)) +
  geom_density_ridges_gradient(scale = 3, rel_min_height = 0.01) +
  scale_fill_viridis(name = "pseudotime", option = "inferno") +
  xlim(0, 45) +
  theme_ridges() + 
  theme(legend.position = "none")


Idents(NCG_germ4) <- "stage"
DimPlot(NCG_germ4)
NCG_germ4_transplant <- subset(NCG_germ4, idents = c("hPGCLCs", "c91+d180", "c41+120"), invert = T)
NCG_germ4_culture <- subset(NCG_germ4, idents = c("c91+d180", "c41+120"), invert = F)
NCG_germ4_noPGCLC <- subset(NCG_germ4, idents = "hPGCLCs", invert = T)
# FeaturePlot(NCG_germ4, features = "pseudotime", split.by = "stage")
FeaturePlot(NCG_germ4, features = "pseudotime", split.by = "stage", combine = TRUE) & NoLegend() & 
  NoAxes() &
  scale_color_viridis(option = "plasma") & DarkTheme() + theme(axis.line.x.bottom=element_line(color="black"), axis.line.y.left = element_line(color = "black"))


DimPlot(NCG_germ4_noPGCLC, group.by = "cell.type")

pt.matrix.ncg <- NCG_germ4_noPGCLC@assays$RNA@data
pseudotime <- as.data.frame(NCG_germ4_noPGCLC$pseudotime)
celltype.names <- as.data.frame(NCG_germ4_noPGCLC$cell.type)

dim(pt.matrix.ncg)

Keren.gene.list <- c("POU5F1", "NANOG", "TCL1A", "PDPN", "PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1","DDX4",
                     "UTF1", "KIT", "ASB9", "NANOS2", "PIWIL4", "MORC1", "MAGEB2","MAGEA3", "TDRD1", "RHOXF1", "MAGEC2", "GFRA1", "ZBTB16", "EGR4", "RET", "MAGEA4",
                     "SOHLH1", "DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8", "DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9", "ZCWPW1", "RAD51AP2", "SYCP3",
                     "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1", "MLH3", "PIWIL1", "MLH1", "SPATA8", "CCDC112", "OVOL1", "OVOL2",
                     "CCDC42", "ACR", "CDRT15", "TEKT1", "H1FOO", "TJP3", "FAM24A", "CATSPER3", "SPACA3", "SPACA1", "ACRV1", "SAXO1", "CAPZA3", "HEMGN", "CA2", "PRM3",
                     "TNP2", "TSSK6", "SMCP", "SPEM2", "LELP1", "MKI67", "UBE2C", "TOP2A")

genes <- c(Keren.gene.list)
genes <- c(#"DNMT1", "UHRF1", "DNMT3L", "DNMT3A", "DNMT3B" 
           # "TET1", "TET2", "TET3", "AICDA"
           # "EHMT1", "EHMT2","SETDB1", "SUV39H1", "SUV39H2"
           # "EZH1", "EZH2", "SUZ12", "EED"
           "KDM1A", "KDM3A", "KDM3B", "KDM3A", "KDM4B", "KDM4C", "PHF8"
           )

genes <- c("DNMT1", "DNMT3A", "TET1", "EHMT1", 
           "SETDB1", "SUV39H2", "SUZ12", "EZH2",
           "PHF8", "KDM3B")

all.genes <- row.names(pt.matrix.ncg)

genes_not_in_all <- genes[!genes %in% all.genes]

genes <- intersect(proliferation.genes[-330], all.genes)

genes <- unique(c(genes, "TOP2A", "MKI67", "CENPF", "UBE2C"))

# Output genes not present in all.genes
print(genes_not_in_all)


pt.matrix.subset <- data.frame(pt.matrix.ncg[genes,])


t.pt.matrix.all <- t(pt.matrix.subset)
t.pt.matrix.all[1:2, 1:2]

t.pt.matrix.all <- as.data.frame(t.pt.matrix.all)
dim(t.pt.matrix.all)
t.pt.matrix.all$pseudotime <- pseudotime[,1]
t.pt.matrix.all$cell.type <- celltype.names[,1]
dim(t.pt.matrix.all)
# t.pt.matrix.all[1:5, 1:5]
t.pt.matrix.sorted <- t.pt.matrix.all[order(t.pt.matrix.all$pseudotime),]
t.pt.matrix.sorted[1:5, ]


#pseudotime line graphs for figure 3.
###FIGURE 3###


# Assuming t.pt.matrix.sorted is your data frame
# Define gene names excluding pseudotime and cell.type columns
gene_names <- names(t.pt.matrix.sorted)[1:(ncol(t.pt.matrix.sorted) - 2)]

# Melt the data frame to long format

t.pt.matrix.long <- melt(t.pt.matrix.sorted, id.vars = c("pseudotime", "cell.type"))

# Plot using ggplot
p <- ggplot(t.pt.matrix.long, aes(x = pseudotime, y = value, color = variable)) +
  geom_smooth(method = "lm", formula = y ~ splines::bs(x, df = 4), se = FALSE) +
  scale_color_manual(values = rainbow(length(gene_names))) +  # Specify colors
  # scale_y_continuous(limits = c(0, 4)) +
  theme_minimal() +
  theme(legend.position = "right",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        axis.ticks = element_line(color = "black")
  )

print(p)





#with cell values

ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=POU5F1, color = cell.type)) + geom_point()



# Reshape the data to long format
gene_data_long <- tidyr::pivot_longer(t.pt.matrix.sorted, cols = c(TFAP2C, DDX4, NANOG), names_to = "Gene", values_to = "Expression")

# Plot the curves
p1 <- ggplot(gene_data_long, aes(x = pseudotime, y = Expression, color = Gene)) +
  geom_smooth(se = FALSE) +
  labs(title = "Gene Expression Curves", x = "Pseudotime", y = "Expression") +
  scale_color_manual(values = c("TFAP2C" = "green", "DDX4" = "red", "NANOG" = "cyan")) +
  theme(
    axis.line = element_line(color = "black", size = 1),  # Set axis lines to black and thick
    panel.grid.major = element_blank(),  # Remove background grid
    panel.grid.minor = element_blank(),
    panel.background = element_blank()  # Remove background color
  )

gene_data_long <- tidyr::pivot_longer(t.pt.matrix.sorted, cols = c(TFAP2C, DDX4, KIT), names_to = "Gene", values_to = "Expression")

p2 <- ggplot(gene_data_long, aes(x = pseudotime, y = Expression, color = Gene)) +
  geom_smooth(se = FALSE) +
  labs(title = "Gene Expression Curves", x = "Pseudotime", y = "Expression") +
  scale_color_manual(values = c("TFAP2C" = "green", "DDX4" = "red", "KIT" = "cyan")) +
  theme(
    axis.line = element_line(color = "black", size = 1),  # Set axis lines to black and thick
    panel.grid.major = element_blank(),  # Remove background grid
    panel.grid.minor = element_blank(),
    panel.background = element_blank()  # Remove background color
  )

gene_data_long <- tidyr::pivot_longer(t.pt.matrix.sorted, cols = c(MAGEC2, DDX4, GFRA1), names_to = "Gene", values_to = "Expression")

p3 <- ggplot(gene_data_long, aes(x = pseudotime, y = Expression, color = Gene)) +
  geom_smooth(se = FALSE) +
  labs(title = "Gene Expression Curves", x = "Pseudotime", y = "Expression") +
  scale_color_manual(values = c("MAGEC2" = "green", "DDX4" = "red", "GFRA1" = "cyan")) +
  theme(
    axis.line = element_line(color = "black", size = 1),  # Set axis lines to black and thick
    panel.grid.major = element_blank(),  # Remove background grid
    panel.grid.minor = element_blank(),
    panel.background = element_blank()  # Remove background color
  )

gene_data_long <- tidyr::pivot_longer(t.pt.matrix.sorted, cols = c(MAGEC2, DDX4, GFRA1), names_to = "Gene", values_to = "Expression")

p3 <- ggplot(gene_data_long, aes(x = pseudotime, y = Expression, color = Gene)) +
  geom_smooth(se = FALSE) +
  labs(title = "Gene Expression Curves", x = "Pseudotime", y = "Expression") +
  scale_color_manual(values = c("MAGEC2" = "green", "DDX4" = "red", "GFRA1" = "cyan")) +
  theme(
    axis.line = element_line(color = "black", size = 1),  # Set axis lines to black and thick
    panel.grid.major = element_blank(),  # Remove background grid
    panel.grid.minor = element_blank(),
    panel.background = element_blank()  # Remove background color
  )


gene_data_long <- tidyr::pivot_longer(t.pt.matrix.sorted, cols = c(TFAP2C, DDX4, UTF1), names_to = "Gene", values_to = "Expression")

p4 <- ggplot(gene_data_long, aes(x = pseudotime, y = Expression, color = Gene)) +
  geom_smooth(se = FALSE) +
  labs(title = "Gene Expression Curves", x = "Pseudotime", y = "Expression") +
  scale_color_manual(values = c("TFAP2C" = "green", "DDX4" = "red", "UTF1" = "cyan")) +
  theme(
    axis.line = element_line(color = "black", size = 1),  # Set axis lines to black and thick
    panel.grid.major = element_blank(),  # Remove background grid
    panel.grid.minor = element_blank(),
    panel.background = element_blank()  # Remove background color
  )


gene_data_long <- tidyr::pivot_longer(t.pt.matrix.sorted, cols = c(SYCP3, DDX4, MAGEA3), names_to = "Gene", values_to = "Expression")

p5 <- ggplot(gene_data_long, aes(x = pseudotime, y = Expression, color = Gene)) +
  geom_smooth(se = FALSE) +
  labs(title = "Gene Expression Curves", x = "Pseudotime", y = "Expression") +
  scale_color_manual(values = c("SYCP3" = "green", "DDX4" = "red", "MAGEA3" = "cyan")) +
  theme(
    axis.line = element_line(color = "black", size = 1),  # Set axis lines to black and thick
    panel.grid.major = element_blank(),  # Remove background grid
    panel.grid.minor = element_blank(),
    panel.background = element_blank()  # Remove background color
  )


plot_grid(p1, p2, p3, p4, p2, p5, ncol = 1)




#ACTUAL COLOURS FOR UMAP FROM ILLUSTRATOR 
cols = c(#"#a6cee2", 
  "#2179b4",  "#b4d88b", "#36a047", "#f6999a", "#f15a29", "#fdbf6f", "#f7941d", "#cab3d6", "#6b3f98")

p1 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=POU5F1, color = cell.type)) + geom_point() + scale_color_manual(values = cols) + theme(legend.position = "none")
p2 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=NANOG, color = cell.type)) + geom_point() + scale_color_manual(values = cols) + theme(legend.position = "none")
p3 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=UTF1, color = cell.type)) + geom_point()  + scale_color_manual(values = cols)+ theme(legend.position = "none")
p4 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=KIT, color = cell.type)) + geom_point()  + scale_color_manual(values = cols)+ theme(legend.position = "none")
p5 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=NANOS2, color = cell.type)) + geom_point() + scale_color_manual(values = cols) + theme(legend.position = "none")
p6 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=SOHLH1, color = cell.type)) + geom_point()  + scale_color_manual(values = cols)+ theme(legend.position = "none")
p7 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=MAGEB2, color = cell.type)) + geom_point()  + scale_color_manual(values = cols)+ theme(legend.position = "none")
p8 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=DMRT1, color = cell.type)) + geom_point() + scale_color_manual(values = cols) + theme(legend.position = "none")
p9 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=STRA8, color = cell.type)) + geom_point()  + scale_color_manual(values = cols)+ theme(legend.position = "none")
p10 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=DAZL, color = cell.type)) + geom_point() + scale_color_manual(values = cols)+ theme(legend.position = "none")
p11 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=PRDM9, color = cell.type)) + geom_point() + scale_color_manual(values = cols)+ theme(legend.position = "none")
p12 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=SYCP1, color = cell.type)) + geom_point() + scale_color_manual(values = cols) + theme(legend.position = "none")
p13 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=ACRV1, color = cell.type)) + geom_point() + scale_color_manual(values = cols) + theme(legend.position = "none")

### Extended Data Fig 3f

plot_gene_expression <- function(gene_name) {
  cols = c(#"#a6cee2", 
    "#2179b4",  "#b4d88b", "#36a047", "#f6999a", "#f15a29", "#fdbf6f", "#f7941d", "#cab3d6", "#6b3f98")
  
  p <- ggplot(t.pt.matrix.sorted, aes(x = pseudotime, y = !!sym(gene_name), color = cell.type)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ splines::bs(x, df = 4), se = FALSE, color = "black") +  # Add spline line
    scale_color_manual(values = cols) +
    scale_y_continuous(limits = c(0, 4)) +  # Set y-axis limits
    theme_minimal() +  # Remove background color
    theme(legend.position = "none",
          panel.grid.major = element_blank(),  # Remove gridlines
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),  # Remove background color
          axis.line = element_line(color = "black"),  # Set axis lines to black and solid
          axis.ticks = element_line(color = "black")  # Set tick marks to black
    )
  
  return(p)
}



p1 <- plot_gene_expression("POU5F1")
p2 <- plot_gene_expression("NANOG")
p3 <- plot_gene_expression("UTF1")
p4 <- plot_gene_expression("KIT")
p5 <- plot_gene_expression("NANOS2")
p6 <- plot_gene_expression("SOHLH1")
p7 <- plot_gene_expression("MAGEB2")
p8 <- plot_gene_expression("DMRT1")
p9 <- plot_gene_expression("STRA8")
p10 <- plot_gene_expression("DAZL")
p11 <- plot_gene_expression("PRDM9")
p12 <- plot_gene_expression("SYCP1")
p13 <- plot_gene_expression("ACRV1")
p14 <- plot_gene_expression("MAGEC2")

# plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, p13, nrow = 3)

plot_grid(p1, p4, p5, p7, p12, p13, nrow = 2)


###


DimPlot(NCG_germ4, split.by = "stage", group.by = "cell.type") &
  # scale_color_viridis(option="plasma") & 
  DarkTheme()

#extended data fig 3h
#work out the proportions by sample. 
DimPlot(NCG_germ4, group.by = "stage")
individual.samples <- SplitObject(NCG_germ4, split.by = "stage")

for(current.loop in 1:length(individual.samples)){
  # current.loop <- 1
  Idents(individual.samples[[current.loop]]) <- "stage"
  print(table(Idents(individual.samples[[current.loop]])))
  Idents(individual.samples[[current.loop]]) <- "cell.type"
  print(table(Idents(individual.samples[[current.loop]])))
}




# Figure showing proportions of samples.



DimPlot(NCG_germ4_noPGCLC)
DimPlot(NCG_germ4_transplant, split.by = "stage")



FeaturePlot(NCG_germ4_noPGCLC, features = c("TFAP2C", "DDX4", "KIT"), ncol = 3)  & NoLegend() & NoAxes() & DarkTheme() & scale_color_viridis(option = "magma")

DimPlot(NCG_germ4, label = T, group.by = "cell.type") + NoLegend() + NoAxes()
DimPlot(NCG_germ4, label = T, group.by = "cell.type")+ NoAxes()




Idents(NCG_germ4) <- "stage"
DimPlot(NCG_germ4)
NCG_germ4_4m <- subset(NCG_germ4, idents = c("d123", "d127", "d120"))
NCG_germ4_6m <- subset(NCG_germ4, idents = c("d180", "d181"))
NCG_germ4_8m <- subset(NCG_germ4, idents = c("d251", "d245"))
Idents(NCG_germ4_4m) <- "cell.type"
Idents(NCG_germ4_6m) <- "cell.type"
Idents(NCG_germ4_8m) <- "cell.type"
p1 <- DimPlot(NCG_germ4_4m) + NoAxes() + NoLegend()
p2 <- DimPlot(NCG_germ4_6m) + NoAxes() + NoLegend()
p3 <- DimPlot(NCG_germ4_8m) + NoAxes() + NoLegend()
plot_grid(p1, p2, p3, ncol = 3)


# ViolinPlots

Idents(NCG_germ4) <- "cell.type"


VlnPlot(NCG_germ4, features = c("MAGEC2", "DDX4", "GFRA1"), group.by = "cell.type")
# 
# ###

###


plot_gene_expression <- function(genes_to_plot) {
  # Extract expression values for the specified genes
  gene_expression <- GetAssayData(NCG_germ4, slot = "data")[genes_to_plot, ]
  
  # Combine cell type information and gene expression data
  gene_data <- data.frame(CellType = NCG_germ4$cell.type)
  for (gene in genes_to_plot) {
    gene_data[gene] <- gene_expression[gene, ]
  }
  
  # Pivot the data to long format
  gene_data_long <- tidyr::pivot_longer(gene_data, cols = -CellType, names_to = "Gene", values_to = "Expression")
  
  # Create a custom order for genes
  gene_data_long$Gene <- factor(gene_data_long$Gene, levels = genes_to_plot)
  
  # Define custom colors for each gene
  # gene_colors <- c("#6abd45", "#97282c", "#57bfc5")
  gene_colors <- c("#97282c", "#57bfc5")
  
  ggplot(gene_data_long, aes(x = CellType, y = Expression, fill = Gene)) +
    geom_violin(position = "dodge", scale = "width", trim = FALSE) +
    # labs(title = "Gene Expression Violin Plot", x = "Cell Type", y = "Expression") +
    scale_fill_manual(values = gene_colors) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),  # Remove x-axis labels
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_line(color = "black"),  # Set axis lines to black and solid
      panel.border = element_blank()  # Remove panel border
    ) +
    scale_y_continuous(limits = c(0, NA))  # Set y-axis limit with 0 as the minimum
}


# Now you can call this function with any set of genes to update the plot


### ### ### ### ### ###
###Figure 3B (right) ###
### ### ### ### ### ###

Idents(NCG_germ4) <- "cell.type"
DimPlot(NCG_germ4)

genes_to_plot <- c("TFAP2C", "DDX4", "NANOG")
FeaturePlot(NCG_germ4, features = genes_to_plot)
plot_gene_expression(genes_to_plot)
DimPlot(NCG_germ4, cols = c("magenta3","magenta3","grey","grey","grey","grey","grey","grey","grey","grey")) + NoAxes() + NoLegend()

genes_to_plot <- c("TFAP2C", "DDX4", "KIT")
FeaturePlot(NCG_germ4, features = genes_to_plot)
DimPlot(NCG_germ4, cols = c("grey","grey","magenta3","magenta3","magenta3","grey","grey","grey","grey","grey")) + NoAxes() + NoLegend()
plot_gene_expression(genes_to_plot)

genes_to_plot <- c("MAGEC2", "DDX4", "GFRA1")
FeaturePlot(NCG_germ4, features = genes_to_plot)
DimPlot(NCG_germ4, cols = c("grey","grey","grey","grey","grey","magenta3","magenta3","grey","grey","grey")) + NoAxes() + NoLegend()
plot_gene_expression(genes_to_plot)

genes_to_plot <- c("MAGEC2", "DDX4", "UTF1")
FeaturePlot(NCG_germ4, features = genes_to_plot)
DimPlot(NCG_germ4, cols = c("grey","grey","grey","grey","grey","grey","magenta3","grey","grey","grey")) + NoAxes() + NoLegend()
plot_gene_expression(genes_to_plot)


genes_to_plot <- c("TFAP2C", "DDX4", "KIT")
FeaturePlot(NCG_germ4, features = genes_to_plot)
DimPlot(NCG_germ4, cols = c("grey","grey","grey","grey","grey","grey","grey","magenta3","magenta3","grey")) + NoAxes() + NoLegend()
plot_gene_expression(genes_to_plot)

genes_to_plot <- c("SYCP3", "DDX4", "MAGEA3")
FeaturePlot(NCG_germ4, features = genes_to_plot)
DimPlot(NCG_germ4, cols = c("grey","grey","grey","grey","grey","grey","grey","grey","grey","magenta3")) + NoAxes() + NoLegend()
plot_gene_expression(genes_to_plot)



genes_to_plot <- c("MEIOB", "SYCP1")
FeaturePlot(NCG_germ4, features = genes_to_plot)
DimPlot(NCG_germ4, cols = c("grey","grey","grey","grey","grey","grey","grey","grey","grey","magenta3")) + NoAxes() + NoLegend()
plot_gene_expression(genes_to_plot)






# Extract the expression values for the specified genes
gene_expression <- GetAssayData(NCG_germ4, slot = "data")[c("DDX4", "TFAP2C"), ]

# Create a data frame with the gene expression data and cell metadata
gene_data <- data.frame(
  CellType = NCG_germ4$cell.type,
  DDX4 = gene_expression["DDX4", ],
  TFAP2C = gene_expression["TFAP2C", ]
)

# Plot a scatter plot of gene expression
ggplot(gene_data, aes(x = TFAP2C, y = DDX4, color = CellType)) +
  geom_point(position = position_jitter(width = 0.1, height = 0.1)) +
  labs(title = "Expression of TFAP2C vs DDX4", x = "TFAP2C expression", y = "DDX4 expression") +
  theme_minimal() +
  theme(
    axis.line = element_line(color = "black"),  # Set axis lines to black
    panel.background = element_blank(),  # Remove background color
    panel.grid.major = element_blank(),  # Remove gridlines
    panel.grid.minor = element_blank()  # Remove gridlines
  )



eoin.palette.light <- c("antiquewhite3",
  "lightcyan3", "pink3", "brown",  "sienna1", "aquamarine3","dodgerblue1", "yellowgreen", "lightgoldenrod2","darkorchid1" )

eoin.palette.dark <- c("gray28",
                        "cyan4", "pink3", "brown",  "sienna3", "aquamarine3","dodgerblue3", "darkgreen", "mediumpurple3","darkorchid1" )




### FIGURE 2 ### 

DimPlot(NCG_germ4, cols = eoin.palette.dark, group.by = "cell.type")

### Figure 2B ###
DefaultAssay(NCG_germ4) <- "RNA"

p2 <- FeaturePlot(NCG_germ4, features = c("DDX4"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "red", "orangered", "yellow", "white"))
p1 <- FeaturePlot(NCG_germ4, features = c("TFAP2C"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "darkgreen", "green", "yellow", "white"))
p3 <- FeaturePlot(NCG_germ4, features = c("PIWIL4"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "cyan4", "cyan", "white"))

plot_grid(p1, p2, p3, ncol = 3)

p4 <- FeaturePlot(NCG_germ4, features = c("MAGEB2"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "blue2", "deepskyblue", "skyblue"))
p5 <- FeaturePlot(NCG_germ4, features = c("SOHLH2"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "magenta4", "hotpink", "white"))
p6 <- FeaturePlot(NCG_germ4, features = c("SYCP1"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "brown", "burlywood"))


### FIGURE 2A###

cell.type.colors <- c("#a6cee2", "#2179b4", "#b4d88b", "#36a047", "#f6999a", "#f15a29", "#fdbf6f", "#f7941d", "#cab3d6", "#6b3f98")
Idents(NCG_germ4) <- "cell.type"
DimPlot(NCG_germ4, cols = cell.type.colors)


setwd("~/Desktop/proliferation_genes/")

genes<-c("TOP2A", "MKI67", "CENPF", "UBE2C")

for(current.gene in 1:length(genes)){
  # current.gene <- 1
  print(FeaturePlot(NCG_germ4, features = genes[current.gene], order = F, min.cutoff = 0) & DarkTheme() & scale_colour_viridis(option = "C"))
  ggsave(filename = paste0("featureplot_", genes[current.gene], ".png"))
  print(plot_gene_expression(genes[current.gene])+ggtitle(genes[current.gene]))
  ggsave(filename = paste0("pseudotimeplot_", genes[current.gene], ".png"))
  print(VlnPlot(NCG_germ4, features = genes[current.gene], cols = cell.type.colors, pt.size = 0))
  ggsave(filename = paste0("violinplot_", genes[current.gene], ".png"))
}

FeaturePlot(NCG_germ4, features = genes, order = F, min.cutoff = 0) & DarkTheme() & scale_colour_viridis(option = "C")

###Figure 2D###

Stacked_VlnPlot(NCG_germ4, features = genes, x_lab_rotate = TRUE,
                # theme(text=element_text(family="Comic Sans MS")),
                colors_use = cell.type.colors,
                pt.size = 0)
#banana
p1 <-   print(plot_gene_expression(genes[1]))
p2 <-   print(plot_gene_expression(genes[2]))
p3 <-   print(plot_gene_expression(genes[3]))
p4 <-   print(plot_gene_expression(genes[4]))

plot_grid(p1, p2, p3, p4, ncol = 1)

#Subset out the spermatogonia and recluster



DimPlot(NCG_germ4)
spermatogonia_xrTestis <- subset(NCG_germ4, idents = "Spermatogonia")
spermatogonia_xrTestis<- FindVariableFeatures(object = spermatogonia_xrTestis, nfeatures = 2000, verbose = FALSE)
# DefaultAssay(object = current.integrated) <- "integrated"
spermatogonia_xrTestis <- ScaleData(object = spermatogonia_xrTestis, vars.to.regress = c("nCount_RNA", "percent.mito"), verbose = FALSE)
spermatogonia_xrTestis <- RunPCA(object = spermatogonia_xrTestis,
                             npcs = 30,
                             verbose = TRUE)
spermatogonia_xrTestis <- FindNeighbors(object = spermatogonia_xrTestis)
spermatogonia_xrTestis <- FindClusters(object = spermatogonia_xrTestis,
                                   resolution = 1.1,
                                   #algorithm = 1,
                                   verbose = TRUE
)

spermatogonia_xrTestis <- RunUMAP(object = spermatogonia_xrTestis, reduction = "pca",
                              dims = 1:30)

DimPlot(spermatogonia_xrTestis)
FeaturePlot(spermatogonia_xrTestis, features = c("ID4", "NANOS2", "TCN2", "EGR4", "PLPPR5", "TSPAN33"
                                                 ), ncol = 3, order = T) & DarkTheme() & scale_color_viridis(option = "C")



#These are Seurat's default list of cell cycle genes

s.genes <- (cc.genes.updated.2019$s.genes)
g2m.genes <- (cc.genes.updated.2019$g2m.genes)

#Cell cycle scoring:
DefaultAssay(spermatogonia_xrTestis) <- "RNA"
spermatogonia_xrTestis <- CellCycleScoring(spermatogonia_xrTestis, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


fp1 <- FeaturePlot(spermatogonia_xrTestis, features = c("S.Score"), min.cutoff = 0) + labs(title = "S score")
fp2 <- FeaturePlot(spermatogonia_xrTestis, features = c("G2M.Score"), min.cutoff = 0) + labs(title = expression(bold(paste(G[2],"-M score",sep=""))))
wrap_plots(fp1, fp2)

DimPlot(spermatogonia_xrTestis, group.by = "Phase")



cell.type.colors <- c("#a6cee2", "#2179b4", "#b4d88b", "#36a047", "#f6999a", "#f15a29", "#fdbf6f", "#f7941d", "#cab3d6", "#6b3f98")
DimPlot(NCG_germ4, split.by = "stage",  cols = cell.type.colors)


#MAKE A PSEUDOTIME COLOR PLOT

n <- 39.05160*10

image(
  1:n, 1, as.matrix(1:n),
  col = viridis(n, option = "inferno"),
  xlab = "viridis", ylab = " ", 
  xaxt = "n",  # Remove default x-axis ticks
  yaxt = "n",  # Remove default y-axis ticks
  bty = "n"    # Remove box around the plot
)

# Add custom x-axis ticks
axis(1, at = seq(0, 40*10, by = 50), labels = seq(0, 40, by = 5))




### redo the individual cell plot marking the cells of interest ###
### FIGURE 3F & G ###

#detect cells with the matching barcode found in fastq files
cells.NCG <- Cells(NCG_germ4)

cell.1 <- "CCACGTTGTTATGTCG"
# cell.2 <- "TAACACGCACCATTCC"
# cell.3 <- "TGCGGGTAGTCTACCA"
cell.4 <- "GGGAGATTCGCAGTTA"
# cell.5 <- "GTCGCGAGTCTGTGAT"
# cell.6 <- "CGTAGTATCAACACCA"
cell.7 <- "AAGCGTTGTGCAATAA"
cell.8 <- "CAACCAATCCCTTGGT" #cell8
table(str_detect(cells.NCG, cell.8))

cell.of.interest1 <- cells.NCG[str_detect(cells.NCG, cell.1)]
# cell.of.interest3 <- cells.NCG[str_detect(cells.NCG, cell.3)]
cell.of.interest4 <- cells.NCG[str_detect(cells.NCG, cell.4)]
cell.of.interest7 <- cells.NCG[str_detect(cells.NCG, cell.7)]
cell.of.interest8 <- cells.NCG[str_detect(cells.NCG, cell.8)]

cells.of.interest.list <- c(cell.of.interest1,
                            cell.of.interest4, cell.of.interest7, cell.of.interest8)

p1 <- DimPlot(NCG_germ4, cells.highlight = cell.of.interest1, size = 2) & ggtitle(cell.of.interest1)
p2 <- DimPlot(NCG_germ4, cells.highlight = cell.of.interest4, size = 2) & ggtitle(cell.of.interest4)
p3 <- DimPlot(NCG_germ4, cells.highlight = cell.of.interest7, size = 2) & ggtitle(cell.of.interest7)
p4 <- DimPlot(NCG_germ4, cells.highlight = cell.of.interest8, size = 2) & ggtitle(cell.of.interest8)

plot_grid(p1, p2, p3, p4)

q1 <- DimPlot(NCG_germ4, cells.highlight = cells.of.interest.list, size = 1, cols = c("gray", "red"))  & NoLegend()
q2 <- FeaturePlot(NCG_germ4, features = "MEIOB", cols = c("gray",  "red"))
q3 <- FeaturePlot(NCG_germ4, features = c("MEIOB"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "red", "orangered", "yellow", "white"))
q4 <- FeaturePlot(NCG_germ4, features = c("MEIOB"), order = F, min.cutoff = 0) & scale_colour_gradientn(colors = c("gray","red", "darkred"))
plot_grid(q1, q2, q3, q4)
plot_grid(q1, q4)





### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### PART 2          ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### In vivo human   ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###



### Read in all human in vivo samples ###



HS51_hu <- Read10X(data.dir = "/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/Hs51/filtered_feature_bc_matrix")
HS30H_hu <- Read10X_h5("/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/HS30H/filtered_feature_bc_matrix.h5")
HS27_hu <- Read10X_h5("/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/HS27/filtered_feature_bc_matrix.h5")
HS26_hu <- Read10X_h5("/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/HS26/filtered_feature_bc_matrix.h5")
HS25_hu <- Read10X(data.dir = "/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/HS25/filtered_feature_bc_matrix")
HS22H_hu <- Read10X(data.dir = "/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/HS22H/filtered_feature_bc_matrix")
HS18H_hu <- Read10X(data.dir = "/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/HS18H/filtered_feature_bc_matrix")
HS14H_hu <- Read10X(data.dir = "/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/HS14H/filtered_feature_bc_matrix")
HS12H_hu <- Read10X(data.dir = "/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/Hs12H/filtered_feature_bc_matrix")
HS3S_hu <- Read10X(data.dir = "/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/HS3S/filtered_feature_bc_matrix")

HS46_hu <- Read10X_h5("/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/HS46/filtered_feature_bc_matrix.h5")
HS43_hu <- Read10X_h5("/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/HS43/filtered_feature_bc_matrix.h5")
HS41FT_hu <- Read10X_h5("/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/HS41FT/filtered_feature_bc_matrix.h5")
HS40_hu <- Read10X_h5("/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/HS40/filtered_feature_bc_matrix.h5")
HS31_hu <- Read10X_h5("/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/HS31/filtered_feature_bc_matrix.h5")

HS49_1_hu <- Read10X_h5("/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/HS49_1/filtered_feature_bc_matrix.h5")
HS26H_hu <- Read10X_h5("/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/HS26H/filtered_feature_bc_matrix.h5")

HS578_multi.data <- Read10X_h5("/Users/ewhelan/Documents/PROJECTS/ksasaki/scRNAseq/Keren/Hs578_multi/filtered_feature_bc_matrix.h5")
HS578_multi_hu <- HS578_multi.data$`Gene Expression`


human.list2 <- c(HS51_hu, HS30H_hu, HS27_hu, HS26_hu, HS25_hu, HS22H_hu, HS18H_hu, HS14H_hu, HS12H_hu, HS3S_hu, HS46_hu, HS43_hu, HS41FT_hu, HS40_hu, HS31_hu,
                 HS49_1_hu, HS26H_hu, xrTestis6M_RAinj_hu) #, HS578_multi_hu)

names.data2 <- c("HS51", "HS30H", "HS27", "HS26", "HS2S", "HS22H", "HS18H", "HS14H", "HS12H", "HS3S", "HS46", "HS43", "HS41FT", "HS40", "HS31", 
                 "HS49_1", "HS26H") #, "HS578_multi")

for(a in 1:length(human.list2)){
  print(length(Cells(human.list2[[a]])))
}

inflection.point2 <- c(1:length(human.list2))
log_lib_size_at_inflection2 <- c(1:length(human.list2))


par(mfrow = c(3,3))

for(y in 1:length(human.list2)){
  #y=1
  print(names.data2[y])
  expression.df <- as.data.frame(human.list2[y])
  
  umi_per_barcode <- colSums(expression.df)
  barcode_rank <- rank(-umi_per_barcode)
  # plot(barcode_rank, umi_per_barcode, 
  #      #ylim=c(1,2000),
  #      xlim=c(1,10000))
  
  log_lib_size <- log10(umi_per_barcode)
  #plot(barcode_rank, log_lib_size, xlim=c(1,10000))
  
  o <- order(barcode_rank)
  log_lib_size <- log_lib_size[o]
  barcode_rank <- barcode_rank[o]
  
  rawdiff <- diff(log_lib_size)/diff(barcode_rank)
  inflection <- which(rawdiff == min(rawdiff[500:4000], na.rm=TRUE))
  
  
  
  plot(barcode_rank, log_lib_size, xlim=c(1,8000),
       pch = 20,
       main = paste(y, names.data2[y], "cells", inflection, sep="_", "UMIs", round(10^log_lib_size[inflection]))
  )
  
  abline(v=inflection, col="blue", lwd=2)
  abline(h=log_lib_size[inflection], col="green", lwd=2)
  
  
  inflection.point2[y] <- inflection
  log_lib_size_at_inflection2[y] <- log_lib_size[inflection]
}

min.UMIs2 <- 10^log_lib_size_at_inflection2

#manually set minimum umi numbers

min.UMIs2[1] <- 10^3.7 #HS51
min.UMIs2[2] <- 10^3.9 #HS30H
min.UMIs2[3] <- 10^3.5 #HS27
min.UMIs2[4] <- 10^3.8 #HS26
min.UMIs2[5] <- 10^3.6 #HS25
min.UMIs2[6] <- 10^3.7 #HS22H
min.UMIs2[7] <- 10^4 #HS18H
min.UMIs2[8] <- 10^3.6 #HS14H
min.UMIs2[9] <- 10^3 #HS12H
min.UMIs2[10] <- 10^3.5 #HS3S
min.UMIs2[11] <- 10^3.6 #HS46
min.UMIs2[12] <- 10^3.5 #HS43
min.UMIs2[13] <- 10^3.6 #HS41FT
min.UMIs2[14] <- 10^3.3 #HS40
min.UMIs2[15] <- 10^3.6 #HS31

min.UMIs2[16] <- 10^3.8#"HS49_1", 
min.UMIs2[17] <- 10^3.8#"HS26H", 
# min.UMIs2[18] <- 10^3.8#"HS578_multi"

#other cutoffs

min.features <- 100
max.features <- 3000
max.mito <- 20

HS51 <- CreateSeuratObject(HS51_hu, project = "HS51")
HS30H <- CreateSeuratObject(HS30H_hu, project = "HS30H")
HS27 <- CreateSeuratObject(HS27_hu, project = "HS27")
HS26 <- CreateSeuratObject(HS26_hu, project = "HS26")
HS25 <- CreateSeuratObject(HS25_hu, project = "HS2S")
HS22H <- CreateSeuratObject(HS22H_hu, project = "HS22H")
HS18H <- CreateSeuratObject(HS18H_hu, project = "HS18H")
HS14H <- CreateSeuratObject(HS14H_hu, project = "HS14H")
HS12H <- CreateSeuratObject(HS12H_hu, project = "HS12H")
HS3S <- CreateSeuratObject(HS3S_hu, project = "HS3S")
HS3S <- CreateSeuratObject(HS3S_hu, project = "HS3S")
HS46 <- CreateSeuratObject(HS46_hu, project = "HS46")
HS43 <- CreateSeuratObject(HS43_hu, project = "HS43")
HS41FT <- CreateSeuratObject(HS41FT_hu, project = "HS41FT")
HS40 <- CreateSeuratObject(HS40_hu, project = "HS40")
HS31 <- CreateSeuratObject(HS31_hu, project = "HS31")

HS49_1 <- CreateSeuratObject(HS49_1_hu, project = "HS49_1")
HS26H <- CreateSeuratObject(HS26H_hu, project = "HS26H")
HS578_multi <- CreateSeuratObject(HS578_multi_hu, project = "HS578_multi")

list.of.seurats <- c( HS51, # "12 yo"
                      HS30H, # "GW7
                      HS27, # "GW18 + 5d"
                      HS26, # "GW8 + 4d"
                      HS25, # "5 yo"
                      HS22H, # "GW7
                      HS18H, # "GW6"
                      HS14H, # "GW8"
                      HS12H, # "GW10 + 5d"
                      HS3S, # "8 yo"
                      HS46, # "GW10"
                      HS43, # "GW14 + 2d"
                      HS41FT, # "GW37"
                      HS40, # "GW37 + 4m"
                      HS31, # "GW17 + 2d"
                      HS49_1, #GW30 + 10d
                      HS26H)#, #GW8 + 4d
# HS578_multi) #26w6d


stage.data <-  c("12 yo",
                 "GW7",
                 "GW18 + 5d",
                 "GW8 + 4d",
                 "5 yo",
                 "GW7",
                 "GW6",
                 "GW8",
                 "GW10 + 5d",
                 "8 yo",
                 "GW10",
                 "GW14 + 2d",
                 "GW37",
                 "GW37 + 4m",
                 "GW17 + 2d",
                 "GW30 + 10d",
                 "GW8 + 4d")

# HS578_multi) #26w6d


for(current.sample in 1:length(list.of.seurats)){
  # current.sample <- 1
  current.seurat <- list.of.seurats[[current.sample]]
  current.seurat$replicate <- names.data2[current.sample]
  current.seurat$experiment <- "human in vivo"
  current.seurat$stage <- stage.data[current.sample]
  current.seurat[['percent.mito']] <- PercentageFeatureSet(current.seurat, pattern = "^MT-")
  list.of.seurats[[current.sample]] <- current.seurat
}

for(current.sample in 1:(length(list.of.seurats))){
  # current.sample <- 5
  current.seurat <- list.of.seurats[[current.sample]]
  current.seurat <- subset(x = current.seurat, 
                           # subset = nFeature_RNA > min.features & nFeature_RNA < max.features & percent.mito < max.mito & nCount_RNA > min.UMIs[current.sample])
                           subset = percent.mito < max.mito & nCount_RNA > min.UMIs2[current.sample]) # 
  
  current.seurat <- NormalizeData(current.seurat)
  current.seurat <- CellCycleScoring(current.seurat, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)
  list.of.seurats[[current.sample]] <- current.seurat
}

number.of.samples <- length(list.of.seurats)
number.of.cells.per.sample <-  data.frame(
  SampleName = c(1:number.of.samples),
  CellNumber = c(1:number.of.samples),
  GeneNumber = c(1:number.of.samples),
  SampleNumber = c(1:number.of.samples),
  SampleOrigin = c(1:number.of.samples),
  medianUMI = c(1:number.of.samples),
  medianGenes = c(1:number.of.samples),
  medianMito =c(1:number.of.samples),
  cutoff = c(1:number.of.samples)
)

for(a in 1:length(list.of.seurats)){
  number.of.cells.per.sample[a,4]<-a
  number.of.cells.per.sample[a,1]<-names.data2[a]
  number.of.cells.per.sample[a,2]<-length(Cells(list.of.seurats[[a]]))
  number.of.cells.per.sample[a,3]<-length(rownames(list.of.seurats[[a]]))
  number.of.cells.per.sample[a,5]<-stage.data[a]
  number.of.cells.per.sample[a,6]<-median(list.of.seurats[[a]]$nCount_RNA)
  number.of.cells.per.sample[a,7]<-median(list.of.seurats[[a]]$nFeature_RNA)
  number.of.cells.per.sample[a,8]<-median(list.of.seurats[[a]]$percent.mito)
  number.of.cells.per.sample[a,9]<-min.UMIs2[a]
}

print(number.of.cells.per.sample)

Cells.merged <- merge(
  x = list.of.seurats[[1]],
  y = list.of.seurats[2:length(list.of.seurats)])

VlnPlot(object = Cells.merged,
        #log = FALSE,
        pt.size = 0,
        group.by = "orig.ident",
        #y.max = 2000,
        #ncol = 3,
        #cols = c("blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue"),
        features = c("nFeature_RNA", "nCount_RNA", "percent.mito")
)



current.list <- list.of.seurats #[c(1:(length(list.of.seurats)-1))]



### STANDARD INTEGRATION ###


setwd("/Users/ewhelan/Documents/scRNAseq/Rscripts/savefiles")

load("current.anchors3_cc_separate.Robj")
current.integrated <- IntegrateData(anchorset = current.anchors, dims = 1:30)


DefaultAssay(object = current.integrated) <- "integrated"
current.integrated <- ScaleData(object = current.integrated, verbose = FALSE)
current.integrated <- RunPCA(object = current.integrated,
                             npcs = 30,
                             verbose = TRUE)
current.integrated <- FindNeighbors(object = current.integrated)
current.integrated <- FindClusters(object = current.integrated,
                                   resolution = 1.1,
                                   #algorithm = 1,
                                   verbose = TRUE
)

current.integrated <- RunUMAP(object = current.integrated, reduction = "pca",
                              dims = 1:30)



HS.combined <- current.integrated
DefaultAssay(HS.combined) <- "RNA"


HS.combined <- CellCycleScoring(HS.combined, s.features = cc.genes$s.genes, g2m.features = cc.genes$g2m.genes, set.ident = TRUE)

p1 <- DimPlot(object = HS.combined, reduction = "umap", group.by = "orig.ident", raster=FALSE)

p2 <- FeaturePlot(HS.combined, reduction = "umap", features = c("SDC4", "KIT", "SYCP1", "PRM1"), order = T, raster=FALSE)

print(p1|p2)

DimPlot(HS.combined)
# save(HS.combined, file = "HS.combined.noCCregression.Robj")
load("HS.combined.noCCregression.Robj")

### re-integrate to remove cell cycle effects

#Try splitting and re-integrating. Uhh, see if this works. 
DefaultAssay(HS.combined) <- "RNA"
Idents(HS.combined) <- "Phase"
HS_S <- subset(HS.combined, idents = "S")
HS_G1 <- subset(HS.combined, idents = "G1")
HS_G2M <- subset(HS.combined, idents = "G2M")

n.dims <- 25

current.list <- c(HS_S, HS_G1, HS_G2M)
for (i in 1:length(x = current.list)) {
  current.list[[i]] <- NormalizeData(object = current.list[[i]], verbose = FALSE)
  current.list[[i]] <- FindVariableFeatures(object = current.list[[i]], nfeatures = 3000, verbose = FALSE)
}

current.anchors <- FindIntegrationAnchors(object.list = current.list, dims = 1:n.dims)
# save(current.anchors, file="current.anchors.Robj")
current.integrated <- IntegrateData(anchorset = current.anchors, dims = 1:n.dims)
# save(current.integrated, file="current.integrated.Robj")


DefaultAssay(object = current.integrated) <- "integrated"
current.integrated <- ScaleData(object = current.integrated, verbose = FALSE)
current.integrated <- RunPCA(object = current.integrated,
                             npcs = n.dims,
                             verbose = TRUE)
current.integrated <- FindNeighbors(object = current.integrated)
current.integrated <- FindClusters(object = current.integrated,
                                   resolution = 1.1,
                                   #algorithm = 1,
                                   verbose = TRUE
)
# current.integrated <- RunTSNE(object = current.integrated, dims = 1:20)
current.integrated <- RunUMAP(object = current.integrated, reduction = "pca",
                              dims = 1:n.dims)

DimPlot(object = current.integrated, reduction = "umap", group.by = "Phase", label = T)
DimPlot(object = current.integrated, reduction = "umap", group.by = "experiment")

HS_reclustered_CC_integrated <- current.integrated
save(HS_reclustered_CC_integrated, file="HS_reclustered_CC_integrated.Robj")
load("HS_reclustered_CC_integrated.Robj")


DefaultAssay(HS_reclustered_CC_integrated) <- "RNA"
DimPlot(HS_reclustered_CC_integrated, group.by = "Phase")
FeaturePlot(HS_reclustered_CC_integrated, order = T, features = germ.cells, min.cutoff = 0)



### gene lists ###

Keren.gene.list <- c("POU5F1", "NANOG", "TCL1A", "PDPN", "PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1",
                     "UTF1", "KIT", "ASB9", "NANOS2", "PIWIL4", "MORC1", "MAGEB2", "TDRD1", "RHOXF1", "MAGEC2", "ZBTB16", "EGR4", "RET", "MAGEA4",
                     "SOHLH1", "DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8", "DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9", "ZCWPW1", "RAD51AP2", "SYCP3",
                     "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1", "MLH3", "PIWIL1", "MLH1", "SPATA8", "CCDC112", "OVOL1", "OVOL2",
                     "CCDC42", "ACR", "CDRT15", "TEKT1", "H1FOO", "TJP3", "FAM24A", "CATSPER3", "SPACA3", "SPACA1", "ACRV1", "SAXO1", "CAPZA3", "HEMGN", "CA2", "PRM3",
                     "TNP2", "TSSK6", "SMCP", "SPEM2", "LELP1", "MKI67", "UBE2C", "TOP2A")

gene.list.PGC <- c("POU5F1", "NANOG", "TCL1A", "PDPN", "PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1",
                   "UTF1", "KIT", "ASB9", "NANOS2")
gene.list.spermatogonia <- c(                 
  "PIWIL4", "MORC1", "MAGEB2", "TDRD1", "RHOXF1", "MAGEC2", "ZBTB16", "EGR4", "RET")
gene.list.diff.spermatgonia <- c(
  
  "MAGEA4",
  "SOHLH1", "DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8")

gene.list.early.meiosis <- c(
  "DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9", "ZCWPW1", "RAD51AP2", "SYCP3",
  "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1", "MLH3")

gene.list.late.meiosis <- c(
  "PIWIL1", "MLH1", "SPATA8", "CCDC112", "OVOL1", "OVOL2",
  "CCDC42", "ACR", "CDRT15", "TEKT1", "H1FOO", "TJP3", "FAM24A", "CATSPER3", "SPACA3", "SPACA1", "ACRV1", "SAXO1")
gene.list.spermatids <- c("CAPZA3", "HEMGN", "CA2", "PRM3",
                          "TNP2", "TSSK6", "SMCP", "SPEM2", "LELP1", "MKI67", "UBE2C", "TOP2A")


gene.list.sertoli <- c("CLU", "AMHR2", "SOX9", "CTSL")
genes.macrophage <- c("APOE", "DAB2", "CD74", "ADGRE1")
genes.endothelial <- c("VWF", "TIE1", "TEK")
genes.myoid <- c("ACTA2", "MYH11", "MYL6", "PDGFRB")
genes.leydig <- c("CYP17A1", "CYP11A1", "STAR", "HSD3B1")
genes.T <- c("CD8A", "CD4", "CD3E", "CD28")
genes.telocytes <- c("DCN", "GSN", "TCF21")

germ.cells <- c("POU5F1", "NANOG", "ZBTB16", "SOHLH1", "PRDM9", "ACRV1", "SPACA1", "PRM1", "TNP2")
somatic.cells <- c("SOX9", "APOE", "VWF", "ACTA2", "STAR", "DCN", "TCF21", "ACTA2", "DAB2")

DimPlot(HS_reclustered_CC_integrated,label = T, group.by = "seurat_clusters")
DimPlot(HS_reclustered_CC_integrated,label = T, group.by = "stage")
FeaturePlot(HS_reclustered_CC_integrated, features = gene.list.spermatids, min.cutoff = 0, order = T)
FeaturePlot(HS_reclustered_CC_integrated, order = T, features = gene.list.spermatids, min.cutoff = 0)
DimPlot(HS_reclustered_CC_integrated, group.by = "Phase")



apoptotic.genes <- c("ADD1","AIFM3","ANKH","ANXA1","APP","ATF3","AVPR1A","BAX","BCAP31","BCL10","BCL2L1",
                     "BCL2L10","BCL2L11","BCL2L2","BGN","BID","BIK","BIRC3","BMF","BMP2","BNIP3L","BRCA1",
                     "BTG2","BTG3","CASP1","CASP2","CASP3","CASP4","CASP6","CASP7","CASP8","CASP9","CAV1",
                     "CCNA1","CCND1","CCND2","CD14","CD2","CD38","CD44","CD69","CDC25B","CDK2","CDKN1A","CDKN1B",
                     "CFLAR","CLU","CREBBP","CTH","CTNNB1","CYLD","DAP","DAP3","DCN","DDIT3","DFFA","DIABLO","DNAJA1","DNAJC3","DNM1L","DPYD","EBP",
                     "EGR3","EMP1","ENO2","ERBB2","ERBB3","EREG","ETF1","F2","F2R","FAS","FASLG","FDXR","FEZ1","GADD45A","GADD45B","GCH1","GNA15",#"GPX1",
                     "GPX3","GPX4","GSN","GSR","GSTM1","GUCY2D",#"H1-0",
                     "HGF","HMGB2","HMOX1","HSPB1","IER3","IFITM3","IFNB1","IFNGR1","IGF2R","IGFBP6",
                     "IL18","IL1A","IL1B","IL6","IRF1","ISG20","JUN","KRT18","LEF1","LGALS3","LMNA","PLPPR4","LUM","MADD","MCL1","MGMT","MMP2","NEDD9",
                     "NEFH","PAK1","PDCD4","PDGFRB","PEA15","PLAT","PLCB2","PMAIP1","PPP2R5B","PPP3R1","PPT1","PRF1","PSEN1","PSEN2","PTK2","RARA",
                     "RELA","RETSAT","RHOB","RHOT2","RNASEL","ROCK1","SAT1","SATB1","SC5D","SLC20A1","SMAD7","SOD1","SOD2","SPTAN1","SQSTM1","TAP1",
                     "TGFB2","TGFBR3","TIMP1","TIMP2","TIMP3","TNF","TNFRSF12A","TNFSF10","TOP2A","TSPO","TXNIP","VDAC2","WEE1","XIAP")


gene.list.PGC.E <- c("POU5F1", "NANOS3", "MIF", "ACTG2")
gene.list.PGC.L <- c("SOX15", "L1TD1", "PDPN", "TUBA1C")
gene.list.M <- c("KIT", "CCND2", "EEF1G", "CCPG1" )
gene.list.M2T1 <- c("MAP4K4", "ASB9", "PRXL2A", "TLE5")
gene.list.T1 <- c("DCAF4L1", "ID1", "ID3", "MORC1")
gene.list.spermatogonia <- c("UTF1", "PAGE2", "RHOXF1", "MAGEB2")
gene.list.diff.spermatgonia.E <- c("SOHLH1", "KCNQ2", "MAGEC2", "MAGEA4")
gene.list.diff.spermatgonia.L   <- c("DMRT1", "SOHLH2", "ESX1", "CTCFL")
gene.list.preleptotene <- c("RAD51AP2", "SCML1", "TEX19", "PRDM9")
gene.list.leptotene <- c("TEX101", "SYCP1", "DMRTC2", "TDRG1")
gene.list.zygotene <- c("CTAGE1", "EQTN", "SPACA1", "GOLGA6L2")
gene.list.late.meiosis <- c(
  "PIWIL1", "MLH1", "SPATA8", "CCDC112", "OVOL1", "OVOL2",
  "CCDC42", "ACR", "CDRT15", "TEKT1", "H1FOO", "TJP3", "FAM24A", "CATSPER3", "SPACA3", "SPACA1", "ACRV1", "SAXO1")
gene.list.spermatids <- c("CAPZA3", "HEMGN", "CA2", "PRM3",
                          "TNP2", "TSSK6", "SMCP", "SPEM2", "LELP1", "MKI67", "UBE2C", "TOP2A")



gene.list.hermann.SSC <- c("ID4", "NANOS2", "TCN2")
gene.list.hermann.progenitor <- c("EGR4", "PLPPR5", "TSPAN33")
gene.list.hermann.earlydiff <- c("ASB9", "NANOS3", "TMEM55A")
gene.list.hermann.latediff <- c("CDT1", "ISOC1", "NMT2")
gene.list.hermann.prelep <- c("CT55", "DHRS13", "HSD17B14")
gene.list.hermann.lepzyg <- c("HIST3H3", "MEIOB", "TEX101")
gene.list.hermann.pachytene <- c("C9orf57", "CETN3", "MGAT4D")
gene.list.hermann.dipsecondary <- c("ARMC1", "PHOSPHO2", "TPPP3")
gene.list.hermann.earlyround <- c("C17orf98", "ENPP2", "LRRC3B")
gene.list.hermann.midround <- c("PRSS37", "PRSS58", "TP53TG5")
gene.list.hermann.lateround <- c("C17orf74", "FSCN3", "PHOSPHO1")
FeaturePlot(HS_germ_cells, features = gene.list.hermann.lateround, min.cutoff = 0, order = T)

HS_reclustered_CC_integrated <- AddModuleScore(HS_reclustered_CC_integrated,
                                               features = list(apoptotic.genes),
                                               name="Apoptotic")
HS_reclustered_CC_integrated <- AddModuleScore(HS_reclustered_CC_integrated,
                                               features = list(gene.list.PGC.E),
                                               name="PGCE")
HS_reclustered_CC_integrated <- AddModuleScore(HS_reclustered_CC_integrated,
                                               features = list(gene.list.PGC.L),
                                               name="PGCL")
HS_reclustered_CC_integrated <- AddModuleScore(HS_reclustered_CC_integrated,
                                               features = list(gene.list.M),
                                               name="M")
HS_reclustered_CC_integrated <- AddModuleScore(HS_reclustered_CC_integrated,
                                               features = list(gene.list.M2T1),
                                               name="M2T")
HS_reclustered_CC_integrated <- AddModuleScore(HS_reclustered_CC_integrated,
                                               features = list(gene.list.T1),
                                               name="T")
HS_reclustered_CC_integrated <- AddModuleScore(HS_reclustered_CC_integrated,
                                               features = list(gene.list.spermatogonia),
                                               name="spg")
HS_reclustered_CC_integrated <- AddModuleScore(HS_reclustered_CC_integrated,
                                               features = list(gene.list.diff.spermatgonia.E),
                                               name="diffspgE")
HS_reclustered_CC_integrated <- AddModuleScore(HS_reclustered_CC_integrated,
                                               features = list(gene.list.diff.spermatgonia.L),
                                               name="diffspgL")
HS_reclustered_CC_integrated <- AddModuleScore(HS_reclustered_CC_integrated,
                                               features = list(gene.list.preleptotene),
                                               name="preleptotene")
HS_reclustered_CC_integrated <- AddModuleScore(HS_reclustered_CC_integrated,
                                               features = list(gene.list.leptotene),
                                               name="leptotene")
HS_reclustered_CC_integrated <- AddModuleScore(HS_reclustered_CC_integrated,
                                               features = list(gene.list.zygotene),
                                               name="zygotene")
HS_reclustered_CC_integrated <- AddModuleScore(HS_reclustered_CC_integrated,
                                               features = list(gene.list.late.meiosis),
                                               name="pachytene_2ndary")
HS_reclustered_CC_integrated <- AddModuleScore(HS_reclustered_CC_integrated,
                                               features = list(gene.list.spermatids),
                                               name="spermatids")



plot_grid(ncol = 5,
          
          FeaturePlot(HS_reclustered_CC_integrated,
                      features = "Apoptotic1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_reclustered_CC_integrated,
                      features = "PGCE1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_reclustered_CC_integrated,
                      features = "PGCL1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_reclustered_CC_integrated,
                      features = "M1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_reclustered_CC_integrated,
                      features = "M2T1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_reclustered_CC_integrated,
                      features = "T1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_reclustered_CC_integrated,
                      features = "spg1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_reclustered_CC_integrated,
                      features = "diffspgE1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_reclustered_CC_integrated,
                      features = "diffspgL1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_reclustered_CC_integrated,
                      features = "preleptotene1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_reclustered_CC_integrated,
                      features = "leptotene1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_reclustered_CC_integrated,
                      features = "zygotene1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_reclustered_CC_integrated,
                      features = "pachytene_2ndary1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_reclustered_CC_integrated,
                      features = "spermatids1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme()
          
)




HS_germ_cells <- subset(HS_reclustered_CC_integrated, idents = c(25, 44, 15, 30, 31, 36, 33, 12), invert = F) #subset 1

DefaultAssay(HS_germ_cells) <- "integrated"
HS_germ_cells <- RunPCA(object = HS_germ_cells,
                        npcs = 30,
                        verbose = TRUE)
HS_germ_cells <- FindNeighbors(object = HS_germ_cells)
HS_germ_cells <- FindClusters(object = HS_germ_cells,
                              resolution = 2,
                              #algorithm = 1,
                              verbose = TRUE
)

HS_germ_cells <- RunUMAP(object = HS_germ_cells, reduction = "pca",
                         dims = 1:30)

DimPlot(HS_germ_cells, group.by = "stage")
DimPlot(HS_germ_cells, label = T)
DefaultAssay(HS_germ_cells) <- "RNA"
FeaturePlot(HS_germ_cells, features = gene.list.spermatids, min.cutoff = 0, order = T)


HS_germ_cells <- AddModuleScore(HS_germ_cells,
                                features = list(apoptotic.genes),
                                name="Apoptotic")
HS_germ_cells <- AddModuleScore(HS_germ_cells,
                                features = list(gene.list.PGC.E),
                                name="PGCE")
HS_germ_cells <- AddModuleScore(HS_germ_cells,
                                features = list(gene.list.PGC.L),
                                name="PGCL")
HS_germ_cells <- AddModuleScore(HS_germ_cells,
                                features = list(gene.list.M),
                                name="M")
HS_germ_cells <- AddModuleScore(HS_germ_cells,
                                features = list(gene.list.M2T1),
                                name="M2T")
HS_germ_cells <- AddModuleScore(HS_germ_cells,
                                features = list(gene.list.T1),
                                name="T")
HS_germ_cells <- AddModuleScore(HS_germ_cells,
                                features = list(gene.list.spermatogonia),
                                name="spg")
HS_germ_cells <- AddModuleScore(HS_germ_cells,
                                features = list(gene.list.diff.spermatgonia.E),
                                name="diffspgE")
HS_germ_cells <- AddModuleScore(HS_germ_cells,
                                features = list(gene.list.diff.spermatgonia.L),
                                name="diffspgL")
HS_germ_cells <- AddModuleScore(HS_germ_cells,
                                features = list(gene.list.preleptotene),
                                name="preleptotene")
HS_germ_cells <- AddModuleScore(HS_germ_cells,
                                features = list(gene.list.leptotene),
                                name="leptotene")
HS_germ_cells <- AddModuleScore(HS_germ_cells,
                                features = list(gene.list.zygotene),
                                name="zygotene")
HS_germ_cells <- AddModuleScore(HS_germ_cells,
                                features = list(gene.list.late.meiosis),
                                name="late.meiosis")
HS_germ_cells <- AddModuleScore(HS_germ_cells,
                                features = list(gene.list.spermatids),
                                name="spermatids")



plot_grid(ncol = 5,
          
          FeaturePlot(HS_germ_cells,
                      features = "Apoptotic1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_germ_cells,
                      features = "PGCE1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_germ_cells,
                      features = "PGCL1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_germ_cells,
                      features = "M1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells,
                      features = "M2T1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells,
                      features = "T1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_germ_cells,
                      features = "spg1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_germ_cells,
                      features = "diffspgE1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_germ_cells,
                      features = "diffspgL1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_germ_cells,
                      features = "preleptotene1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_germ_cells,
                      features = "leptotene1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_germ_cells,
                      features = "zygotene1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_germ_cells,
                      features = "late.meiosis1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme(),
          
          FeaturePlot(HS_germ_cells,
                      features = "spermatids1", label = F, repel = FALSE) +
            scale_color_viridis(option="plasma") + DarkTheme()
          
)

HS_germ_cells <- AddModuleScore(HS_germ_cells, features = list(gene.list.hermann.SSC), name="hermann.sscs")
HS_germ_cells <- AddModuleScore(HS_germ_cells, features = list(gene.list.hermann.progenitor), name="hermann.progenitor")
HS_germ_cells <- AddModuleScore(HS_germ_cells, features = list(gene.list.hermann.earlydiff), name="hermann.earlydiff")
HS_germ_cells <- AddModuleScore(HS_germ_cells, features = list(gene.list.hermann.latediff), name="hermann.latediff")
HS_germ_cells <- AddModuleScore(HS_germ_cells, features = list(gene.list.hermann.prelep), name="hermann.prelep")
HS_germ_cells <- AddModuleScore(HS_germ_cells, features = list(gene.list.hermann.lepzyg), name="hermann.lep_zyg")
HS_germ_cells <- AddModuleScore(HS_germ_cells, features = list(gene.list.hermann.pachytene), name="hermann.pachytene")
HS_germ_cells <- AddModuleScore(HS_germ_cells, features = list(gene.list.hermann.dipsecondary), name="hermann.diplotene_secondary")
HS_germ_cells <- AddModuleScore(HS_germ_cells, features = list(gene.list.hermann.earlyround), name="hermann.earlyround")
HS_germ_cells <- AddModuleScore(HS_germ_cells, features = list(gene.list.hermann.midround), name="hermann.midround")
HS_germ_cells <- AddModuleScore(HS_germ_cells, features = list(gene.list.hermann.lateround), name="hermann.lateround")

plot_grid(ncol = 5,
          FeaturePlot(HS_germ_cells,features = "hermann.sscs1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells, features = "hermann.progenitor1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells, features = "hermann.earlydiff1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells, features = "hermann.latediff1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells, features = "hermann.prelep1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells, features = "hermann.lep_zyg1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells, features = "hermann.pachytene1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells, features = "hermann.diplotene_secondary1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells, features = "hermann.earlyround1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells, features = "hermann.midround1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells, features = "hermann.lateround1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme()
)


list_of_cols <- rep("gray", each = 33)
list_of_cols[21+1] <- "blue"
DimPlot(HS_germ_cells, cols = list_of_cols, label = T)

subset(HS_germ_cells, idents = c(21, 0, 1, 9, 6, 16, 18), invert = TRUE) -> HS_germ_cells2

DefaultAssay(HS_germ_cells2) <- "integrated"
HS_germ_cells2 <- ScaleData(object = HS_germ_cells2, verbose = FALSE)
HS_germ_cells2 <- RunPCA(object = HS_germ_cells2,
                         npcs = 20,
                         verbose = TRUE)
HS_germ_cells2 <- FindNeighbors(object = HS_germ_cells2)
HS_germ_cells2 <- FindClusters(object = HS_germ_cells2,
                               resolution = 2.5,
                               #algorithm = 1,
                               verbose = TRUE
)
# current.integrated <- RunTSNE(object = current.integrated, dims = 1:20)
HS_germ_cells2 <- RunUMAP(object = HS_germ_cells2, reduction = "pca",
                          dims = 1:20)

DimPlot(HS_germ_cells2, label = T)


subset(HS_germ_cells2, idents = c(18, 21, 29, 31), invert = TRUE) -> HS_germ_cells_final
DefaultAssay(HS_germ_cells_final) <- "integrated"
HS_germ_cells_final <- ScaleData(object = HS_germ_cells_final, verbose = FALSE)
HS_germ_cells_final <- RunPCA(object = HS_germ_cells_final,
                              npcs = 16,
                              verbose = TRUE)
HS_germ_cells_final <- FindNeighbors(object = HS_germ_cells_final)
HS_germ_cells_final <- FindClusters(object = HS_germ_cells_final,
                                    resolution = 1.5,
                                    #algorithm = 1,
                                    verbose = TRUE
)
# current.integrated <- RunTSNE(object = current.integrated, dims = 1:20)
HS_germ_cells_final <- RunUMAP(object = HS_germ_cells_final, reduction = "pca",
                               dims = 1:16)

DimPlot(HS_germ_cells_final, label = T)

DefaultAssay(HS_germ_cells_final) <- "RNA"


plot_grid(ncol = 5,
          FeaturePlot(HS_germ_cells_final,features = "hermann.sscs1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final, features = "hermann.progenitor1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final, features = "hermann.earlydiff1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final, features = "hermann.latediff1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final, features = "hermann.prelep1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final, features = "hermann.lep_zyg1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final, features = "hermann.pachytene1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final, features = "hermann.diplotene_secondary1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final, features = "hermann.earlyround1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final, features = "hermann.midround1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final, features = "hermann.lateround1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme()
)


plot_grid(ncol = 6,
          DimPlot(HS_germ_cells_final, group.by = "monocle_clusters", label = T) + NoLegend(),
          FeaturePlot(HS_germ_cells_final,features = "PGCE1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final,features = "PGCL1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final,features = "M1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final,features = "M2T1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final,features = "T1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          # FeaturePlot(HS_germ_cells_final,features = "hermann.sscs1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final,features = "spg1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final,features = "hermann.progenitor1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final,features = "hermann.earlydiff1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          # FeaturePlot(HS_germ_cells_final,features = "diffspgE1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final,features = "hermann.latediff1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          # FeaturePlot(HS_germ_cells_final,features = "diffspgL1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final,features = "preleptotene1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final, features = "hermann.prelep1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final,features = "leptotene1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final,features = "zygotene1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final, features = "hermann.lep_zyg1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final, features = "hermann.pachytene1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final, features = "hermann.diplotene_secondary1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          # FeaturePlot(HS_germ_cells_final, features = "hermann.earlyround1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          # FeaturePlot(HS_germ_cells_final, features = "hermann.midround1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          # FeaturePlot(HS_germ_cells_final, features = "hermann.lateround1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_germ_cells_final, features = "spermatids1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme()
)




### monocle pseudotime ###

# cds <- as.cell_data_set(NCG_germ4subset)
cds <- as.cell_data_set(HS_germ_cells_final)


cds <- cluster_cells(cds, reduction_method = c("UMAP"), #cluster_method = "louvain",
                     resolution = 1.5e-3, 
                     k = 8,
                     partition_qval = 0.05
)

plot_cells(cds, show_trajectory_graph = FALSE)
cds <- learn_graph(cds, 
                   close_loop = F,
                   # learn_graph_control=list(ncenter=500, nn.method = "nn2", minimal_branch_len=10), #250, 18 is pretty good
                   use_partition = F)
print(plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)) #+ ggtitle(a))

# Figure 2C #######################################
cds <- order_cells(cds)#, root_cells = max.avp)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)+ NoLegend()

cds.pseudotime <- pseudotime(cds, reduction_method = "UMAP")
cds.clusters <- clusters(cds, reduction_method = "UMAP")



#now add the pseudotime to the original object

HS_germ_cells_final$pseudotime <- cds.pseudotime
HS_germ_cells_final$monocle_clusters <- cds.clusters




### label cells ###

Idents(HS_germ_cells_final) <- HS_germ_cells_final[["monocle_clusters"]]
DimPlot(HS_germ_cells_final, group.by = "monocle_clusters", label = T)

new.cluster.ids <- character(length(levels(Idents(HS_germ_cells_final))))

PGCE.ids <- c(17, 4, 21, 14)
PGCL.ids <- c(12, 2)
M1.ids <- c(22, 13, 10)
M2T1.ids <- c(6, 3)
T1.ids <-c(5, 37)
spg.ids <- c(15, 25)
diffE.ids <- c(7)
diffL.ids <-c(28)
prelep.ids <- c(23)
lepzyg.ids <- c(27, 20, 24)
pachy.ids <-c(8, 32)
dipsec.ids <-c(16, 34, 33, 19)
spermatid.ids <- c(1, 9,35, 30, 18, 29, 36, 11, 31, 26)

new.cluster.ids[PGCE.ids] <- "PGCs (early)"
new.cluster.ids[PGCL.ids] <- "PGCs (late)"
new.cluster.ids[M1.ids] <- "M prospermatogonia"
new.cluster.ids[M2T1.ids] <- "M/T1 prospermatogonia"
new.cluster.ids[T1.ids] <- "T1 prospermatogonia"
new.cluster.ids[spg.ids] <- "Spermatogonia"
new.cluster.ids[diffE.ids] <- "Diff. spermatogonia (early)"
new.cluster.ids[diffL.ids] <- "Diff. spermatogonia (late)"
new.cluster.ids[prelep.ids] <- "Preleptotene spermatocytes"
new.cluster.ids[lepzyg.ids] <- "Leptotene/Zygotene"
new.cluster.ids[pachy.ids] <- "Pachytene spermatocytes"
new.cluster.ids[dipsec.ids] <- "Diplotene/2ndary spermatocytes"
new.cluster.ids[spermatid.ids] <- "Spermatids"



new.cluster.ids

names(x = new.cluster.ids) <- levels(x = HS_germ_cells_final)
HS_germ_cells_final <- RenameIdents(object = HS_germ_cells_final, new.cluster.ids)
HS_germ_cells_final$cell.type <- Idents(HS_germ_cells_final)

HS_germ_cells_final$cell.type <- factor(x = HS_germ_cells_final$cell.type, levels = c("PGCs (early)", 
                                                                                      "PGCs (late)",
                                                                                      "M prospermatogonia",
                                                                                      "M/T1 prospermatogonia",
                                                                                      "T1 prospermatogonia",
                                                                                      "Spermatogonia",
                                                                                      "Diff. spermatogonia (early)",
                                                                                      "Diff. spermatogonia (late)",
                                                                                      "Preleptotene spermatocytes",
                                                                                      "Leptotene/Zygotene",
                                                                                      "Pachytene spermatocytes",
                                                                                      "Diplotene/2ndary spermatocytes",
                                                                                      "Spermatids"))


Idents(HS_germ_cells_final) <- HS_germ_cells_final$cell.type
DimPlot(HS_germ_cells_final, group.by = "cell.type", label = T)
# save(HS_germ_cells_final, file="HS_germ_cells_final_2023.07.05.Robj")
# save(HS_germ_cells_final, file="HS_germ_cells_final_2023.07.06.Robj")
# save(HS_germ_cells_final, file="HS_germ_cells_final_2023.07.19.Robj")
load("HS_germ_cells_final_2023.07.19.Robj")

#Plot out designations with a feature of interest.
p1 <- DimPlot(HS_germ_cells_final, label = F, reduction = "umap", group.by = "cell.type")
p2 <- DimPlot(NCG_germ4, label = F, reduction = "umap", group.by = "cell.type")

FeaturePlot(HS_germ_cells_final, features = c("GAGE12H", "NANOS2", "XIST", "HELLPAR"))
DimPlot(HS_germ_cells_final, group.by = "cell.type", label = T)
DimPlot(NCG_germ4, group.by = "cell.type")
plot_grid(p1, p2)




### DEGs for Ingenuity ###

DimPlot(HS_germ_cells_final)

HS_germ_cells_final_subset <- subset(HS_germ_cells_final, idents = c("Leptotene/Zygotene", "Pachytene spermatocytes", "Diplotene/2ndary spermatocytes", "Spermatids"), invert = T)


# all.markers.subset <- filter(all.markers_NCG, avg_log2FC > 1)

PGCs.E.markers.all <- FindMarkers(HS_germ_cells_final_subset, ident.1 = "PGCs (early)")
PGCs.L.markers.all <- FindMarkers(HS_germ_cells_final_subset, ident.1 = "PGCs (late)")
M.ProSpermatogonia.markers.all <- FindMarkers(HS_germ_cells_final_subset, ident.1 = "M prospermatogonia")
M2T1.markers.all <- FindMarkers(HS_germ_cells_final_subset, ident.1 = "M/T1 prospermatogonia")
T.ProSpermatogonia.markers.all <- FindMarkers(HS_germ_cells_final_subset, ident.1 = "T1 prospermatogonia")
Spg.markers.all <- FindMarkers(HS_germ_cells_final_subset, ident.1 = "Spermatogonia")
Diff.Spg.E.markers.all <- FindMarkers(HS_germ_cells_final_subset, ident.1 = "Diff. spermatogonia (early)")
Diff.Spg.L.markers.all <- FindMarkers(HS_germ_cells_final_subset, ident.1 = "Diff. spermatogonia (late)")
Pre.Leptotene.markers.all <- FindMarkers(HS_germ_cells_final_subset, ident.1 = "Preleptotene spermatocytes")

setwd("/Users/ewhelan/Desktop/NCG/NCG_marker_genes")
write.csv(PGCs.E.markers.all, "HS_PGCs.E.markers.all.csv")
write.csv(PGCs.L.markers.all, "HS_PGCs.L.markers.all.csv")
write.csv(M.ProSpermatogonia.markers.all, "HS_M.ProSpermatogonia.markers.all.csv")
write.csv(M2T1.markers.all, "HS_M2T1.markers.all.csv")
write.csv(T.ProSpermatogonia.markers.all, "HS_T.ProSpermatogonia.markers.all.csv")
write.csv(Spg.markers.all, "HS_Spg.markers.all.csv")
write.csv(Diff.Spg.E.markers.all, "HS_Diff.Spg.E.markers.all.csv")
write.csv(Diff.Spg.L.markers.all, "HS_Diff.Spg.L.markers.all.csv")
write.csv(Pre.Leptotene.markers.all, "HS_Pre.Leptotene.markers.all.csv")


### ### ### ### ### ### ###

#ITEGRATE IN VIVO AND IN VITRO

### ### ### ### ### ### ###

FeaturePlot(HS_germ_cells_final, features = c("TNP1", "TNP2", "PRM1", "PRM2"))


DefaultAssay(HS_germ_cells_final) <- "RNA"
DefaultAssay(NCG_germ4) <- "RNA"
DefaultAssay(NCG_germ_macaque4) <- "RNA"

DimPlot(HS_germ_cells_final, group.by = "monocle_clusters")
DimPlot(NCG_germ4, group.by = "seurat_clusters")
HS_germ_cells_final$original_clusters_HS <- HS_germ_cells_final$monocle_clusters
HS_germ_cells_final$original_clusters <- HS_germ_cells_final$monocle_clusters

NCG_germ4$original_clusters_HS <- NCG_germ4$seurat_clusters
NCG_germ4$original_clusters <- NCG_germ4$seurat_clusters

HS_germ_cells_final$origin <- "in.vivo"
NCG_germ4$origin <- "in.vitro"

NCG_germ_macaque4$origin <- "in.vitro"
NCG_germ_macaque4$original_clusters <- NCG_germ_macaque4$seurat_clusters

#split by phase before integration to remove cell cycle effects

Idents(HS_germ_cells_final) <- "Phase"
HS_germ_cells_final_S <- subset(HS_germ_cells_final, idents = "S")
HS_germ_cells_final_G1 <- subset(HS_germ_cells_final, idents = "G1")
HS_germ_cells_final_G2M <- subset(HS_germ_cells_final, idents = "G2M")

Idents(NCG_germ4) <- "Phase"
NCG_germ4_S <- subset(NCG_germ4, idents = "S")
NCG_germ4_G1 <- subset(NCG_germ4, idents = "G1")
NCG_germ4_G2M <- subset(NCG_germ4, idents = "G2M")


current.list <- c(HS_germ_cells_final_S, HS_germ_cells_final_G1, HS_germ_cells_final_G2M, 
                  NCG_germ4_S, NCG_germ4_G1, NCG_germ4_G2M)


n.dims <- 30
for (i in 1:length(x = current.list)) {
  current.list[[i]] <- NormalizeData(object = current.list[[i]], verbose = FALSE)
  current.list[[i]] <- FindVariableFeatures(object = current.list[[i]], nfeatures = 3000, verbose = FALSE)
}

current.anchors <- FindIntegrationAnchors(object.list = current.list, dims = 1:n.dims)
# save(current.anchors, file="current.anchors.Robj")
current.integrated <- IntegrateData(anchorset = current.anchors, dims = 1:n.dims)
# save(current.integrated, file="current.integrated.Robj")

DefaultAssay(object = current.integrated) <- "integrated"
current.integrated <- ScaleData(object = current.integrated, vars.to.regress = c("nCount_RNA", "percent.mito"), verbose = FALSE)
current.integrated <- RunPCA(object = current.integrated,
                             npcs = n.dims,
                             verbose = TRUE)
current.integrated <- FindNeighbors(object = current.integrated)
current.integrated <- FindClusters(object = current.integrated,
                                   resolution = 0.8,
                                   #algorithm = 1,
                                   verbose = TRUE
)
# current.integrated <- RunTSNE(object = current.integrated, dims = 1:20)
current.integrated <- RunUMAP(object = current.integrated, reduction = "pca",
                              dims = 1:20)
DimPlot(current.integrated, split.by = "species")

for(loop.number in 1:10){
  
  current.integrated <- RunUMAP(object = current.integrated, reduction = "pca", n.neighbors = 20, 
                                n.components = 3,
                                dims = 1:25)
  p1 <- (DimPlot(current.integrated, label = T,
                 group.by = "species")+ggtitle(paste("n.neighbors = ", (5*loop.number))))
  p2 <- FeaturePlot(current.integrated, #gene.list
                    features = c("PRM1"),#, "SYCP1", "SYCP3", "MEIOC", "PRDM9", "DMC1", "SPO11"), 
                    order = T, 
                    # split.by = "species2",
                    min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")
  print(plot_grid(p1, p2))
}



###Try again with Seurat v5's integration system

setwd("~/Documents/scRNAseq/Rscripts/savefiles")
load("NCG_germ_macaque4_2024.01.03.Robj")
load("HS_germ_cells_final_2023.07.19.Robj") 
load("NCG_germ4_2023.07.23.Robj") #this is the final in vitro object
load("HS_NCG_2023.07.19.Robj") #this is the combined in vivo and in vitro object

DimPlot(NCG_germ_macaque4, group.by = "species")
DimPlot(HS_germ_cells_final, group.by = "species")
DimPlot(NCG_germ4, group.by = "species")

HS_NCG[["species"]] <- "human"
HS_germ_cells_final[["species"]] <- "human"

Idents(HS_NCG) <- "experiment"
HS_NCG <- RenameIdents(object = HS_NCG, `human in vivo` = "human_in_vivo")
HS_NCG <- RenameIdents(object = HS_NCG, `transplant` = "transplant_human")
HS_NCG[["experiment"]] <- Idents(HS_NCG)
DimPlot(HS_NCG)


HS_NCG[["source"]] <- HS_NCG[["experiment"]]
Idents(HS_NCG) <- "source"
table(Idents(HS_NCG))
HS_NCG <- RenameIdents(object = HS_NCG, `hPGCLC` = "transplant_human")
HS_NCG[["source"]] <- Idents(HS_NCG)
DimPlot(HS_NCG)


NCG_germ_macaque4[["experiment"]] <- "transplant_macaque"

NCG_germ_macaque4[["source"]] <- "macaque_xrTestis"
HS_germ_cells_final[["source"]] <- "human_inVivo"
NCG_germ4[["source"]] <- "human_xrTestis"

DimPlot(NCG_germ_macaque4, group.by = "Phase")
DimPlot(HS_germ_cells_final, group.by = "Phase")
DimPlot(NCG_germ4, group.by = "Phase")



### ### ### 
# Map macaque onto human
# # # # # #


#now try and map macaque onto human


#first, re-do the human


HS_germ_cells_final2 <- HS_germ_cells_final
RunUMAP(HS_germ_cells_final2, dims = 1:16, return.model = T)
DimPlot(HS_germ_cells_final2)

HS_germ_cells_final2[["RNA"]] <- split(HS_germ_cells_final2[["RNA"]], f = HS_germ_cells_final2$orig.ident)


# HS_germ_cells_final2[["RNA"]] <- split(HS_germ_cells_final2[["RNA"]], f = HS_germ_cells_final2$Phase)
HS_germ_cells_final2 <- NormalizeData(HS_germ_cells_final2)
HS_germ_cells_final2 <- FindVariableFeatures(HS_germ_cells_final2)
HS_germ_cells_final2 <- ScaleData(HS_germ_cells_final2)
HS_germ_cells_final2 <- RunPCA(HS_germ_cells_final2)
# HS_germ_cells_final2 <- FindNeighbors(HS_germ_cells_final2, dims = 1:30, reduction = "pca")
# HS_germ_cells_final2 <- FindClusters(HS_germ_cells_final2, resolution = 2, cluster.name = "unintegrated_clusters")
# HS_germ_cells_final2 <- RunUMAP(HS_germ_cells_final2, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
# DimPlot(HS_germ_cells_final2, reduction = "umap.unintegrated", group.by = c("experiment", "Phase"))
# FeaturePlot(HS_germ_cells_final2, features = "TFAP2C", reduction = "umap.unintegrated")

HS_germ_cells_final2 <- IntegrateLayers(object = HS_germ_cells_final2, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                        verbose = FALSE)

# re-join layers after integration
HS_germ_cells_final2[["RNA"]] <- JoinLayers(HS_germ_cells_final2[["RNA"]])

HS_germ_cells_final2 <- FindNeighbors(HS_germ_cells_final2, reduction = "integrated.cca", dims = 1:30)
HS_germ_cells_final2 <- FindClusters(HS_germ_cells_final2, resolution = 1)
for(loop.number in 15:30){
  HS_germ_cells_final2 <- RunUMAP(HS_germ_cells_final2, dims = 1:loop.number, reduction = "integrated.cca")
  
  print(DimPlot(HS_germ_cells_final2, reduction = "umap", group.by = c("seurat_clusters", "Phase"))+ggtitle(paste0("loop.number: ", loop.number)))
}

HS_germ_cells_final2 <- RunUMAP(HS_germ_cells_final2, dims = 1:25, reduction = "integrated.cca", return.model = TRUE)

FeaturePlot(HS_germ_cells_final2, reduction = "umap", features = c("TFAP2C", "DDX4", "PIWIL4", "MEIOB", "ACRV1", "PRM1"), ncol = 3,order = T) & DarkTheme() & scale_color_viridis(option = "C")


human.query <- NCG_germ4
human.reference <- HS_germ_cells_final2
DefaultAssay(human.query) <- "RNA"
human.query[["xrTestis.cell.types"]] <- human.query[["cell.type"]] 
human.reference[["inVivo.cell.types"]] <- human.reference[["cell.type"]] 
DefaultAssay(human.reference) <- "RNA"
DimPlot(human.reference, group.by = "inVivo.cell.types")

# human.reference <- NormalizeData(human.reference)
# human.reference <- FindVariableFeatures(human.reference)
# human.reference <- ScaleData(human.reference)
# human.reference <- RunPCA(human.reference)
# human.reference <- FindNeighbors(human.reference, dims = 1:30)
# human.reference <- FindClusters(human.reference)
# human.reference <- RunUMAP(human.reference, dims = 1:20, return.model=TRUE)
DimPlot(human.reference, group.by = c("experiment", "species"))


human.query <- NormalizeData(human.query)
pancreas.anchors <- FindTransferAnchors(reference = human.reference, query = human.query, dims = 1:30,
                                        reference.reduction = "pca")
predictions <- TransferData(anchorset = pancreas.anchors, refdata = human.reference$cell.type, dims = 1:30)
human.query <- AddMetaData(human.query, metadata = predictions)


human.query$prediction.match <- human.query$predicted.id == human.query$cell.type
table(human.query$prediction.match)
table(human.query$predicted.id)


human.reference <- RunUMAP(human.reference, dims = 1:30, reduction = "integrated.cca", return.model = TRUE)
human.query <- MapQuery(anchorset = pancreas.anchors, reference = human.reference, query = human.query,
                        refdata = list(celltype = "cell.type"), reference.reduction = "pca", reduction.model = "umap")



human.query <- TransferData(anchorset = pancreas.anchors, reference = human.reference, query = human.query,
                            refdata = list(celltype = "cell.type"))
human.query <- IntegrateEmbeddings(anchorset = pancreas.anchors, reference = human.reference, query = human.query,
                                   new.reduction.name = "ref.pca")
human.query <- ProjectUMAP(query = human.query, query.reduction = "ref.pca", reference = human.reference,
                           reference.reduction = "pca", reduction.model = "umap")

# p1 <- DimPlot(human.reference, reduction = "umap", group.by = "cell.type", label = TRUE, label.size = 3,
#               repel = TRUE) + NoLegend() + ggtitle("Reference annotations") + xlim(c(-15,15)) + ylim(c(-15,15))
# p2 <- DimPlot(human.query, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
#               label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Human query transferred labels")+ xlim(c(-15,15)) + ylim(c(-15,15))
# p3 <- DimPlot(human.query, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
#               label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Human query original labels")+ xlim(c(-15,15)) + ylim(c(-15,15))
# p1 + p2


#now repeat for macaque


macaque.query <- NCG_germ_macaque4
human.reference <- HS_germ_cells_final2
DefaultAssay(macaque.query) <- "RNA"
DefaultAssay(human.reference) <- "RNA"
macaque.query[["xrTestis.cell.types"]] <- macaque.query[["cell.type"]] 

# human.reference <- NormalizeData(human.reference)
# human.reference <- FindVariableFeatures(human.reference)
# human.reference <- ScaleData(human.reference)
# human.reference <- RunPCA(human.reference)
# human.reference <- FindNeighbors(human.reference, dims = 1:30)
# human.reference <- FindClusters(human.reference)
# human.reference <- RunUMAP(human.reference, dims = 1:20, return.model=TRUE)
DimPlot(human.reference, group.by = c("experiment", "species"))


macaque.query <- NormalizeData(macaque.query)
pancreas.anchors <- FindTransferAnchors(reference = human.reference, query = macaque.query, dims = 1:30,
                                        reference.reduction = "pca")
predictions <- TransferData(anchorset = pancreas.anchors, refdata = human.reference$cell.type, dims = 1:30)
macaque.query <- AddMetaData(macaque.query, metadata = predictions)


macaque.query$prediction.match <- macaque.query$predicted.id == macaque.query$cell.type
table(macaque.query$prediction.match)
table(macaque.query$predicted.id)


human.reference <- RunUMAP(human.reference, dims = 1:30, n.components = 3L, reduction = "integrated.cca", return.model = TRUE)
macaque.query <- MapQuery(anchorset = pancreas.anchors, reference = human.reference, query = macaque.query,
                          refdata = list(celltype = "cell.type"), reference.reduction = "pca", reduction.model = "umap")



macaque.query <- TransferData(anchorset = pancreas.anchors, reference = human.reference, query = macaque.query,
                              refdata = list(celltype = "cell.type"))
macaque.query <- IntegrateEmbeddings(anchorset = pancreas.anchors, reference = human.reference, query = macaque.query,
                                     new.reduction.name = "ref.pca")
macaque.query <- ProjectUMAP(query = macaque.query, query.reduction = "ref.pca", reference = human.reference,
                             reference.reduction = "pca", reduction.model = "umap")


p1 <- DimPlot(human.reference, reduction = "umap", group.by = "inVivo.cell.types", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Reference annotations") + xlim(c(-15,15)) + ylim(c(-15,15))
p2 <- DimPlot(human.query, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Human query transferred labels")+ xlim(c(-15,15)) + ylim(c(-15,15))
p3 <- DimPlot(human.query, reduction = "ref.umap", group.by = "xrTestis.cell.types", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Human query original labels")+ xlim(c(-15,15)) + ylim(c(-15,15))

# p3 <- DimPlot(human.reference, reduction = "umap", group.by = "cell.type", label = TRUE, label.size = 3,
#               repel = TRUE) + NoLegend() + ggtitle("Reference annotations")+ xlim(c(-15,15)) + ylim(c(-15,15))
p4 <- DimPlot(macaque.query, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Macaque query transferred labels")+ xlim(c(-15,15)) + ylim(c(-15,15))
p1 + p2 + p3 + p4




DimPlot(HS_NCG)
macaque.query <- NCG_germ_macaque4
human.reference <- HS_NCG
DefaultAssay(macaque.query) <- "RNA"
DefaultAssay(human.reference) <- "RNA"

human.reference <- RunUMAP(object = human.reference, reduction = "pca", dims = 1:31, return.model = T)
Idents(human.reference) <- "experiment"
human.reference <- subset(human.reference, idents = "human in vivo") #c("hPGCLC", "transplant"))
DimPlot(human.reference)


# # run standard anlaysis workflow
human.reference <- NormalizeData(human.reference)
human.reference <- FindVariableFeatures(human.reference)

human.reference <- ScaleData(human.reference)
human.reference <- RunPCA(human.reference, n.pcs = 30)
human.reference <- FindNeighbors(human.reference, dims = 1:30, reduction = "pca")
human.reference <- FindClusters(human.reference, resolution = 2, cluster.name = "unintegrated_clusters")

setwd("~/Desktop/human.output")

for(loop.number in 5:30){
  human.reference <- RunUMAP(human.reference, dims = 1:loop.number, reduction = "pca", reduction.name = "umap.unintegrated")
  
  plot <- DimPlot(human.reference, reduction = "umap.unintegrated") + ggtitle(loop.number)
  
  # Generate a dynamic file name based on the gene name
  file_name <- paste0("dims", loop.number, "_human.clusters.png")  # Adjust the file extension as needed
  
  # Save the plot to a file using the dynamic file name
  ggsave(file = file_name, plot = plot, width = 5, height = 5)  # Replace X and Y with actual dimensions
  
  plot <- FeaturePlot(human.reference, features = c("TFAP2C", "DDX4", "MEIOB", "ACRV1")) & DarkTheme() & scale_color_viridis(option = "C")
  
  # Generate a dynamic file name based on the gene name
  file_name <- paste0("dims", loop.number, "_human.genes.png")  # Adjust the file extension as needed
  
  # Save the plot to a file using the dynamic file name
  ggsave(file = file_name, plot = plot, width = 5, height = 5)  # Replace X and Y with actual dimensions
  
}

FeaturePlot(human.reference, features = c("TFAP2C", "DDX4"))

human.reference <- RunUMAP(human.reference, dims = 1:22, reduction = "pca", reduction.name = "umap.unintegrated")


#now do the pseudotime for this object:




human.reference[["UMAP"]] <- human.reference[["umap.unintegrated"]]

cds <- as.cell_data_set(human.reference, reduction = "UMAP")


cds <- cluster_cells(cds, reduction_method = c("UMAP"), #cluster_method = "louvain",
                     resolution = 1.5e-4, 
                     k = 8,
                     partition_qval = 0.05
)

plot_cells(cds, show_trajectory_graph = FALSE)
cds <- learn_graph(cds, 
                   close_loop = F,
                   # learn_graph_control=list(ncenter=500, nn.method = "nn2", minimal_branch_len=10), #250, 18 is pretty good
                   use_partition = F)
print(plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)) #+ ggtitle(a))

# Figure 2C #######################################
cds <- order_cells(cds)#, root_cells = max.avp)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)+ NoLegend()

cds.pseudotime <- pseudotime(cds, reduction_method = "UMAP")
cds.clusters <- clusters(cds, reduction_method = "UMAP")



#now add the pseudotime to the original object

human.reference$pseudotime <- cds.pseudotime
human.reference$monocle_clusters <- cds.clusters



macaque.query <- NormalizeData(macaque.query)
pancreas.anchors <- FindTransferAnchors(reference = human.reference, query = macaque.query, dims = 1:22,
                                        reference.reduction = "pca")
predictions <- TransferData(anchorset = pancreas.anchors, refdata = human.reference$cell.type, dims = 1:22)
macaque.query <- AddMetaData(macaque.query, metadata = predictions)


macaque.query$prediction.match <- macaque.query$predicted.id == macaque.query$cell.type
table(macaque.query$prediction.match)
table(macaque.query$predicted.id)


# human.reference <- RunUMAP(human.reference, dims = 1:30, n.components = 3L, reduction = "integrated.cca", return.model = TRUE)
macaque.query <- MapQuery(anchorset = pancreas.anchors, reference = human.reference, query = macaque.query,
                          refdata = list(celltype = "cell.type"), reference.reduction = "pca", reduction.model = "umap")



macaque.query <- TransferData(anchorset = pancreas.anchors, reference = human.reference, query = macaque.query,
                              refdata = list(celltype = "cell.type"))
macaque.query <- IntegrateEmbeddings(anchorset = pancreas.anchors, reference = human.reference, query = macaque.query,
                                     new.reduction.name = "ref.pca")
macaque.query <- ProjectUMAP(query = macaque.query, query.reduction = "ref.pca", reference = human.reference,
                             reference.reduction = "pca", reduction.model = "umap")


# p1 <- DimPlot(human.reference, reduction = "umap", group.by = "inVivo.cell.types", label = TRUE, label.size = 3,
#               repel = TRUE) + NoLegend() + ggtitle("Reference annotations") + xlim(c(-15,15)) + ylim(c(-15,15))
# p2 <- DimPlot(human.query, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
#               label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Human query transferred labels")+ xlim(c(-15,15)) + ylim(c(-15,15))
# p3 <- DimPlot(human.query, reduction = "ref.umap", group.by = "xrTestis.cell.types", label = TRUE,
#               label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Human query original labels")+ xlim(c(-15,15)) + ylim(c(-15,15))
# 
# # p3 <- DimPlot(human.reference, reduction = "umap", group.by = "cell.type", label = TRUE, label.size = 3,
# #               repel = TRUE) + NoLegend() + ggtitle("Reference annotations")+ xlim(c(-15,15)) + ylim(c(-15,15))
# p4 <- DimPlot(macaque.query, reduction = "ref.umap", group.by = "predicted.id", label = TRUE,
#               label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Macaque query transferred labels")+ xlim(c(-15,15)) + ylim(c(-15,15))
# p1 + p2 + p3 + p4


p1 <- DimPlot(human.reference) + xlim(c(-15,15)) + ylim(c(-15,15)) + NoLegend()
p2 <- DimPlot(macaque.query, group.by = "original_clusters", label = T) + xlim(c(-15,15)) + ylim(c(-15,15)) + NoLegend()

p1 + p2


FeaturePlot(human.reference, features = c("DDX4", "MEIOB", "ACRV1", "PRM1"))

merged.seurat <- merge(human.reference, macaque.query, merge.dr = T)
DimPlot(merged.seurat)





#orange


eoin.palette.light <- c("antiquewhite3",
                        "lightcyan3", "pink3", "brown",  "sienna1", "aquamarine3","dodgerblue1", "yellowgreen", "lightgoldenrod2","darkorchid1", "darkgreen", "darkblue")

eoin.palette.dark <- c("gray28",
                       "cyan4", "pink3", "brown",  "sienna3", "aquamarine3","dodgerblue3", "darkgreen", "mediumpurple3","darkorchid1", "yellowgreen", "lightgoldenrod2", "darkblue")

eoin.palette.dark2 <- c("antiquewhite3","gray28",
                        "cyan4", "pink3", "brown",  "sienna3", "aquamarine3","dodgerblue3", "darkgreen", "mediumpurple3","darkorchid1", "yellowgreen", "lightgoldenrod2", "darkblue")


cell.type.colors.full <- c("#2179b4", "#b4d88b", "#36a047", "#f6999a", "#f15a29", "#fdbf6f", "#f7941d", "#cab3d6", "#6b3f98",
                           "#8c3f30", "#cdc7b1", "#2a2d7c")

p1 <- DimPlot(human.reference, reduction = "umap", group.by = "cell.type.combined", cols = cell.type.colors.full)+ xlim(c(-10,15)) + ylim(c(-10,15)) + NoLegend() + ggtitle("Human reference")
p2 <- DimPlot(macaque.query, reduction = "ref.umap", group.by = "cell.type", cols = cell.type.colors.full)+ xlim(c(-10,15)) + ylim(c(-10,15)) + NoLegend() + ggtitle("Macaque xrTestis")

p1 + p2


DimPlot(HS_NCG, group.by = "experiment")
Idents(HS_NCG) <- "experiment"
NCG_only <- subset(HS_NCG, idents = c("hPGCLC", "transplant"))

p3 <- DimPlot(NCG_only, group.by = "cell.type.combined", cols = cell.type.colors.full)+ xlim(c(-10,15)) + ylim(c(-10,15)) + NoLegend() + ggtitle("Human xrTestis")


### ### ### ###
## Figure 6c ##
### ### ### ###
p1 + p3 + p2


plot_grid(p1, p3, p2, ncol = 3)

save(human.reference, file = "human.reference.Robj")
save(macaque.query, file = "macaque.query.aligned.to.human.reference.Robj")
load("human.reference.Robj")
load("macaque.query.aligned.to.human.reference.Robj")

#banana4


# # #

# merge smples

human.reference2 <- human.reference
macaque.query2 <- macaque.query

macaque.query2@reductions$umap2 <- macaque.query@reductions$ref.umap
human.reference2@reductions$umap2 <- human.reference2@reductions$umap

merged.seurat <- merge(human.reference2, macaque.query2, merge.dr = NA)
DimPlot(merged.seurat, reduction = "umap2.1", group.by = "species")

# # # 








DimPlot(current.integrated, label = T,
        group.by = "NCG_germ_macaque3_temp_clusters")
# group.by = "experiment"
# group.by = "species2"
# group.by = "NCG_germ_macaque3_temp_clusters"
# group.by = "macaque3_clusters")

DimPlot(current.integrated)

p1 <- DimPlot(object = current.integrated, reduction = "umap", group.by = "Phase", label = T)
p2 <- DimPlot(object = current.integrated, reduction = "umap", group.by = "cell.type", label = T)

plot_grid(p1, p2)


FeaturePlot(current.integrated, features = 'nFeature_RNA', min.cutoff = 0)



DimPlot(current.integrated, label = T,
        split.by = "species",
        group.by = "species2")
FeaturePlot(current.integrated, #gene.list
            features = gene.list, order = T, 
            # split.by = "species",
            min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")

FeaturePlot(current.integrated, #gene.list
            features = c("MEIOB"),#, "SYCP1", "SYCP3", "MEIOC", "PRDM9", "DMC1", "SPO11"), 
            order = T, 
            # split.by = "species2",
            min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")


DimPlot(current.integrated, group.by = "NCG_germ_macaque3_temp_clusters", label = T)
# DimPlot(current.integrated, label = T)

# HS_NCG <- current.integrated
# HS_NCG_straight_integrated <- current.integrated

NCG_NCG_monkey <- current.integrated
setwd("/Users/ewhelan/Documents/scRNAseq/Rscripts/savefiles")
# save(HS_NCG_monkey, file = "HS_NCG_NCGmonkey_all_merged.2023.12.03.Robj")
save(NCG_NCG_monkey, file = "NCG_NCGmonkey_old_integration_method.2024.01.03.Robj")

# save(HS_NCG_original, file = "HS_NCG_original.2023.07.07.Robj")
load("HS_NCG_original.2023.07.07.Robj")
load("HS_NCG_NCGmonkey_all_merged.2023.12.03.Robj")
HS_NCG <- HS_NCG_original


p1 <- (DimPlot(HS_NCG_monkey, label = T,
               group.by = "species2"))
p2 <- FeaturePlot(HS_NCG_monkey, #gene.list
                  features = c("MEIOB", "STRA8"),#, "SYCP1", "SYCP3", "MEIOC", "PRDM9", "DMC1", "SPO11"), 
                  order = T, 
                  # split.by = "species2",
                  min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")
print(plot_grid(p1, p2))

FeaturePlot(HS_NCG_monkey, #gene.list
            features = c("PIWIL4", "STRA8", "SOHLH2", "REC8"),#, "SYCP1", "SYCP3", "MEIOC", "PRDM9", "DMC1", "SPO11"), 
            order = T, 
            split.by = "species2",
            min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")


DefaultAssay(HS_NCG) <- "RNA"

p1 <- FeaturePlot(HS_NCG, features = 'Apoptotic1', min.cutoff = 0)
p2 <- DimPlot(HS_NCG, label = T)
plot_grid(p1, p2)
#remove cluster 60/72 as showing a lot of apoptotic cells and cluster 96 as low-quality cells
HS_NCG <- subset(HS_NCG, idents = c(60, 72, 96), invert = T)
# for(x in 20:40){
DefaultAssay(object = HS_NCG) <- "integrated"
n.pcs <- 31 #34 was pretty good! #31 is also good #36 #35
HS_NCG <- RunPCA(object = HS_NCG,npcs = n.pcs,verbose = TRUE)
HS_NCG <- FindNeighbors(object = HS_NCG)
HS_NCG <- FindClusters(object = HS_NCG, resolution = 1.1, verbose = TRUE)
HS_NCG <- RunUMAP(object = HS_NCG, reduction = "pca", dims = 1:n.pcs)
print(DimPlot(HS_NCG, group.by = "cell.type", label = T)+ggtitle(n.pcs))
# }




DefaultAssay(HS_NCG) <- "RNA"

HS_NCG <- AddModuleScore(HS_NCG,features = list(apoptotic.genes),name="Apoptotic")
HS_NCG <- AddModuleScore(HS_NCG,features = list(gene.list.PGC.E),name="PGCE")
HS_NCG <- AddModuleScore(HS_NCG,features = list(gene.list.PGC.L),name="PGCL")
HS_NCG <- AddModuleScore(HS_NCG,features = list(gene.list.M),name="M")
HS_NCG <- AddModuleScore(HS_NCG,features = list(gene.list.M2T1),name="M2T")
HS_NCG <- AddModuleScore(HS_NCG,features = list(gene.list.T1),name="T")
HS_NCG <- AddModuleScore(HS_NCG,features = list(gene.list.spermatogonia),name="spg")
HS_NCG <- AddModuleScore(HS_NCG,features = list(gene.list.diff.spermatgonia.E),name="diffspgE")
HS_NCG <- AddModuleScore(HS_NCG,features = list(gene.list.diff.spermatgonia.L),name="diffspgL")
HS_NCG <- AddModuleScore(HS_NCG,features = list(gene.list.preleptotene),name="preleptotene")
HS_NCG <- AddModuleScore(HS_NCG,features = list(gene.list.leptotene),name="leptotene")
HS_NCG <- AddModuleScore(HS_NCG,features = list(gene.list.zygotene),name="zygotene")
HS_NCG <- AddModuleScore(HS_NCG,features = list(gene.list.late.meiosis),name="pachytene_2ndary")
HS_NCG <- AddModuleScore(HS_NCG,features = list(gene.list.spermatids),name="spermatids")

HS_NCG <- AddModuleScore(HS_NCG, features = list(gene.list.hermann.SSC), name="hermann.sscs")
HS_NCG <- AddModuleScore(HS_NCG, features = list(gene.list.hermann.progenitor), name="hermann.progenitor")
HS_NCG <- AddModuleScore(HS_NCG, features = list(gene.list.hermann.earlydiff), name="hermann.earlydiff")
HS_NCG <- AddModuleScore(HS_NCG, features = list(gene.list.hermann.latediff), name="hermann.latediff")
HS_NCG <- AddModuleScore(HS_NCG, features = list(gene.list.hermann.prelep), name="hermann.prelep")
HS_NCG <- AddModuleScore(HS_NCG, features = list(gene.list.hermann.lepzyg), name="hermann.lep_zyg")
HS_NCG <- AddModuleScore(HS_NCG, features = list(gene.list.hermann.pachytene), name="hermann.pachytene")
HS_NCG <- AddModuleScore(HS_NCG, features = list(gene.list.hermann.dipsecondary), name="hermann.diplotene_secondary")
HS_NCG <- AddModuleScore(HS_NCG, features = list(gene.list.hermann.earlyround), name="hermann.earlyround")
HS_NCG <- AddModuleScore(HS_NCG, features = list(gene.list.hermann.midround), name="hermann.midround")
HS_NCG <- AddModuleScore(HS_NCG, features = list(gene.list.hermann.lateround), name="hermann.lateround")


plot_grid(ncol = 4,
          DimPlot(HS_NCG, group.by = "seurat_clusters", label = T) + NoLegend(),
          FeaturePlot(HS_NCG,features = "PGCE1", label = F, min.cutoff = 0,repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG,features = "PGCL1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG,features = "M1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG,features = "M2T1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG,features = "T1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG,features = "spg1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG,features = "hermann.progenitor1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG,features = "hermann.earlydiff1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG,features = "diffspgE1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG,features = "hermann.latediff1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG,features = "diffspgL1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG,features = "preleptotene1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG, features = "hermann.prelep1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG,features = "leptotene1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG,features = "zygotene1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG, features = "hermann.lep_zyg1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG, features = "hermann.pachytene1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG, features = "hermann.diplotene_secondary1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          # FeaturePlot(HS_NCG, features = "hermann.earlyround1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          # FeaturePlot(HS_NCG, features = "hermann.midround1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          # FeaturePlot(HS_NCG, features = "hermann.lateround1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme(),
          FeaturePlot(HS_NCG, features = "spermatids1", label = F, repel = FALSE) + scale_color_viridis(option="plasma") + DarkTheme()
)


### monocle pseudotime ###

# cds <- as.cell_data_set(NCG_germ4subset)
cds <- as.cell_data_set(HS_NCG)


cds <- cluster_cells(cds, reduction_method = c("UMAP"), #cluster_method = "louvain",
                     resolution = 1.7e-4, 
                     k = 8,
                     partition_qval = 0.05
)

plot_cells(cds, show_trajectory_graph = FALSE)

cds <- learn_graph(cds, 
                   close_loop = F,
                   learn_graph_control=list(ncenter=450,
                                            # nn.method = "nn2",
                                            minimal_branch_len=9), #450, 8 is pretty good
                   use_partition = F)
print(plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)) #+ ggtitle(a))

# Figure 2C #######################################
cds <- order_cells(cds)#, root_cells = max.avp)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)+ NoLegend()

cds.pseudotime <- pseudotime(cds, reduction_method = "UMAP")
cds.clusters <- clusters(cds, reduction_method = "UMAP")



#now add the pseudotime to the original object

HS_NCG$pseudotime <- cds.pseudotime
HS_NCG$monocle_clusters <- cds.clusters


plot_grid(p1, p2)
Idents(HS_NCG) <- "cell.type"
DimPlot(HS_NCG, 
        # split.by = "experiment", 
        group.by = "cell.type", label = T)


list_of_cols <- rep("gray", each = 14)
list_of_cols[6] <- "blue"
list_of_cols[5] <- "magenta"
list_of_cols[7] <- "turquoise"
Idents(HS_NCG) <- "cell.type"
DimPlot(HS_NCG, 
        split.by = "experiment",
        cols = list_of_cols,
        label = T)

DimPlot(HS_NCG, group.by = "Phase")

# save(HS_NCG, file = "HS_NCG.early.cutoff.T1.2023.07.07.Robj")


# set new cell idents

p1 <- DimPlot(HS_NCG, group.by = "seurat_clusters", label = T)
p2 <- DimPlot(HS_NCG, group.by = "monocle_clusters", label = T)
plot_grid(p1, p2)


Idents(HS_NCG) <- HS_NCG$experiment
DimPlot(HS_NCG)
table(Idents(HS_NCG))
HS_NCG <- SeuratObject::RenameIdents(HS_NCG, 'hPGCLC' = 'transplant')

Idents(HS_NCG) -> HS_NCG$experiment2
FeaturePlot(HS_NCG, features = "pseudotime", split.by = "experiment2") & scale_color_viridis(option = "C")


Idents(HS_NCG) <- HS_NCG$monocle_clusters
cell.type.colors <- c("#2179b4", "#b4d88b", "#36a047", "#f6999a", "#f15a29", "#fdbf6f", "#f7941d", "#cab3d6", "#6b3f98")




new.cluster.ids <- character(length(levels(Idents(HS_NCG))))

PGCE.ids <- c(6, 17)
PGCL.ids <- c(1, 4, 15, 14)
M1.ids <- c(5, 18, 8)
M2T1.ids <- c(3)
T1.ids <-c(10)
spg.ids <- c(7, 12, 2)
diffE.ids <- c(13)
diffL.ids <-c(9)
prelep.ids <- c(16)
# lepzyg.ids <- c(27, 20, 24)
pachy.ids <-c(22, 21)
dipsec.ids <-c(20,19)
spermatid.ids <- c(11)

new.cluster.ids[PGCE.ids] <- "PGCs (early)"
new.cluster.ids[PGCL.ids] <- "PGCs (late)"
new.cluster.ids[M1.ids] <- "M prospermatogonia"
new.cluster.ids[M2T1.ids] <- "M/T1 prospermaogonia"
new.cluster.ids[T1.ids] <- "T1 prospermatogonia"
new.cluster.ids[spg.ids] <- "Spermatogonia"
new.cluster.ids[diffE.ids] <- "Diff. spermatogonia (early)"
new.cluster.ids[diffL.ids] <- "Diff. spermatogonia (late)"
new.cluster.ids[prelep.ids] <- "Prelep/Lep/Zygotene"
# new.cluster.ids[lepzyg.ids] <- "Leptotene/Zygotene"
new.cluster.ids[pachy.ids] <- "Pachytene spermatocytes"
new.cluster.ids[dipsec.ids] <- "Diplotene/2ndary spermatocytes"
new.cluster.ids[spermatid.ids] <- "Spermatids"



new.cluster.ids

names(x = new.cluster.ids) <- levels(x = HS_NCG)
HS_NCG <- RenameIdents(object = HS_NCG, new.cluster.ids)
HS_NCG$cell.type.combined <- Idents(HS_NCG)

HS_NCG$cell.type.combined <- factor(x = HS_NCG$cell.type.combined, levels = c("PGCs (early)", 
                                                                              "PGCs (late)",
                                                                              "M prospermatogonia",
                                                                              "M/T1 prospermaogonia",
                                                                              "T1 prospermatogonia",
                                                                              "Spermatogonia",
                                                                              "Diff. spermatogonia (early)",
                                                                              "Diff. spermatogonia (late)",
                                                                              # "Preleptotene spermatocytes",
                                                                              "Prelep/Lep/Zygotene",
                                                                              "Pachytene spermatocytes",
                                                                              "Diplotene/2ndary spermatocytes",
                                                                              "Spermatids"))


Idents(HS_NCG) <- HS_NCG$cell.type.combined
DimPlot(HS_NCG, 
        group.by = "cell.type.combined",
        # group.by = "stage", 
        # split.by = "experiment",
        label = T) + NoLegend()
# save(HS_NCG, file="HS_NCG_2023.07.05.Robj")
# save(HS_NCG, file="HS_NCG_2023.07.06.Robj")


#basis for Figure 4
DimPlot(HS_NCG,label = T, group.by = "cell.type.combined", cols = c(cell.type.colors, "coral4", "cornsilk3", "blue4"))


# save(HS_NCG, file="HS_NCG_2023.07.19.Robj")
load("HS_NCG_2023.07.17.Robj")



Idents(HS_NCG) <- "experiment"
DimPlot(HS_NCG)
# HS_NCG_only_HS <- subset(HS_NCG, idents = "human_in_vivo")
HS_NCG_only_HS <- subset(HS_NCG, idents = "human in vivo")
# HS_NCG_only_NCG <- subset(HS_NCG, idents = c("transplant_human", "hPGCLC"))
HS_NCG_only_NCG <- subset(HS_NCG, idents = c("transplant", "hPGCLC"))

#figure out where undifferentiated and differentiated spermatogonia are
FeaturePlot(HS_NCG, features = c("KIT", "STRA8", "RET", "ETV5"), order = T)
FeaturePlot(HS_NCG, features = c("MAGEA4", "NANOS2"), order = T)

FeaturePlot(HS_NCG, features = c("REC8", "TEX101"), split.by = "experiment", order = T)






#PSEUDOTIME PLOT
p1 <- FeaturePlot(HS_NCG_only_HS, features = "pseudotime") + scale_color_viridis(option="plasma") + NoLegend() + ggtitle("in vivo")
p2 <- FeaturePlot(HS_NCG_only_NCG, features = "pseudotime") + scale_color_viridis(option="plasma") + NoLegend() + ggtitle("in vitro")
p3 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
                 label_branch_points = FALSE)+ NoLegend() + ggtitle("Pseudotime")

plot_grid(p3, p1, p2, ncol = 3)


p1 <- DimPlot(HS_NCG_only_HS, group.by = "stage", label = F)
p2 <- DimPlot(HS_NCG_only_NCG, group.by = "stage")
plot_grid(p1, p2)

#make a heatmap

#First, find markers

DimPlot(HS_NCG, split.by = "experiment")

Idents(HS_NCG) <- HS_NCG$cell.type.combined
Idents(HS_NCG_only_NCG) <- HS_NCG_only_NCG$cell.type.combined
Idents(HS_NCG_only_HS) <- HS_NCG_only_HS$cell.type.combined

all.markers <- FindAllMarkers(HS_NCG)
all.markers_NCG <- FindAllMarkers(HS_NCG_only_NCG)
all.markers_HS <- FindAllMarkers(HS_NCG_only_HS)

# write.csv(all.markers, file = "all.markers.HS_NCG.csv")
# write.csv(all.markers, file = "all.markers.only_NCG.csv")

DimPlot(HS_NCG_only_NCG)

HS_NCG_only_NCG_subset <- subset(HS_NCG_only_NCG, idents = c("Pachytene spermatocytes", "Diplotene/2ndary spermatocytes", "Spermatids"), invert = T)
HS_NCG_only_HS_subset <- subset(HS_NCG_only_HS, idents = c("Pachytene spermatocytes", "Diplotene/2ndary spermatocytes", "Spermatids"), invert = T)
# setwd("/Users/ewhelan/Documents/scRNAseq/Rscripts/savefiles")

# save(HS_NCG_only_NCG_subset, file = "HS_NCG_only_NCG_subset.Robj")
# save(HS_NCG_only_HS_subset, file = "HS_NCG_only_HS_subset.Robj")

load("HS_NCG_only_NCG_subset.Robj")
load("HS_NCG_only_HS_subset.Robj")

# all.markers.subset <- filter(all.markers_NCG, avg_log2FC > 1)


### ### ### ### ###
###HEATMAP NCG ###
### ### ### ### ###


PGCs.E.markers.all.ncg <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "PGCs (early)")
PGCs.L.markers.all.ncg <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "PGCs (late)")
M.ProSpermatogonia.markers.all.ncg <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "M prospermatogonia")
M2T1.markers.all.ncg <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "M/T1 prospermaogonia")
T.ProSpermatogonia.markers.all.ncg <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "T1 prospermatogonia")
Spg.markers.all.ncg <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "Spermatogonia")
Diff.Spg.E.markers.all.ncg <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "Diff. spermatogonia (early)")
Diff.Spg.L.markers.all.ncg <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "Diff. spermatogonia (late)")
Pre.Leptotene.markers.all.ncg <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "Prelep/Lep/Zygotene")

# setwd("/Users/ewhelan/Desktop/xrTestisPaper/NCG/NCG_marker_genes")
# write.csv(PGCs.E.markers.all.ncg, file="HS_NCG_only_NCG_subset_PGCs.E.markers.all.ncg.csv")
# write.csv(PGCs.L.markers.all.ncg, file="HS_NCG_only_NCG_subset_PGCs.L.markers.all.ncg.csv")
# write.csv(M.ProSpermatogonia.markers.all.ncg, file="HS_NCG_only_NCG_subset_M.ProSpermatogonia.markers.all.ncg.csv")
# write.csv(M2T1.markers.all.ncg, file="HS_NCG_only_NCG_subset_M2T1.markers.all.ncg.csv")
# write.csv(T.ProSpermatogonia.markers.all.ncg, file="HS_NCG_only_NCG_subset_T.ProSpermatogonia.markers.all.ncg.csv")
# write.csv(Spg.markers.all.ncg, file="HS_NCG_only_NCG_subset_Spg.markers.all.ncg.csv")
# write.csv(Diff.Spg.E.markers.all.ncg, file="HS_NCG_only_NCG_subset_Diff.Spg.E.markers.all.ncg.csv")
# write.csv(Diff.Spg.L.markers.all.ncg, file="HS_NCG_only_NCG_subset_Diff.Spg.L.markers.all.ncg.csv")
# write.csv(Pre.Leptotene.markers.all.ncg, file="HS_NCG_only_NCG_subset_Pre.Leptotene.markers.all.ncg.csv")


genes.per.cluster <- 300
PGCs.E.markers.ordered.ncg <- PGCs.E.markers.all.ncg %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
PGCs.L.markers.ordered.ncg <- PGCs.L.markers.all.ncg %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
M.ProSpermatogonia.markers.ordered.ncg <- M.ProSpermatogonia.markers.all.ncg %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
M2T1.markers.ordered.ncg <- M2T1.markers.all.ncg %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
T.ProSpermatogonia.markers.ordered.ncg <- T.ProSpermatogonia.markers.all.ncg %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
Spg.markers.ordered.ncg <- Spg.markers.all.ncg %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
Diff.Spg.E.markers.ordered.ncg <- Diff.Spg.E.markers.all.ncg %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
Diff.Spg.L.markers.ordered.ncg <- Diff.Spg.L.markers.all.ncg %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
Pre.Leptotene.markers.ordered.ncg <- Pre.Leptotene.markers.all.ncg %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)

PGCs.E.markers.ordered.ncg$cell.type <- "PGCS.E" 
PGCs.L.markers.ordered.ncg$cell.type <- "PGCs.L" 
M.ProSpermatogonia.markers.ordered.ncg$cell.type <- "M" 
M2T1.markers.ordered.ncg$cell.type <- "M2T1" 
T.ProSpermatogonia.markers.ordered.ncg$cell.type <- "T"
Spg.markers.ordered.ncg$cell.type <- "Spg"
Diff.Spg.E.markers.ordered.ncg$cell.type <- "Diff.Spg.E"
Diff.Spg.L.markers.ordered.ncg$cell.type <- "Diff.Spg.L"
Pre.Leptotene.markers.ordered.ncg$cell.type <- "Prelep.lep"

topNmarkers.ncg <- rbind(PGCs.E.markers.ordered.ncg, PGCs.L.markers.ordered.ncg, M.ProSpermatogonia.markers.ordered.ncg, M2T1.markers.ordered.ncg,
                         T.ProSpermatogonia.markers.ordered.ncg, Spg.markers.ordered.ncg, Diff.Spg.E.markers.ordered.ncg, Diff.Spg.L.markers.ordered.ncg, Pre.Leptotene.markers.ordered.ncg)
topNmarkers.ncg.df <- as.data.frame(topNmarkers.ncg)
topNmarkers.ncg.df$gene <- rownames(topNmarkers.ncg.df)

# Heatmap of select genes by cluster ----

Idents(HS_NCG_only_NCG) <- HS_NCG_only_NCG$cell.type.combined
DimPlot(HS_NCG_only_NCG)

data.pgcE.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "PGCs (early)")])
data.pgcL.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "PGCs (late)")])
data.m.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "M prospermatogonia")])
data.mt1.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "M/T1 prospermaogonia")])
data.t1.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "T1 prospermatogonia")])
data.undiff.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "Spermatogonia")])
data.diffE.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "Diff. spermatogonia (early)")])
data.diffL.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "Diff. spermatogonia (late)")])
data.Spermatocytes.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "Prelep/Lep/Zygotene")])



av.counts.pgcE.ncg <- apply(data.pgcE.ncg, 1, mean)
av.counts.pgcL.ncg <- apply(data.pgcL.ncg, 1, mean)
av.counts.m.ncg <- apply(data.m.ncg, 1, mean)
av.counts.mt1.ncg <- apply(data.mt1.ncg, 1, mean)
av.counts.t1.ncg <- apply(data.t1.ncg, 1, mean)
av.counts.undiff.ncg <- apply(data.undiff.ncg, 1, mean)
av.counts.diffE.ncg <- apply(data.diffE.ncg, 1, mean)
av.counts.diffL.ncg <- apply(data.diffL.ncg, 1, mean)
av.counts.Spermatocytes.ncg <- apply(data.Spermatocytes.ncg, 1, mean)


av.counts.HS_NCG_only_NCG <- cbind(av.counts.pgcE.ncg, av.counts.pgcL.ncg, av.counts.m.ncg, av.counts.mt1.ncg,av.counts.t1.ncg,
                                   av.counts.undiff.ncg, av.counts.diffE.ncg, av.counts.diffL.ncg, av.counts.Spermatocytes.ncg)


av.counts.HS_NCG_only_NCG.df <- as.data.frame(av.counts.HS_NCG_only_NCG)

gene.name.list.ncg <- rownames(av.counts.HS_NCG_only_NCG.df)

av.counts.HS_NCG_only_NCG.df$gene.name <- gene.name.list.ncg

# genes.of.interest <- Keren.gene.list[c(1:48)] #this is just Keren's list of genes
# genes.of.interest <- (all.markers.subset$gene) #this is a log-fold cutoff
genes.of.interest.ncg <- rownames(topNmarkers.ncg) #this is just taking the top N markers

av.counts.genes.of.interest.ncg <- av.counts.HS_NCG_only_NCG.df %>%
  dplyr::filter(gene.name %in% genes.of.interest.ncg)


av.counts.genes.of.interest.ordered.ncg <- av.counts.genes.of.interest.ncg[match(genes.of.interest.ncg, av.counts.genes.of.interest.ncg$gene.name),]

integer.col.ncg <- ncol(av.counts.genes.of.interest.ordered.ncg)

av.counts.genes.of.interest.ordered.ncg <- av.counts.genes.of.interest.ordered.ncg[,-integer.col.ncg]

av.counts.genes.of.interest.matrix1.ncg <- as.matrix(av.counts.genes.of.interest.ordered.ncg)
colnames(av.counts.genes.of.interest.matrix1.ncg) <- c("PGCs (early)", "PGCs (late)", "M prospermatogonia", 
                                                       "M/T1 prospermaogonia", "T1 prospermatogonia", "Spermatogonia", 
                                                       "Diff. spermatogonia (early)", "Diff. spermatogonia (late)", "Prelep/Lep/Zygotene")

av.counts.genes.of.interest.matrix.no.na.ncg <- na.omit(av.counts.genes.of.interest.matrix1.ncg)



l.mat <- rbind(c(0,4), c(0, 3), c(2,1), c(0,0))
l.mat

# Figure 3A ----

#NOTE select matrix1 or matrix2 below to change between unselected and EpCAM


# tiff(file = "plot.tiff",width = 5000, height = 100000)

heatmap.ncg <- heatmap.2(av.counts.genes.of.interest.matrix.no.na.ncg, Colv = FALSE,
                         Rowv = FALSE,
                         trace = "none",
                         # col =  colorRampPalette(colors=c("blue","white","red"))(100),
                         col = inferno(100),
                         density.info="none",
                         srtCol = -60, #angle of labels
                         adjCol = c(0,0), #offset of labels, should be centered now
                         cexRow=1.0, cexCol=1.3,
                         lmat= l.mat,
                         lwid = c(0.5, 5),
                         lhei=c(1, 1, 5, 1),
                         scale = "row"
                         
)

dev.off()
# 
# ### DOT PLOT XRTESTIS ###
# 
# genes.of.interest <- c("POU5F1", 
#                        "NANOG",
#                        "TCL1A")
# 
# dotplot2 <- ggplot(data = all.markers, na.rm = T,
#                    aes(x = cluster, y = (gene), 
#                        color = avg_log2FC, size = pct.1))+ 
#   geom_point() +
#   # scale_color_gradient(high = "red" , low = "gray") +
#   scale_color_viridis(option = "inferno")+
#   scale_size(range = c(-7, 7)) +
#   theme_bw() + 
#   ylab("Gene") + 
#   xlab("Cell Type") + 
#   ggtitle("xrTestis") +
#   scale_y_discrete(position = "right")+
#   theme(axis.text.x = element_text( #face = "bold",
#     color = "black", size = 14, angle = -60, hjust = 0
#   ),
#   axis.text.y = element_text( #face = "bold",
#     color = "black",  hjust = 1)
#   ) +
#   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
#   theme(legend.position = "left")
# 





### ### ### ### ###
###HEATMAP HS ###
### ### ### ### ###


PGCs.E.markers.all.hs <- FindMarkers(HS_NCG_only_HS_subset, ident.1 = "PGCs (early)")
PGCs.L.markers.all.hs <- FindMarkers(HS_NCG_only_HS_subset, ident.1 = "PGCs (late)")
M.ProSpermatogonia.markers.all.hs <- FindMarkers(HS_NCG_only_HS_subset, ident.1 = "M prospermatogonia")
M2T1.markers.all.hs <- FindMarkers(HS_NCG_only_HS_subset, ident.1 = "M/T1 prospermaogonia")
T.ProSpermatogonia.markers.all.hs <- FindMarkers(HS_NCG_only_HS_subset, ident.1 = "T1 prospermatogonia")
Spg.markers.all.hs <- FindMarkers(HS_NCG_only_HS_subset, ident.1 = "Spermatogonia")
Diff.Spg.E.markers.all.hs <- FindMarkers(HS_NCG_only_HS_subset, ident.1 = "Diff. spermatogonia (early)")
Diff.Spg.L.markers.all.hs <- FindMarkers(HS_NCG_only_HS_subset, ident.1 = "Diff. spermatogonia (late)")
Pre.Leptotene.markers.all.hs <- FindMarkers(HS_NCG_only_HS_subset, ident.1 = "Prelep/Lep/Zygotene")

# setwd("/Users/ewhelan/Desktop/xrTestisPaper/NCG/NCG_marker_genes")
# write.csv(PGCs.E.markers.all.hs, file="HS_NCG_only_HS_subset_PGCs.E.markers.all.hs.csv")
# write.csv(PGCs.L.markers.all.hs, file="HS_NCG_only_HS_subset_PGCs.L.markers.all.hs.csv")
# write.csv(M.ProSpermatogonia.markers.all.hs, file="HS_NCG_only_HS_subset_M.ProSpermatogonia.markers.all.hs.csv")
# write.csv(M2T1.markers.all.hs, file="HS_NCG_only_HS_subset_M2T1.markers.all.hs.csv")
# write.csv(T.ProSpermatogonia.markers.all.hs, file="HS_NCG_only_HS_subset_T.ProSpermatogonia.markers.all.hs.csv")
# write.csv(Spg.markers.all.hs, file="HS_NCG_only_HS_subset_Spg.markers.all.hs.csv")
# write.csv(Diff.Spg.E.markers.all.hs, file="HS_NCG_only_HS_subset_Diff.Spg.E.markers.all.hs.csv")
# write.csv(Diff.Spg.L.markers.all.hs, file="HS_NCG_only_HS_subset_Diff.Spg.L.markers.all.hs.csv")
# write.csv(Pre.Leptotene.markers.all.hs, file="HS_NCG_only_HS_subset_Pre.Leptotene.markers.all.hs.csv")

genes.per.cluster <- 300
PGCs.E.markers.ordered.hs <- PGCs.E.markers.all.hs %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
PGCs.L.markers.ordered.hs <- PGCs.L.markers.all.hs %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
M.ProSpermatogonia.markers.ordered.hs <- M.ProSpermatogonia.markers.all.hs %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
M2T1.markers.ordered.hs <- M2T1.markers.all.hs %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
T.ProSpermatogonia.markers.ordered.hs <- T.ProSpermatogonia.markers.all.hs %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
Spg.markers.ordered.hs <- Spg.markers.all.hs %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
Diff.Spg.E.markers.ordered.hs <- Diff.Spg.E.markers.all.hs %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
Diff.Spg.L.markers.ordered.hs <- Diff.Spg.L.markers.all.hs %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
Pre.Leptotene.markers.ordered.hs <- Pre.Leptotene.markers.all.hs %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)

PGCs.E.markers.ordered.hs$cell.type <- "PGCS.E" 
PGCs.L.markers.ordered.hs$cell.type <- "PGCs.L" 
M.ProSpermatogonia.markers.ordered.hs$cell.type <- "M" 
M2T1.markers.ordered.hs$cell.type <- "M2T1" 
T.ProSpermatogonia.markers.ordered.hs$cell.type <- "T"
Spg.markers.ordered.hs$cell.type <- "Spg"
Diff.Spg.E.markers.ordered.hs$cell.type <- "Diff.Spg.E"
Diff.Spg.L.markers.ordered.hs$cell.type <- "Diff.Spg.L"
Pre.Leptotene.markers.ordered.hs$cell.type <- "Prelep.lep"

topNmarkers.hs <- rbind(PGCs.E.markers.ordered.hs, PGCs.L.markers.ordered.hs, M.ProSpermatogonia.markers.ordered.hs, M2T1.markers.ordered.hs,
                        T.ProSpermatogonia.markers.ordered.hs, Spg.markers.ordered.hs, Diff.Spg.E.markers.ordered.hs, Diff.Spg.L.markers.ordered.hs, Pre.Leptotene.markers.ordered.hs)



# Heatmap of select genes by cluster ----

Idents(HS_NCG_only_HS) <- HS_NCG_only_HS$cell.type.combined
DimPlot(HS_NCG_only_HS)

data.pgcE.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "PGCs (early)")])
data.pgcL.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "PGCs (late)")])
data.m.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "M prospermatogonia")])
data.mt1.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "M/T1 prospermaogonia")])
data.t1.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "T1 prospermatogonia")])
data.undiff.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "Spermatogonia")])
data.diffE.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "Diff. spermatogonia (early)")])
data.diffL.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "Diff. spermatogonia (late)")])
data.Spermatocytes.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "Prelep/Lep/Zygotene")])


av.counts.pgcE.hs <- apply(data.pgcE.hs, 1, mean)
av.counts.pgcL.hs <- apply(data.pgcL.hs, 1, mean)
av.counts.m.hs <- apply(data.m.hs, 1, mean)
av.counts.mt1.hs <- apply(data.mt1.hs, 1, mean)
av.counts.t1.hs <- apply(data.t1.hs, 1, mean)
av.counts.undiff.hs <- apply(data.undiff.hs, 1, mean)
av.counts.diffE.hs <- apply(data.diffE.hs, 1, mean)
av.counts.diffL.hs <- apply(data.diffL.hs, 1, mean)
av.counts.Spermatocytes.hs <- apply(data.Spermatocytes.hs, 1, mean)


av.counts.HS_NCG_only_HS <- cbind(av.counts.pgcE.hs, av.counts.pgcL.hs, av.counts.m.hs, av.counts.mt1.hs,av.counts.t1.hs,
                                  av.counts.undiff.hs, av.counts.diffE.hs, av.counts.diffL.hs, av.counts.Spermatocytes.hs)


av.counts.HS_NCG_only_HS.df <- as.data.frame(av.counts.HS_NCG_only_HS)

gene.name.list.hs <- rownames(av.counts.HS_NCG_only_HS.df)

av.counts.HS_NCG_only_HS.df$gene.name <- gene.name.list.hs

# genes.of.interest <- Keren.gene.list[c(1:48)] #this is just Keren's list of genes
# genes.of.interest <- (all.markers.subset$gene) #this is a log-fold cutoff
genes.of.interest.hs <- rownames(topNmarkers.hs) #this is just taking the top N markers

av.counts.genes.of.interest.hs <- av.counts.HS_NCG_only_HS.df %>%
  dplyr::filter(gene.name %in% genes.of.interest.hs)


av.counts.genes.of.interest.ordered.hs <- av.counts.genes.of.interest.hs[match(genes.of.interest.hs, av.counts.genes.of.interest.hs$gene.name),]

integer.col.hs <- ncol(av.counts.genes.of.interest.ordered.hs)

av.counts.genes.of.interest.ordered.hs <- av.counts.genes.of.interest.ordered.hs[,-integer.col.hs]

av.counts.genes.of.interest.matrix1.hs <- as.matrix(av.counts.genes.of.interest.ordered.hs)
colnames(av.counts.genes.of.interest.matrix1.hs) <- c("PGCs (early)", "PGCs (late)", "M prospermatogonia", 
                                                      "M/T1 prospermaogonia", "T1 prospermatogonia", "Spermatogonia", 
                                                      "Diff. spermatogonia (early)", "Diff. spermatogonia (late)", "Prelep/Lep/Zygotene")

av.counts.genes.of.interest.matrix.no.na.hs <- na.omit(av.counts.genes.of.interest.matrix1.hs)



l.mat <- rbind(c(0,4), c(0, 3), c(2,1), c(0,0))
l.mat

# Figure 3A ----

#NOTE select matrix1 or matrix2 below to change between unselected and EpCAM


# tiff(file = "plot.tiff",width = 5000, height = 100000)

heatmap.hs <- heatmap.2(av.counts.genes.of.interest.matrix.no.na.hs, Colv = FALSE,
                        Rowv = FALSE,
                        trace = "none",
                        # col =  colorRampPalette(colors=c("blue","white","red"))(100),
                        col = inferno(100),
                        density.info="none",
                        srtCol = -60, #angle of labels
                        adjCol = c(0,0), #offset of labels, should be centered now
                        cexRow=1.0, cexCol=1.3,
                        lmat= l.mat,
                        lwid = c(0.5, 5),
                        lhei=c(1, 1, 5, 1),
                        scale = "row"
                        
)

dev.off()


# 
# library(reshape2)
# library(ggplot2)
# library(scales)
# library(hrbrthemes)
# dat <- melt(av.counts.genes.of.interest.matrix1)
# 
# ggplot(dat, aes(x=Var2,y=Var1))+ 
#   geom_tile(aes(fill = value),colour='black') +
#   scale_fill_gradientn(colours=c("green","gray","red"),
#                        values  = rescale(c(min(dat$value), 1000, max(dat$value)))) +
#   coord_equal() 
# 
# dev.off()


# pineapple

genes <- c("POU5F1", 
           "NANOG", 
           # "TCL1A", 
           "PDPN", 
           # "PRDM1", 
           "NANOS3", 
           "TFAP2C", 
           # "SOX17", 
           # "L1TD1",
           # "UTF1", 
           "KIT",
           "ASB9", 
           # "DGKK", 
           "CITED2", 
           # "GAGE12H",
           "NANOS2", 
           "PIWIL4", 
           # "MORC1", 
           # "MAGEB2", 
           # "TDRD1", 
           "RHOXF1", 
           "MAGEC2", 
           # "ZBTB16", 
           "EGR4", 
           "SIX1", 
           # "EOMES", 
           # "TEKT1", 
           # "MAGEB2", 
           # "DMRT3",
           # "MSL3", 
           # "DUSP6", 
           # "ETV5", 
           # "RET", 
           # "MAGEA4",
           "SOHLH1", 
           # "DMRT1", 
           # "SOHLH2", 
           # "ESX1", 
           # "DMRTB1", 
           # "STRA8", 
           # "DAZL",
           # "SCML2",
           # "GAL3ST1",
           
           # "PRDM9",
           # "ZCWPW1",
           
           # "SYCP3",
           
           # "PEG10",
           
           # "ATR",
           "HORMAD1",
           "MEIOC", 
           "TEX101", 
           
           
           "MEIOB",
           # "SYCP1",
           "DMC1", 
           "MYBL1", 
           "SPO11",
           
           "HORMAD2",
           "RAD51AP2",
           "TEX19")

length(genes)
Idents(HS_NCG_only_NCG) <- "cell.type.combined"
# DimPlot(HS_NCG_only_NCG)
HS_NCG_only_NCG_subset <- subset(HS_NCG_only_NCG, idents = c("Pachytene spermatocytes", "Diplotene/2ndary spermatocytes", "Spermatids"), invert = T)
# DimPlot(HS_NCG_only_NCG_subset)
# all.markers <- FindAllMarkers(HS_NCG_only_NCG_subset, min.pct = 0, logfc.threshold = 0, features = genes)
PGCs.E.markers <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "PGCs (early)", min.pct = 0, logfc.threshold = 0, features = genes)
PGCs.L.markers <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "PGCs (late)", min.pct = 0, logfc.threshold = 0, features = genes)
M.ProSpermatogonia.markers <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "M prospermatogonia", min.pct = 0, logfc.threshold = 0, features = genes)
M2T1.markers <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "M/T1 prospermaogonia", min.pct = 0, logfc.threshold = 0, features = genes)
T.ProSpermatogonia.markers <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "T1 prospermatogonia", min.pct = 0, logfc.threshold = 0, features = genes)
Spg.markers <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "Spermatogonia", min.pct = 0, logfc.threshold = 0, features = genes)
Diff.Spg.E.markers <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "Diff. spermatogonia (early)", min.pct = 0, logfc.threshold = 0, features = genes)
Diff.Spg.L.markers <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "Diff. spermatogonia (late)", min.pct = 0, logfc.threshold = 0, features = genes)
Pre.Leptotene.markers <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "Prelep/Lep/Zygotene", min.pct = 0, logfc.threshold = 0, features = genes)


PGCs.E.markers.ordered <- PGCs.E.markers %>% 
  arrange(avg_log2FC) %>% 
  slice_head(n = 50)


# graph.data <- graph.data %>%
#   arrange(Celltype)

order.of.genes <- as.vector(genes)

all.markers_NCG <- FindAllMarkers(HS_NCG_only_NCG_subset, min.pct = 0, logfc.threshold = 0, features = genes)
all.markers_HS <- FindAllMarkers(HS_NCG_only_HS_subset, min.pct = 0, logfc.threshold = 0, features = genes)


# graph.data$Pathway2 <- factor(graph.data$Pathway, levels = rev(unique(graph.data$Pathway)))
all.markers$gene <- factor(all.markers$gene, levels = rev(genes))
all.markers$cluster <- factor(all.markers$cluster, levels = c("PGCs (early)", 
                                                              "PGCs (late)", 
                                                              "M prospermatogonia",
                                                              "M/T1 prospermaogonia", 
                                                              "T1 prospermatogonia", 
                                                              "Spermatogonia",
                                                              "Diff. spermatogonia (early)", 
                                                              "Diff. spermatogonia (late)", 
                                                              "Prelep/Lep/Zygotene"))

all.markers_HS$gene  <- factor(all.markers_HS$gene, levels = rev(genes))
all.markers_HS$cluster <- factor(all.markers_HS$cluster, levels = c("PGCs (early)", 
                                                                    "PGCs (late)", 
                                                                    "M prospermatogonia",
                                                                    "M/T1 prospermaogonia", 
                                                                    "T1 prospermatogonia", 
                                                                    "Spermatogonia",
                                                                    "Diff. spermatogonia (early)", 
                                                                    "Diff. spermatogonia (late)", 
                                                                    "Prelep/Lep/Zygotene"))

all.markers_NCG$gene <- factor(all.markers_NCG$gene, levels = rev(genes))
all.markers_NCG$cluster <- factor(all.markers_NCG$cluster, levels = c("PGCs (early)", 
                                                                      "PGCs (late)", 
                                                                      "M prospermatogonia",
                                                                      "M/T1 prospermaogonia", 
                                                                      "T1 prospermatogonia", 
                                                                      "Spermatogonia",
                                                                      "Diff. spermatogonia (early)", 
                                                                      "Diff. spermatogonia (late)", 
                                                                      "Prelep/Lep/Zygotene"))


dotplot_xrtestis <- ggplot(data = all.markers_NCG, na.rm = T,
                           aes(x = cluster, y = (gene), 
                               color = avg_log2FC, size = pct.1))+ 
  geom_point() +
  # scale_color_gradient(high = "red" , low = "blue") +
  # scale_color_gradient2(high = "#f15a29" , mid = "#fdbf6f", low = "#6b3f98", midpoint = 0, limits = c(-4, 4), oob = squish) +
  # scale_color_scico(palette = "bam", limits = c(-4, 4), oob = squish)+
  # scale_color_scico(palette = "broc", limits = c(-4, 4), oob = squish)+
  scale_color_scico(palette = "vikO", limits = c(-4, 4), midpoint = 0, oob = squish)+
  # scale_color_viridis(option = "inferno", limits = c(-4, 4), oob = squish)+
  # scale_color_distiller(palette = "RdYlBu", limicts = c(-4, 4), oob = squish)+
  scale_size(range = c(-7, 7), limits = c(0, 1)) +
  theme_bw() + 
  ylab("Gene") + 
  xlab("Cell Type") + 
  ggtitle("xrTestis") +
  scale_y_discrete(position = "right")+
  theme(axis.text.x = element_text( #face = "bold",
    color = "black", size = 14, angle = -60, hjust = 0
  ),
  axis.text.y = element_text( #face = "bold",
    color = "black",  hjust = 1)
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "left")


dotplot_xrtestis

dotplot_xrtestis_black <- ggplot(data = all.markers_NCG, na.rm = T,
                                 aes(x = cluster, y = (gene), 
                                     color = avg_log2FC, size = pct.1))+ 
  geom_point() +
  # scale_color_gradient(high = "red" , low = "blue") +
  scale_color_viridis(option = "inferno", limits = c(-4, 4), oob = squish)+
  scale_size(range = c(-7, 7), limits = c(0, 1)) +
  theme_bw() + 
  ylab("Gene") + 
  xlab("Cell Type") + 
  ggtitle("xrTestis") +
  scale_y_discrete(position = "right")+
  theme(axis.text.x = element_text( #face = "bold",
    color = "black", size = 14, angle = -60, hjust = 0
  ),
  axis.text.y = element_text( #face = "bold",
    color = "black",  hjust = 1)
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "left")+
  theme(panel.background = element_rect(fill = "black"))




dotplot_inVivo <- ggplot(data = all.markers_HS, na.rm = T,
                         aes(x = cluster, y = (gene),
                             color = avg_log2FC, size = pct.1))+
  geom_point() +
  # scale_color_gradient2(high = "#f15a29" , mid = "#fdbf6f", low = "#6b3f98", midpoint = 0, limits = c(-4, 4), oob = squish) +
  # scale_color_gradient2(high = "red" , mid = "orange", low = "blue", midpoint = 0, limits = c(-4, 4), oob = squish) +
  # scale_color_scico(palette = "bam", limits = c(-4, 4), oob = squish)+
  # scale_color_scico(palette = "broc", limits = c(-4, 4), oob = squish)+
  scale_color_scico(palette = "vikO", limits = c(-4, 4), midpoint = 0, oob = squish)+
  # scale_color_viridis(option = "inferno", limits = c(-4, 4), oob = squish)+
  # scale_color_distiller(palette = "RdYlBu", limits = c(-4, 4), oob = squish)+
  scale_size(range = c(-7, 7), limits = c(0, 1)) +
  theme_bw() +
  ylab("Gene") +
  xlab("Cell Type") +
  ggtitle("In Vivo") +
  scale_y_discrete(position = "right")+
  theme(axis.text.x = element_text( #face = "bold",
    color = "black", size = 14, angle = -60, hjust = 0
  ),
  axis.text.y = element_text( #face = "bold",
    color = "black",  hjust = 1)
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "left")




dotplot_inVivo_black <- ggplot(data = all.markers_HS, na.rm = TRUE,
                               aes(x = cluster, y = gene, 
                                   color = avg_log2FC, size = pct.1))+ 
  geom_point() +
  
  scale_color_viridis(option = "inferno", limits = c(-4, 4), oob = squish)+
  scale_size(range = c(-7, 7), limits = c(0, 1)) +
  theme_bw() + 
  ylab("Gene") + 
  xlab("Cell Type") + 
  ggtitle("In Vivo") +
  scale_y_discrete(position = "right")+
  theme(axis.text.x = element_text(
    color = "black", size = 14, angle = -60, hjust = 0
  ),
  axis.text.y = element_text(
    color = "black",  hjust = 1)
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "left") +
  theme(panel.background = element_rect(fill = "black"))  # Add black background

plot_grid(dotplot_xrtestis,dotplot_inVivo, ncol = 2)
plot_grid(dotplot_inVivo_black, dotplot_xrtestis_black, ncol = 2)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

#pairwise comparison



all.markers$cluster <- factor(all.markers$cluster, levels = c("PGCs (early)", 
                                                              "PGCs (late)", 
                                                              "M prospermatogonia",
                                                              "M/T1 prospermaogonia", 
                                                              "T1 prospermatogonia", 
                                                              "Spermatogonia",
                                                              "Diff. spermatogonia (early)", 
                                                              "Diff. spermatogonia (late)", 
                                                              "Prelep/Lep/Zygotene"))

###

p.val.cutoff <- 0.05
log2fc.cutoff <- 0.4

###

# HS_NCG_only_NCG_subset_original -> HS_NCG_only_NCG_subset

# HS_NCG_only_NCG_subset <- HS_NCG_only_HS_subset


PGCE.vs.PGCL <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "PGCs (early)", ident.2 = "PGCs (late)", logfc.threshold = log2fc.cutoff, min.pct = 0)
PGCL.vs.M <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "PGCs (late)", ident.2 = "M prospermatogonia", logfc.threshold = log2fc.cutoff, min.pct = 0)
M.vs.MT1 <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "M prospermatogonia", ident.2 = "M/T1 prospermaogonia", logfc.threshold = log2fc.cutoff, min.pct = 0)
MT1.vs.T1 <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "M/T1 prospermaogonia", ident.2 = "T1 prospermatogonia", logfc.threshold = log2fc.cutoff, min.pct = 0)
T1.vs.Spg <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "T1 prospermatogonia", ident.2 = "Spermatogonia", logfc.threshold = log2fc.cutoff, min.pct = 0)
Spg.vs.DiffSpg <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "Spermatogonia", ident.2 = "Diff. spermatogonia (early)", logfc.threshold = log2fc.cutoff, min.pct = 0)
DiffSpg.vs.DiffSpg <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "Diff. spermatogonia (early)", ident.2 = "Diff. spermatogonia (late)", logfc.threshold = log2fc.cutoff, min.pct = 0)
DiffSpg.vs.Prelep <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "Diff. spermatogonia (late)", ident.2 = "Prelep/Lep/Zygotene", logfc.threshold = log2fc.cutoff, min.pct = 0)
# DiffSpg.vs.Prelep.all <- FindMarkers(HS_NCG_only_NCG_subset, ident.1 = "Diff. spermatogonia (late)", ident.2 = "Prelep/Lep/Zygotene", logfc.threshold = 0, min.pct = 0)


list.of.comparisons <- list(PGCE.vs.PGCL, PGCL.vs.M, M.vs.MT1, MT1.vs.T1, T1.vs.Spg, Spg.vs.DiffSpg, DiffSpg.vs.DiffSpg, DiffSpg.vs.Prelep)
list.of.counts <- list(av.counts.pgcE, av.counts.pgcL, av.counts.m, av.counts.mt1, av.counts.t1, av.counts.undiff, av.counts.diffE, av.counts.diffL, av.counts.Spermatocytes)
list.of.cell.types <- c("PGCs (early)", 
                        "PGCs (late)", 
                        "M prospermatogonia",
                        "M/T1 prospermaogonia", 
                        "T1 prospermatogonia", 
                        "Spermatogonia",
                        "Diff. spermatogonia (early)", 
                        "Diff. spermatogonia (late)", 
                        "Prelep/Lep/Zygotene")

for(number in 1:(length(list.of.cell.types)-1)){
  # number <- 1
  current.comparison <- list.of.comparisons[[number]]
  current.counts <- list.of.counts[[number]]
  current.counts2 <- list.of.counts[[number+1]]
  
  current.comparison$gene <- rownames(current.comparison)
  # current.comparison.filtered <- filter(current.comparison, p_val_adj < p.val.cutoff)
  current.comparison.filtered.up <- filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC > 0)
  current.comparison.filtered.down <- filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC < 0)
  current.comparison.dataframe <- data.frame(current.counts, current.counts2)
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  current.comparison.subset <- filter(current.comparison.dataframe, av.counts.diffL > 0 | av.counts.Spermatocytes>0)
  # current.comparison.sig.up <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.up$gene)
  # current.comparison.sig.down <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.down$gene)
  current.comparison.subset <- na.omit(current.comparison.subset)
  current.comparison.subset$DEGs <- "Not significant"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene] <- "Up"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene] <- "Down"
  print(paste0(list.of.cell.types[number], " vs ", list.of.cell.types[number+1], 
               " Up: ", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]), 
               " Down:", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])))
  
  
  assign(paste0("p", number), (ggplot() + geom_point(data = current.comparison.subset, aes(x = current.counts, y = current.counts2, color = DEGs)) + 
                                 scale_color_manual(values=c('#E69F00','#999999', '#56B4E9')) + 
                                 theme_classic() + labs(title=paste(list.of.cell.types[number], " vs. ", list.of.cell.types[number+1], ", q", p.val.cutoff, ", log2fc", log2fc.cutoff),
                                                        x=list.of.cell.types[number], y = list.of.cell.types[number+1]) + theme(plot.title = element_text(size=9))) 
  )
}

plot_grid(ncol = 4, p1, p2, p3, p4, 
          p5, p6, p7, p8)


###PSEUDOTIME STARTS HERE###

HS_NCG_only_NCG2 <- HS_NCG_only_NCG

cds <- as.cell_data_set(HS_NCG_only_NCG2)


cds <- cluster_cells(cds, reduction_method = c("UMAP"), #cluster_method = "louvain",
                     resolution = 5e-5, 
                     k = 8,
                     partition_qval = 0.05
)

plot_cells(cds, show_trajectory_graph = FALSE)
cds <- learn_graph(cds, 
                   close_loop = F,
                   learn_graph_control=list(ncenter=300, minimal_branch_len=7), #250, 18 is pretty good
                   use_partition = F)
print(plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)) #+ ggtitle(a))

# Figure 2C #######################################
cds <- order_cells(cds)#, root_cells = max.avp)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)+ NoLegend()


cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(HS_NCG_only_NCG2[["RNA"]])
cds_subset <- cds["SYCP3",]
plot_genes_in_pseudotime(cds_subset)

pseudotime.values <- pseudotime(cds, reduction_method = "UMAP")


HS_NCG_only_NCG2[["pseudotime"]] <- pseudotime.values


# c("SYCP1", "MEIOB", "MEIOC", "SPO11")
#"STRA8", "EGR4", "ESX1", "MAGEC2"

FeaturePlot(HS_NCG_only_NCG2, features = c("SYCP1", "MEIOB", "MEIOC", "SPO11"))  & scale_color_viridis(option="magma") & DarkTheme()
DimPlot(HS_NCG_only_NCG2, group.by = "experiment")
DimPlot(HS_NCG_only_NCG2, group.by = "orig.ident")
DimPlot(HS_NCG_only_NCG2, group.by = "seurat_clusters", label = T)
Idents(HS_NCG_only_NCG2) <- "seurat_clusters"
cluster.16 <- subset(HS_NCG_only_NCG2, idents = 16)
Idents(cluster.16) <- "orig.ident"
table(Idents(cluster.16))

setwd("~/Documents/scRNAseq/Rscripts/savefiles")
save(HS_NCG_only_NCG2, file = "HS_NCG_only_NCG2_2023_05_12.Robj")


### heatmap ###


genes <- c("POU5F1", "NANOG", "TCL1A", "PDPN", "PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1",
           "UTF1", "KIT", "ASB9", "NANOS2", "PIWIL4", "MORC1", "MAGEB2", "TDRD1", "RHOXF1", "MAGEC2", "ZBTB16", "EGR4", "RET", "MAGEA4",
           "SOHLH1", "DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8", "DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9", "ZCWPW1", "RAD51AP2", "SYCP3",
           "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1")

# "MLH3", "PIWIL1", "MLH1", "SPATA8", "CCDC112", "OVOL1", "OVOL2",
# "CCDC42", "ACR", "CDRT15", "TEKT1", "H1FOO", "TJP3", "FAM24A", "CATSPER3", "SPACA3", "SPACA1", "ACRV1", "SAXO1", "CAPZA3", "HEMGN", "CA2", "PRM3",
# "TNP2", "TSSK6", "SMCP", "SPEM2", "LELP1", "MKI67", "UBE2C", "TOP2A")



DefaultAssay(HS_NCG_only_NCG2) <- "RNA"


#can I do this with seurat objects? Need to pull out the normalized counts

pt.matrix.all <- HS_NCG_only_NCG2@assays$RNA@data
pseudotime <- as.data.frame(HS_NCG_only_NCG2$pseudotime)

dim(pt.matrix.all)
# gene.names.all <- rownames(pt.matrix.all)

# length(genes)
# length(intersect(genes, gene.names.all))
# setdiff(genes, gene.names.all)

pt.matrix.all[1:4, 1:4]





pt.matrix.subset <- pt.matrix.all[genes,]
t.pt.matrix.all <- t(pt.matrix.subset)
t.pt.matrix.all[1:4, 1:4]


pseudotime[1:4,1]
t.pt.matrix.all <- as.data.frame(t.pt.matrix.all)
dim(t.pt.matrix.all)
t.pt.matrix.all$pseudotime <- pseudotime[,1]
dim(t.pt.matrix.all)
t.pt.matrix.all[1:5, 1:5]
t.pt.matrix.sorted <- t.pt.matrix.all[order(t.pt.matrix.all$pseudotime),]
t.pt.matrix.sorted[1:5, ]


# pt.matrix <- t(as.matrix(t.pt.matrix.sorted[,1:(length(genes))]))
pt.matrix <- t(as.matrix(t.pt.matrix.sorted))
pt.matrix[, 1:5]

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=5)$y}))
# pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=2)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
dim(pt.matrix)
pt.matrix <- pt.matrix[c(1:length(genes)),]
dim(pt.matrix)
# rownames(pt.matrix) <- genes;
# #K means with 6 groups
# htkm <- Heatmap(
#   pt.matrix,
#   name                         = "z-score",
#   col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
#   show_row_names               = TRUE,
#   show_column_names            = FALSE,
#   row_names_gp                 = gpar(fontsize = 6),
#   clustering_distance_rows     = "kendall", #"euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
#   use_raster                   = FALSE,
#   km                           = 6, #originally 6
#   row_title_rot                = 0,
#   cluster_rows                 = FALSE,
#   cluster_row_slices           = FALSE, #false originally
#   cluster_columns              = FALSE)
# 
# print(htkm)

#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 16),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  use_raster                   = FALSE,
  cluster_rows                 = FALSE,
  cluster_row_slices           = FALSE,
  cluster_columns              = FALSE)

print(hthc)





##subset spermatids

DimPlot(HS_NCG)

HS_NCG_spermatocytes <- subset(HS_NCG, idents = c(#"Diff. spermatogonia (late)", 
  "Prelep/Lep/Zygotene", 
  "Pachytene spermatocytes",
  "Diplotene/2ndary spermatocytes",
  "Spermatids"
))


DefaultAssay(object = HS_NCG_spermatocytes) <- "integrated"
HS_NCG_spermatocytes <- ScaleData(object = HS_NCG_spermatocytes, vars.to.regress = c("nCount_RNA", "percent.mito"), verbose = FALSE)
HS_NCG_spermatocytes <- RunPCA(object = HS_NCG_spermatocytes,
                               npcs = 30,
                               verbose = TRUE)
HS_NCG_spermatocytes <- FindNeighbors(object = HS_NCG_spermatocytes)
HS_NCG_spermatocytes <- FindClusters(object = HS_NCG_spermatocytes,
                                     resolution = 10,
                                     #algorithm = 1,
                                     verbose = TRUE
)
# HS_NCG_spermatocytes <- RunTSNE(object = HS_NCG_spermatocytes, dims = 1:20)
HS_NCG_spermatocytes <- RunUMAP(object = HS_NCG_spermatocytes, reduction = "pca",
                                dims = 1:30)

DefaultAssay(HS_NCG_spermatocytes) <- "RNA"
DimPlot(HS_NCG_spermatocytes, group.by = "cell.type", label = T)
DimPlot(HS_NCG_spermatocytes)
FeaturePlot(HS_NCG_spermatocytes, features = c("STRA8", "MEIOB", "TEX101", "PIWIL1", "SPATA8", "CCDC42", "ACRV1", 
                                               "TNP1", "PRM1"))
FeaturePlot(HS_NCG_spermatocytes, features = c("DMRT1", "SOHLH2", "ESX1", "STRA8", "MEIOB", "TEX101", "PIWIL1", "SPATA8", "ACRV1"))

gene.list.diff.spermatgonia <- c(
  
  "MAGEA4",
  "SOHLH1", "DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8")

gene.list.early.meiosis <- c(
  "DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9", "ZCWPW1", "RAD51AP2", "SYCP3",
  "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1", "MLH3")
gene.list.late.meiosis <- c(
  "PIWIL1", "MLH1", "SPATA8", "CCDC112", "OVOL1", "OVOL2",
  "CCDC42", "ACR", "CDRT15", "TEKT1", "H1FOO", "TJP3", "FAM24A", "CATSPER3", "SPACA3", "SPACA1", "ACRV1", "SAXO1")
gene.list.spermatids <- c("CAPZA3", "HEMGN", "CA2", "PRM3",
                          "TNP2", "TSSK6", "SMCP", "SPEM2", "LELP1", "MKI67", "UBE2C", "TOP2A")


FeaturePlot(HS_NCG_spermatocytes, features = c("TEX19", "MLH3", "PRDM9"), split.by = "experiment")







pt.matrix.ncg <- HS_NCG_only_NCG@assays$RNA@data
pt.matrix.hs <- HS_NCG_only_HS@assays$RNA@data
pseudotime.ncg <- as.data.frame(HS_NCG_only_NCG$pseudotime)
pseudotime.hs <- as.data.frame(HS_NCG_only_HS$pseudotime)

celltype.names.ncg <- as.data.frame(HS_NCG_only_NCG$cell.type)
celltype.names.hs <- as.data.frame(HS_NCG_only_HS$cell.type)

dim(pt.matrix.ncg)

Keren.gene.list <- c("POU5F1", "NANOG", "TCL1A", "PDPN", "PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1",
                     "UTF1", "KIT", "ASB9", "NANOS2", "PIWIL4", "MORC1", "MAGEB2", "TDRD1", "RHOXF1", "MAGEC2", "ZBTB16", "EGR4", "RET", "MAGEA4",
                     "SOHLH1", "DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8", "DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9", "ZCWPW1", "RAD51AP2", "SYCP3",
                     "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1", "MLH3", "PIWIL1", "MLH1", "SPATA8", "CCDC112", "OVOL1", "OVOL2",
                     "CCDC42", "ACR", "CDRT15", "TEKT1", "H1FOO", "TJP3", "FAM24A", "CATSPER3", "SPACA3", "SPACA1", "ACRV1", "SAXO1", "CAPZA3", "HEMGN", "CA2", "PRM3",
                     "TNP2", "TSSK6", "SMCP", "SPEM2", "LELP1", "MKI67", "UBE2C", "TOP2A")


genes <- Keren.gene.list

genes <- c("POU5F1", "NANOS2", "MAGEB2", "SOHLH1", "STRA8", "DAZL", "SYCP1", "ACRV1", 
           "TNP2")



### FOR A SUBSET
# pt.matrix.subset.ncg <- data.frame(pt.matrix.ncg[genes,])
# pt.matrix.subset.hs <- data.frame(pt.matrix.hs[genes,])

###FOR ALL GENES
pt.matrix.subset.ncg <- as.matrix(pt.matrix.ncg)
pt.matrix.subset.hs <- as.matrix(pt.matrix.hs)


t.pt.matrix.all.ncg <- t(pt.matrix.subset.ncg)
t.pt.matrix.all.ncg[1:2, 1:2]
t.pt.matrix.all.ncg <- as.data.frame(t.pt.matrix.all.ncg)
dim(t.pt.matrix.all.ncg)
t.pt.matrix.all.ncg$pseudotime <- pseudotime.ncg[,1]
t.pt.matrix.all.ncg$cell.type <- celltype.names.ncg[,1]
dim(t.pt.matrix.all.ncg)
# t.pt.matrix.all.ncg[1:5, 1:5]
t.pt.matrix.sorted.ncg <- t.pt.matrix.all.ncg[order(t.pt.matrix.all.ncg$pseudotime),]
t.pt.matrix.sorted.ncg[1:5, ]

t.pt.matrix.all.hs <- t(pt.matrix.subset.hs)
t.pt.matrix.all.hs[1:2, 1:2]
t.pt.matrix.all.hs <- as.data.frame(t.pt.matrix.all.hs)
dim(t.pt.matrix.all.hs)
t.pt.matrix.all.hs$pseudotime <- pseudotime.hs[,1]
t.pt.matrix.all.hs$cell.type <- celltype.names.hs[,1]
dim(t.pt.matrix.all.hs)
# t.pt.matrix.all.hs[1:5, 1:5]
t.pt.matrix.sorted.hs <- t.pt.matrix.all.hs[order(t.pt.matrix.all.hs$pseudotime),]
t.pt.matrix.sorted.hs[1:5, 1:2]




for(current.loop in 1:length(genes)){
  # current.loop <- 1
  gene.of.interest <- genes[current.loop]
  
  gene.max.hs <- max(t.pt.matrix.sorted.hs[,gene.of.interest])
  gene.max.ncg <- max(t.pt.matrix.sorted.ncg[,gene.of.interest])
  gene.max <- max(gene.max.hs, gene.max.ncg)
  
  pseudotime.max.hs <- max(pseudotime.hs)
  pseudotime.max.ncg <- max(pseudotime.ncg)
  pseudotime.max <- max(pseudotime.max.hs, pseudotime.max.ncg)
  
  a1 <- ggplot(t.pt.matrix.sorted.hs, aes(x=pseudotime, y=.data[[gene.of.interest]], color = cell.type)) +  #the .data is needed to pull the variable name
    geom_point()+ theme(legend.position = "none") + xlim(0, pseudotime.max) + ylim(0, gene.max) +ggtitle("In vivo")
  a2 <- ggplot(t.pt.matrix.sorted.ncg, aes(x=pseudotime, y=.data[[gene.of.interest]], color = cell.type)) + 
    geom_point()+ theme(legend.position = "none")+ xlim(0, pseudotime.max) + ylim(0, gene.max)+ggtitle("In vitro")
  
  q1 <- plot_grid(a1, a2, ncol = 1)
  assign(paste0("p", current.loop), q1)
  
  
}

plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, nrow = 3)

ggplot(t.pt.matrix.sorted.hs, aes(x=pseudotime, y=.data[[gene.of.interest]], color = cell.type)) +  #the .data is needed to pull the variable name
  geom_point()+ xlim(0, pseudotime.max) + ylim(0, gene.max) +ggtitle("In vivo")

# p1 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=POU5F1, color = cell.type)) + geom_point() + theme(legend.position = "none")
# p2 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=NANOG, color = cell.type)) + geom_point() + theme(legend.position = "none")
# p3 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=UTF1, color = cell.type)) + geom_point() + theme(legend.position = "none")
# p4 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=KIT, color = cell.type)) + geom_point() + theme(legend.position = "none")
# p5 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=NANOS2, color = cell.type)) + geom_point() + theme(legend.position = "none")
# p6 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=SOHLH1, color = cell.type)) + geom_point() + theme(legend.position = "none")
# p7 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=MAGEB2, color = cell.type)) + geom_point() + theme(legend.position = "none")
# p8 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=DMRT1, color = cell.type)) + geom_point() + theme(legend.position = "none")
# p9 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=STRA8, color = cell.type)) + geom_point() + theme(legend.position = "none")
# p10 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=DAZL, color = cell.type)) + geom_point() + theme(legend.position = "none")
# p11 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=PRDM9, color = cell.type)) + geom_point() + theme(legend.position = "none")
# p12 <- ggplot(t.pt.matrix.sorted, aes(x=pseudotime, y=SYCP1, color = cell.type)) + geom_point() + theme(legend.position = "none")

plot_grid(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, p11, p12, nrow = 4)



Idents(HS_NCG) <- "experiment"
DimPlot(HS_NCG)
DimPlot(HS_NCG, group.by = "cell.type.combined", split.by = "experiment")

HS_NCG_only_HS <- subset(HS_NCG, idents = "human in vivo")
# HS_NCG_only_HS <- subset(HS_NCG, idents = "human_in_vivo")
HS_NCG_only_NCG <- subset(HS_NCG, idents = c("transplant", "hPGCLC"))
# HS_NCG_only_NCG <- subset(HS_NCG, idents = c("transplant_human", "hPGCLC"))


DimPlot(HS_NCG)
DimPlot(HS_NCG_only_NCG)
DimPlot(HS_NCG_only_HS)



Idents(HS_NCG) <- "cell.type.combined"

HS_NCG_PGCE <- subset(HS_NCG, idents = "PGCs (early)")
HS_NCG_PGCL<- subset(HS_NCG, idents = "PGCs (late)")
HS_NCG_M<- subset(HS_NCG, idents = "M prospermatogonia")
HS_NCG_MT1<- subset(HS_NCG, idents = "M/T1 prospermaogonia") #NOTE SPELLING ERROR IS STILL HERE
HS_NCG_T<- subset(HS_NCG, idents = "T1 prospermatogonia")
HS_NCG_spg<- subset(HS_NCG, idents = "Spermatogonia")
HS_NCG_diffspgE<- subset(HS_NCG, idents = "Diff. spermatogonia (early)")
HS_NCG_diffspgL<- subset(HS_NCG, idents = "Diff. spermatogonia (late)")
# HS_NCG_prelep<- subset(HS_NCG, idents = "Preleptotene spermatocytes")
# HS_NCG_lepzyg<- subset(HS_NCG, idents = "Leptotene/Zygotene")
HS_NCG_prelep <- subset(HS_NCG, idents = "Prelep/Lep/Zygotene")

DimPlot(HS_NCG, group.by = "experiment")
Idents(HS_NCG_PGCE) <- "experiment"
Idents(HS_NCG_PGCL) <- "experiment"
Idents(HS_NCG_M) <- "experiment"
Idents(HS_NCG_MT1) <- "experiment"
Idents(HS_NCG_T) <- "experiment"
Idents(HS_NCG_spg) <- "experiment"
Idents(HS_NCG_diffspgE) <- "experiment"
Idents(HS_NCG_diffspgL) <- "experiment"
Idents(HS_NCG_prelep) <- "experiment"
# Idents(HS_NCG_lepzyg) <- "experiment"

PGCs.E.markers.hs.vs.ncg <- FindMarkers(HS_NCG_PGCE, ident.1 = "human in vivo", ident.2 = "hPGCLC")
PGCs.L.markers.hs.vs.ncg <- FindMarkers(HS_NCG_PGCL, ident.1 = "human in vivo", ident.2 = "transplant")
M.ProSpermatogonia.markers.hs.vs.ncg <- FindMarkers(HS_NCG_M, ident.1 = "human in vivo", ident.2 = "transplant")
M2T1.markers.hs.vs.ncg <- FindMarkers(HS_NCG_MT1, ident.1 = "human in vivo", ident.2 = "transplant")
T.ProSpermatogonia.markers.hs.vs.ncg <- FindMarkers(HS_NCG_T, ident.1 = "human in vivo", ident.2 = "transplant")
Spg.markers.hs.vs.ncg <- FindMarkers(HS_NCG_spg, ident.1 = "human in vivo", ident.2 = "transplant")
Diff.Spg.E.markers.hs.vs.ncg <- FindMarkers(HS_NCG_diffspgE, ident.1 = "human in vivo", ident.2 = "transplant")
Diff.Spg.L.markers.hs.vs.ncg <- FindMarkers(HS_NCG_diffspgL, ident.1 = "human in vivo", ident.2 = "transplant")
Pre.Leptotene.markers.hs.vs.ncg <- FindMarkers(HS_NCG_prelep, ident.1 = "human in vivo", ident.2 = "transplant")
# LepZyg.markers.hs.vs.ncg <- FindMarkers(HS_NCG_lepzyg, ident.1 = "human in vivo", ident.2 = "transplant")
# 
# setwd("/Users/ewhelan/Desktop/NCG/NCG_marker_genes2")
# write.csv(PGCs.E.markers.hs.vs.ncg, file="HS_NCG_PGCs.E.markers.hs.vs.ncg.csv")
# write.csv(PGCs.L.markers.hs.vs.ncg, file="HS_NCG_PGCs.L.markers.hs.vs.ncg.csv")
# write.csv(M.ProSpermatogonia.markers.hs.vs.ncg, file="HS_NCG_M.ProSpermatogonia.markers.hs.vs.ncg.csv")
# write.csv(M2T1.markers.hs.vs.ncg, file="HS_NCG_M2T1.markers.hs.vs.ncg.csv")
# write.csv(T.ProSpermatogonia.markers.hs.vs.ncg, file="HS_NCG_T.ProSpermatogonia.markers.hs.vs.ncg.csv")
# write.csv(Spg.markers.hs.vs.ncg, file="HS_NCG_Spg.markers.hs.vs.ncg.csv")
# write.csv(Diff.Spg.E.markers.hs.vs.ncg, file="HS_NCG_Diff.Spg.E.markers.hs.vs.ncg.csv")
# write.csv(Diff.Spg.L.markers.hs.vs.ncg, file="HS_NCG_Diff.Spg.L.markers.hs.vs.ncg.csv")
# write.csv(Pre.Leptotene.markers.hs.vs.ncg, file="HS_NCG_Pre.Leptotene.markers.hs.vs.ncg.csv")

FeaturePlot(HS_NCG, features = "KRT18", split.by = "experiment")

DimPlot(HS_NCG)
DimPlot(HS_NCG_only_NCG)

DimPlot(HS_NCG_only_HS)
Idents(HS_NCG_only_HS) <- "cell.type"
Idents(HS_NCG_only_NCG) <- "cell.type"

data.PGCE.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "PGCs (early)")])
data.PGCL.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "PGCs (late)")])
data.M.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "M prospermatogonia")])
data.MT1.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "M/T1 propermaogonia")])
data.T1.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "T1 prospermatogonia")])
data.undiff.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "Spermatogonia")])
data.diffE.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "Diff. spermatogonia (early)")])
data.diffL.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "Diff. spermatogonia (late)")])
data.prelep.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "Preleptotene spermatocytes")])
# data.lepzyg.hs <- as.matrix(GetAssayData(HS_NCG_only_HS, slot = "data")[, WhichCells(HS_NCG_only_HS, ident = "Leptotene/Zygotene")])

data.PGCE.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "PGCs (early)")])
data.PGCL.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "PGCs (late)")])
data.M.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "M prospermatogonia")])
data.MT1.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "M/T1 propermaogonia")])
data.T1.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "T1 prospermatogonia")])
data.undiff.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "Spermatogonia")])
data.diffE.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "Diff. spermatogonia (early)")])
data.diffL.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "Diff. spermatogonia (late)")])
data.prelep.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "Preleptotene spermatocytes")])
# data.lepzyg.ncg <- as.matrix(GetAssayData(HS_NCG_only_NCG, slot = "data")[, WhichCells(HS_NCG_only_NCG, ident = "Leptotene/Zygotene")])

av.counts.PGCE.hs <- apply(data.PGCE.hs, 1, mean)
av.counts.PGCL.hs <- apply(data.PGCL.hs, 1, mean)
av.counts.M.hs <- apply(data.M.hs, 1, mean)
av.counts.MT1.hs <- apply(data.MT1.hs, 1, mean)
av.counts.T1.hs <- apply(data.T1.hs, 1, mean)
av.counts.undiff.hs <- apply(data.undiff.hs, 1, mean)
av.counts.diffE.hs <- apply(data.diffE.hs, 1, mean)
av.counts.diffL.hs <- apply(data.diffL.hs, 1, mean)
av.counts.prelep.hs <- apply(data.prelep.hs, 1, mean)

av.counts.PGCE.ncg <- apply(data.PGCE.ncg, 1, mean)
av.counts.PGCL.ncg <- apply(data.PGCL.ncg, 1, mean)
av.counts.M.ncg <- apply(data.M.ncg, 1, mean)
av.counts.MT1.ncg <- apply(data.MT1.ncg, 1, mean)
av.counts.T1.ncg <- apply(data.T1.ncg, 1, mean)
av.counts.undiff.ncg <- apply(data.undiff.ncg, 1, mean)
av.counts.diffE.ncg <- apply(data.diffE.ncg, 1, mean)
av.counts.diffL.ncg <- apply(data.diffL.ncg, 1, mean)
av.counts.prelep.ncg <- apply(data.prelep.ncg, 1, mean)

av.counts.hs <- cbind(av.counts.PGCE.hs, av.counts.PGCL.hs, av.counts.M.hs, av.counts.MT1.hs,
                      av.counts.T1.hs, av.counts.undiff.hs, av.counts.diffE.hs, av.counts.diffL.hs, 
                      av.counts.prelep.hs)

av.counts.ncg <- cbind(av.counts.PGCE.ncg, av.counts.PGCL.ncg, av.counts.M.ncg, av.counts.MT1.ncg,
                       av.counts.T1.ncg, av.counts.undiff.ncg, av.counts.diffE.ncg, av.counts.diffL.ncg, 
                       av.counts.prelep.ncg)


av.counts.hs.df <- as.data.frame(av.counts.hs)
av.counts.ncg.df <- as.data.frame(av.counts.ncg)

gene.name.list.hs <- rownames(av.counts.hs.df)
gene.name.list.ncg <- rownames(av.counts.ncg.df)

length(gene.name.list.hs)
length(gene.name.list.ncg)
length(intersect(gene.name.list.hs, gene.name.list.ncg))

###

#fix the cell names

levels(HS_NCG$cell.type.combined)
HS_NCG$cell.type.combined.short <- HS_NCG$cell.type.combined
DimPlot(HS_NCG, group.by = "cell.type.combined.short", split.by = "experiment")
Idents(HS_NCG) <- "cell.type.combined"
HS_NCG_subset <- subset(HS_NCG, idents = c("PGCs (early)", "Pachytene spermatocytes", "Diplotene/2ndary spermatocytes", "Spermatids"), invert = T)
DimPlot(HS_NCG_subset)
levels(HS_NCG_subset$cell.type.combined.short)

levels(HS_NCG_subset$cell.type.combined.short)[levels(HS_NCG_subset$cell.type.combined.short) == "PGCs (early)"] <- "PGCsE"
levels(HS_NCG_subset$cell.type.combined.short)[levels(HS_NCG_subset$cell.type.combined.short) == "PGCs (late)"] <- "PGCsL"
levels(HS_NCG_subset$cell.type.combined.short)[levels(HS_NCG_subset$cell.type.combined.short) == "M prospermatogonia"] <- "M"
levels(HS_NCG_subset$cell.type.combined.short)[levels(HS_NCG_subset$cell.type.combined.short) == "M/T1 prospermaogonia"] <- "M2T1"
levels(HS_NCG_subset$cell.type.combined.short)[levels(HS_NCG_subset$cell.type.combined.short) == "T1 prospermatogonia"] <- "T1"
levels(HS_NCG_subset$cell.type.combined.short)[levels(HS_NCG_subset$cell.type.combined.short) == "Spermatogonia"] <- "Spg"
levels(HS_NCG_subset$cell.type.combined.short)[levels(HS_NCG_subset$cell.type.combined.short) == "Diff. spermatogonia (early)"] <- "DiffSpgE"
levels(HS_NCG_subset$cell.type.combined.short)[levels(HS_NCG_subset$cell.type.combined.short) == "Diff. spermatogonia (late)"] <- "DiffSpgL"
levels(HS_NCG_subset$cell.type.combined.short)[levels(HS_NCG_subset$cell.type.combined.short) == "Prelep/Lep/Zygotene"] <- "Prelep"

levels(HS_NCG_subset$cell.type.combined.short)
Idents(HS_NCG_subset) <- HS_NCG_subset$cell.type.combined.short
DimPlot(HS_NCG_subset)

Idents(HS_NCG_subset) <- HS_NCG_subset$experiment
HS_NCG_subset$experiment <- factor(x = HS_NCG_subset$experiment, levels = c("human in vivo", "transplant"))

levels(HS_NCG_subset$experiment)
levels(HS_NCG_subset$experiment)[levels(HS_NCG_subset$experiment) == "human in vivo"] <- "inVivo"
levels(HS_NCG_subset$experiment)[levels(HS_NCG_subset$experiment) == "transplant"] <- "xrTestis"


#let's try with the ?AverageExpression function
Idents(HS_NCG_subset) <- HS_NCG_subset$experiment
DimPlot(HS_NCG_subset)

DimPlot(HS_NCG_subset, group.by = "cell.type.combined.short", split.by = "experiment")
average.expression.output <- AverageExpression(HS_NCG_subset, group.by = c("cell.type.combined.short", "experiment"))$RNA
average.expression.output.mx <- as.matrix(average.expression.output)

# #aggregateexpression is recommended
# aggregate.expression.output <- AggregateExpression(object = HS_NCG_noPGCLCs, group.by = c("cell.type.combined", "experiment"))$RNA
# 
# aggregate.expression.output.mx <- as.matrix(aggregate.expression.output)

###

# Select the columns you want to compare
hs_column <- av.counts.hs.df[, -ncol(av.counts.hs.df)]  # Exclude the "gene" column
ncg_column <- av.counts.ncg.df[, -ncol(av.counts.ncg.df)]  # Exclude the "gene" column

# Compute correlation between corresponding columns
correlation_values <- cor(hs_column, ncg_column, use = "pairwise.complete.obs")

# Show correlation values
print(correlation_values)
rsquared <- correlation_values^2
print(rsquared)

inVivoData <- av.counts.hs.df
inVitroData <- av.counts.ncg.df
correlation_matrix <- cor(inVivoData, inVitroData)
print(correlation_matrix)
correltation_matrix_r <- correlation_matrix
correlation_matrix <- correlation_matrix^2



# Convert the correlation matrix to a tibble
correlation_df <- as_tibble(as.data.frame(rsquared), rownames = "Cell_Type_1")

# Factor the levels to maintain the original order
correlation_df$Cell_Type_1 <- factor(correlation_df$Cell_Type_1, levels = rownames(correlation_matrix))
correlation_df <- gather(correlation_df, key = "Cell_Type_2", value = "Correlation", -Cell_Type_1)

# Factor the levels for Cell_Type_2
correlation_df$Cell_Type_2 <- factor(correlation_df$Cell_Type_2, levels = colnames(correlation_matrix))

# Create the heatmap using ggplot2
heatmap_plot <- ggplot(correlation_df, aes(x = Cell_Type_1, y = Cell_Type_2, fill = Correlation)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = round(Correlation, 2)), vjust = 1) +  # Add correlation values
  # scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.8, na.value = NA) +
  # scale_fill_gradient2(low = "#451d5b", mid = "white", high = "red", midpoint = 0.8, na.value = NA) +
  scale_fill_viridis() +
  theme_minimal() +
  labs(title = "Correlation Heatmap",
       x = "Cell Type (In Vivo)",
       y = "Cell Type (In Vitro)") +
  theme(axis.text.x = element_text(angle = -22, hjust = 1))  # Rotate x-axis labels at a 45-degree angle

# Print the heatmap
print(heatmap_plot)


###

av.counts.hs.df$gene <- gene.name.list.hs
av.counts.ncg.df$gene <- gene.name.list.ncg


list.of.comparisons <- list(PGCs.E.markers.hs.vs.ncg, 
                            PGCs.L.markers.hs.vs.ncg, 
                            M.ProSpermatogonia.markers.hs.vs.ncg, M2T1.markers.hs.vs.ncg, T.ProSpermatogonia.markers.hs.vs.ncg, 
                            Spg.markers.hs.vs.ncg, Diff.Spg.E.markers.hs.vs.ncg, Diff.Spg.L.markers.hs.vs.ncg, 
                            Pre.Leptotene.markers.hs.vs.ncg)
# setwd("~/Desktop/NCG/NCG_vs_HS")
# write.csv(PGCs.E.markers.hs.vs.ncg, file="PGCs.E.markers.hs.vs.ncg.csv")
# write.csv(PGCs.L.markers.hs.vs.ncg, file="PGCs.L.markers.hs.vs.ncg.csv")
# write.csv(M.ProSpermatogonia.markers.hs.vs.ncg, file="M.ProSpermatogonia.markers.hs.vs.ncg.csv")
# write.csv(M2T1.markers.hs.vs.ncg, file="M2T1.markers.hs.vs.ncg.csv")
# write.csv(T.ProSpermatogonia.markers.hs.vs.ncg, file="T.ProSpermatogonia.markers.hs.vs.ncg.csv")
# write.csv(Spg.markers.hs.vs.ncg, file="Spg.markers.hs.vs.ncg.csv")
# write.csv(Diff.Spg.E.markers.hs.vs.ncg, file="Diff.Spg.E.markers.hs.vs.ncg.csv")
# write.csv(Diff.Spg.L.markers.hs.vs.ncg, file="Diff.Spg.L.markers.hs.vs.ncg.csv")
# write.csv(Pre.Leptotene.markers.hs.vs.ncg, file="Pre.Leptotene.markers.hs.vs.ncg.csv")


# Select the columns you want to compare
hs_column <- av.counts.hs.df[, -ncol(av.counts.hs.df)]  # Exclude the "gene" column
ncg_column <- av.counts.ncg.df[, -ncol(av.counts.ncg.df)]  # Exclude the "gene" column

# Compute correlation between corresponding columns
correlation_values <- cor(hs_column, ncg_column, use = "pairwise.complete.obs")

# Show correlation values
print(correlation_values)
rsquared <- correlation_values^2
print(rsquared)

list.of.cell.types <- c("PGCs (early)", 
                        "PGCs (late)",
                        "M prospermatogonia",
                        "M/T1 propermaogonia", #NOTE SPELLING ERROR IS STILL HERE
                        "T1 prospermatogonia",
                        "Spermatogonia",
                        "Diff. spermatogonia (early)",
                        "Diff. spermatogonia (late)",
                        "Preleptotene spermatocytes")

list.of.cell.types.short <- c("PGCsE",
                              "PGCsL",
                              "M",
                              "M2T1", 
                              "T1",
                              "Spg",
                              "DiffSpgE",
                              "DiffSpgL",
                              "Prelep")

# "Leptotene/Zygotene") #not including 

p.val.cutoff <- 0.05
# log2fc.cutoff <- 0.585
log2fc.cutoff <- 0.3

viridis.colours <- inferno(5)

### Extended Data Fig 8 ###

genes.of.interest <- c("SMC1B", "VCX", "VCX2", "SMC3", "TEX101", "SYCE",
                       "IL13RA2", "MDK", "TUBB", "RPS5", "ID4", "RB1")

gene.of.interest.PGCE <- c("TFAP2C", "SOX17", "YBX1", "MIF", "POU5F1", 
                           "NANOS3", "TOP2A", "KHDC3L", "MKI67", "DND1", "DAZL")
gene.of.interest.PGCL1 <- c("DDX4", "DAZL", "MKI67",  "TOP2A", "ID4")
gene.of.interest.PGCL2 <- c("ACTG2",  "POU5F1", "S100A13", "CCNB1", "PCSK1N" )
gene.of.interest.M1 <- c("DAZL", "ID4", "DDX43", "DDX4", "E2F8")
gene.of.interest.M2 <- c("NANOG", "ACTG2", "MKI67", "NANOS3", "POU5F1", "STMN1", "CENPA")
gene.of.interest.transitional1 <- c("MDK", "ASB9", "CITED2", "MAP4K4", "TDRD1", "ID3")
gene.of.interest.transitional2 <- c("MDK", "ASB9", "NANOS3", "POU5F1", "TSPAN6", "APOE", "STMN1", "MEG3")
gene.of.interest.T1 <- c("NANOS2", "RHOXF1", "TDRD1", "MAGEB2", "TEX15", "DCAF4L1", "MAGEC2", "PAGE2B",
                         "SIX1", "DDX4", "FOXP1", "TCF3", "VCX", 'EGR4', "VCX3B", "DPPA2")
gene.of.interest.T1_2 <- c("MEG3", "NANOS2", "TPM2")
gene.of.interest.spg1 <- c("PIWIL4", "SIX1", "HES5", "MAGEB2",  "TCF3", "UTF1", "EGR4")
gene.of.interest.spg2 <- c("DPPA2", "SIX1", "FGFR3", "PIWIL4", "EGR4", "ID4", "SOHLH1", "UTF1")
gene.of.interest.earlydiff1 <- c("ASB9", "STRA8","PRSS21", "SYCP3", "MLEC", "PARP4", "REC8", "LARP1", "DNMT1", "PIWIL1", "TOP2B")
gene.of.interest.earlydiff2 <- c("KIT", "PIWIL1", "ASB9", "ETV7", "L1TD1", "SCML1", "MKI67", "TEX15", "HORMAD1", "MEIOC", "SYCP3")
gene.of.interest.latediff1 <- c("SYCP3", "SMC1B", "MEIOC", "MAGEA4", "UBE2C", "ESX1", "MKI67")
gene.of.interest.latediff2 <- c("CENPF", "KIT", "ITGA6", "UBE2C", "ELAVL2",  "RHOXF1")
gene.of.interest.prelep1 <- c("TEX101", "MEIOB", "HORMAD1", "SYCP2", "SYCP1", "ZCWPW1","SCML1", "TOP2A","SYCP3", "SPO11", "PRDM9")


genes.of.interest.PGCE <- c("TFAP2C", "SOX17", "YBX1", "POU5F1", "PRDM1", 
                            "KIT",   "KRT8", "KRT18",  
                            "PCSK1N", #up 
                            "CENPF", "SCG5", 
                            "TOP2A", "DPPA4", "NANOS3", "TOP2A", "KHDC3L",  "DND1") #down
genes.of.interest.PGCL <- c("TFAP2C", "SOX17", "GSTP1", "MIF", "POU5F1", "DDX4", "DAZL", "MKI67",  "TOP2A",
                            "KRT18",  
                            "VIM", "IFI16", "PDPN", "TMED9", "CHEK1", "CHCHD2", 
                            "NANOS3", "TOP2A", "KHDC3L", "MKI67", "DND1", "DAZL") #down

# genes.of.interest.PGCL <- c("TFAP2C", "SOX17", "YBX1", "MIF", "POU5F1", "DDX4", "DAZL", "MKI67",  "TOP2A", "ID4")
genes.of.interest.M <- c("POU5F1", "SOX17", "TFAP2C", "PLCG2", "CRABP1",
                         "DDX4", "DAZL", "MKI67",  "TOP2A", "ID4", "NANOG", "ACTG2", "MKI67", "NANOS3", "PRDX5", 
                         "SIVA1", "MT-ND1", "MT-ND2", "XIST", "PHLDA3", "TMED9", "PDPN")
genes.of.interest.intermediate <- c("PRXL2A", "TLE5", "POU5F1", "SOX17", "TFAP2C", 
                                    "RHOXF1", "PAGEB2", "NANOG", "ID1", "DCAF4L1", "MAGEC2", 
                                    "MKI67", "THY1", "CHEK1", "HEBP2", "SKAP2", 
                                    "CITED2", "MAP4K4", "TDRD1")
genes.of.interest.T1 <- c("NANOS2", "DCAF4L1", "MAGEC2", "PAGE2B", "MT-ND4L", "PRAME", 
                          "DDX4",  "TCF3", "VCX", 'EGR4', "PLPPR3",   "GAGE12H", "APOE", "RPL23", 
                          "XIST",  "ID2", "ID3")
genes.of.interest.spermatogonia.MAIN <- c("FGFR3", "SIX1", "MAGEB2", 
                                          "ID2", "ID3", "ETV7", "ACTB", "ID4", "DDX4", "MORC1", "PIWIL4")

genes.of.interest.spermatogonia <- c("LY6K", "GAGE12H", "JUND", "PRDX5", "ZNF518A",  "DHRS2",
                                     # "DPPA2", , "CTAG2"
                                     # "SIX1", "FGFR3", 
                                     "EGR4",  "SOHLH1", "UTF1", 
                                     "MT-ND5", 
                                     # "MAGEB2", 
                                     "TEX15")

genes.of.interest.earlydiff <- c("KIT","DAZL", "PAGE5", "PAGE2", "ID1", "JUND", "CYP26A1", 
                                 # "CYP26B1", 
                                 "RPL41", "RPL10A", "ZNF428", 
                                 "ASB9", "SYCP3", "SOHLH2", "L1TD1", "MKI67"
)
genes.of.interest.earlydiff.MAIN <- c("STRA8", "PIWIL1", "REC8", "TEX11", "ACTB", 
                                      "SOHLH1", "UTF1", "ID4", "MAGEA4", "DDX4")


genes.of.interest.latediff <- c("SYCP3", "SMC1B", "MEIOC", "MAGEA4", "UBE2C", "MKI67", "TOP2A", 
                                "VCX", "RPL10", "LDHB", "RACK1", "TUBB", "MDK",  "PCNA", 
                                "CENPF", "KIT", "UBE2C", "ELAVL2")
genes.of.interest.prelep <- c("TEX101", "MEIOB", "HORMAD1", "SYCP2", "SYCP1", "ZCWPW1","SCML1", "TOP2A","SYCP3", "SPO11", "PRDM9", 
                              "MEIOC", "TEX19", "RAD51AP2",  
                              "RPS12", "RPL10", "RPS5", "MDK", "TUBB")

gene.list.premigratory <- c("PCDH1", "ID1", "MSX2", "KAL1", "SMAD7", "KRT19",
                            "CYFIP2", "SLC16A2", "EPHB1", "ASS1", "PCSK1N", 
                            "CABLES1", "ADAM19", "WFDC2", "SH3BGRL3", "C1orf85",
                            "LGALS3", "TRO", "FAM54A", "CDT1", "POLA1", "GDF3", "PMEL", "GREM1", 
                            "RIN2", "PPP1R14C", "ENPP2", "EPAS1", "WSCD1", "MYH2", "HAPLN3", 
                            "TPH2", "LONRF3", "GLDC", "SCG5", "SLC7A3", "ATP6V1H", 
                            "TFRAF3IP2", "CTSL2", "T", "SC5DL", "GSN")

# genes.of.interest.T1 <- c("GAGE12H", "APOE", "RPL23", "RPS26")
# genes.of.interest.spermatogonia1 <- c("MT-ND5", "ID3", "ID1", "ID2", "MAGEB2", "TEX15", "CTAG2")
# genes.of.interest.spermatogonia2 <- c("ID2", "ID3", "ETV7",  "ID4", 
#                        "MORC1", "PIWIL4", "ACTB", "DDX4")
# gene.of.interest.earlydiff <- c("REC8","TEX11", "ACTB", "PIWIL1", "STRA8",
#                                 "SOHLH1", "UTF1", "ID4", "MAGEA4", "DDX4")
list.of.lists.of.genes <- list(genes.of.interest.PGCE,
                               genes.of.interest.PGCL,
                               genes.of.interest.M,
                               genes.of.interest.intermediate,
                               genes.of.interest.T1,
                               genes.of.interest.spermatogonia,
                               genes.of.interest.earlydiff,
                               genes.of.interest.latediff,
                               genes.of.interest.prelep
)

length(list.of.comparisons)

for(number in 1:(length(list.of.comparisons))){
  # number <- 5
  print(number)
  # gene.of.interest <- c(list.of.lists.of.genes[[number]])
  # gene.of.interest <- c("FGFR2", "FGFR3")
  gene.of.interest <- list.of.lists.of.genes[[number]]
  # gene.of.interest <- gene.list.premigratory
  current.comparison <- list.of.comparisons[[number]]
  
  # current.counts <- av.counts.hs.df[[number]]
  # current.counts2 <- av.counts.ncg.df[[number]]
  
  current.counts <- av.counts.hs.df[,number]
  current.counts2 <- av.counts.ncg.df[,number]
  
  current.comparison$gene <- rownames(current.comparison)
  # current.comparison.filtered <- filter(current.comparison, p_val_adj < p.val.cutoff)
  current.comparison.filtered.up <- filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC > 0)
  current.comparison.filtered.down <- filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC < 0)
  
  all.significant.genes <- filter(current.comparison, p_val_adj < p.val.cutoff)
  all.significant.genes <- all.significant.genes$gene
  
  current.comparison.dataframe <- data.frame(current.counts, current.counts2)
  current.comparison.dataframe$log2fc <- current.counts - current.counts2
  current.comparison.dataframe$regulation <- "notDifferent"
  # current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  current.comparison.dataframe$gene.name <- rownames(av.counts.hs.df)
  rownames(current.comparison.dataframe) <- rownames(av.counts.hs.df)
  
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc > log2fc.cutoff] <- "higherInVivo"
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc < -log2fc.cutoff] <- "higherInVitro"
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  
  all.genes <- rownames(current.comparison)
  nonsignificant.genes <- all.genes[!(all.genes %in% all.significant.genes)]
  
  
  # current.comparison.dataframe$regulation[!(current.comparison.dataframe$gene.name %in% all.significant.genes)] <- "notDifferent"
  
  
  # current.comparison.subset <- filter(current.comparison.dataframe, av.counts.diffL > 0 | av.counts.Spermatocytes>0)
  current.comparison.subset <- dplyr::filter(current.comparison.dataframe, current.counts > 0 | current.counts2 >0)
  # current.comparison.sig.up <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.up$gene)
  # current.comparison.sig.down <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.down$gene)
  current.comparison.subset <- na.omit(current.comparison.subset)
  current.comparison.subset$DEGs <- "Not significant"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene] <- "higherInVivo"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene] <- "higherInVitro"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% gene.of.interest] <- "GeneOfInterest"
  current.comparison.subset$regulation[current.comparison.subset$gene.name %in% gene.of.interest] <- "GeneOfInterest"
  print(paste0(list.of.cell.types[number], ": xrTestis vs in Vivo", 
               " Up: ", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]), 
               " Down:", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])))
  
  colnames(current.comparison.subset) <- c("HS.counts", "NCG.counts", "log2FC_manual", "regulation", "gene.name", "DEGs")
  
  print(list.of.cell.types[number])
  print(table(current.comparison.subset$regulation))
  
  
  degs_order <- c( "higherInVitro", "GeneOfInterest", "Not significant", "higherInVivo")
  
  
  current.comparison.subset_ordered <- current.comparison.subset %>%
    mutate(DEGs = factor(DEGs, levels = degs_order)) 
  
  # setwd("~/Desktop/xrTestisFigures/inVivo.vs.xrTestis")
  write.csv(current.comparison.subset_ordered, file = paste0(list.of.cell.types.short[number], "_comparison_UPinVivo_DOWNxrTestis_",
                                                             log2fc.cutoff, "_log2FC.csv"))
  table(current.comparison.subset_ordered$regulation)
  
  assign(paste0("r", number), 
         ggplot(data = current.comparison.subset_ordered, aes(x = HS.counts, y = NCG.counts,
                                                              label = gene.name,
                                                              color = regulation)) +
           rasterise(geom_point(data = filter(current.comparison.subset_ordered, regulation == "notDifferent"), alpha = 0.5), dpi = 900) + # Plot "Not Significant" points first
           # geom_point(data = filter(current.comparison.subset_ordered, regulation == "notDifferent"), alpha = 0.5) + # Plot "Not Significant" points first
           geom_point(data = filter(current.comparison.subset_ordered, regulation != "notDifferent")) + # Plot other points
           # geom_text(data = filter(current.comparison.subset_ordered, regulation != "notDifferent"),hjust=-0.1, vjust=-0.1) +
           geom_point(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest")) +
           geom_text(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest"), hjust = -0.1, vjust = -0.1) +
           
           
           # scale_color_manual(values = c(cell.type.colors[number+1], "red4", '#999999', cell.type.colors[number])) +
           scale_color_manual(values = c("black","blue",   "red3",'#999999')) +
           theme_classic() +
           xlim(0, 3.5) + ylim(0, 3.5) + 
           labs(title = paste(list.of.cell.types[number], " xrTestis vs in Vivo", ", q", p.val.cutoff, ", log2fc", log2fc.cutoff),
                x = paste0(list.of.cell.types[number], " in Vivo"), y = paste0(list.of.cell.types[number], " xrTestis")) +
           theme(plot.title = element_text(size = 9),
                 panel.grid.major = element_line(color = "gray", linetype = "dotted"),
                 panel.grid.minor = element_line(color = "lightgray", linetype = "dotted")) +
           geom_abline(intercept = 0, slope = 1, linetype = "dotted")
  )
  
  
}



plot_grid(ncol = 3, r1, r2, r3, r4,
          r5, r6, r7, r8, r9)


#just for T1 prospermatogonia

number <- 5
current.comparison <- list.of.comparisons[[number]]
# current.comparison <- marker.genes.ssc
current.counts <- av.counts.hs.df[,c(number, (length(list.of.cell.types))+1)]

current.counts2 <- av.counts.ncg.df[,c(number, (length(list.of.cell.types))+1)]


current.comparison$gene <- rownames(current.comparison)
# current.comparison.filtered <- filter(current.comparison, p_val_adj < p.val.cutoff)
current.comparison.filtered.up <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC > log2fc.cutoff)
current.comparison.filtered.down <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC < -log2fc.cutoff)
current.comparison.dataframe <- data.frame(current.counts[,1], current.counts2[,1])
rownames(current.comparison.dataframe) <- rownames(current.counts)
current.comparison.subset <- current.comparison.dataframe[rowSums(current.comparison.dataframe[])>0,]
# current.comparison.sig.up <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.up$gene)
# current.comparison.sig.down <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.down$gene)
current.comparison.subset <- na.omit(current.comparison.subset)
current.comparison.subset$gene.name <- rownames(current.comparison.subset)
current.comparison.subset$DEGs <- "Not significant"
current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene] <- "Higher in vivo"
current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene] <- "Higher in xrTestis"
table(current.comparison.subset$DEGs)
print(paste0(number, " ", list.of.cell.types[number], " HS vs NCG", 
             " Up: ", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]), 
             " Down:", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])
))


colnames(current.comparison.subset) <- c("HS.counts", "NCG.counts", "gene.name", "DEGs")
ml = lm(HS.counts~NCG.counts, data = current.comparison.subset) 

# Extracting R-squared parameter from summary 
r.squared.value <- summary(ml)$r.squared 

current.comparison.subset$HS.counts <- (current.comparison.subset$HS.counts)*6
current.comparison.subset$NCG.counts <- (current.comparison.subset$NCG.counts)*6


assign(paste0("p", number), (ggplot() + geom_point(data = current.comparison.subset, aes(x = HS.counts, y = NCG.counts, color = DEGs)) + 
                               # scale_color_manual(values=c('#E69F00','#999999', '#56B4E9')) +
                               scale_color_manual(values=c(viridis.colours[3], viridis.colours[4],'#999999')) +
                               theme_classic() + labs(title=paste(list.of.cell.types[number],", q =", p.val.cutoff, ", log2fc =", log2fc.cutoff, "r^2 = ", round(r.squared.value, 2), 
                                                                  "n up = ", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]),
                                                                  "n down = ", " Down:", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])))
                             # x=list.of.cell.types[number], y = list.of.cell.types[number]) + 
                             #theme(plot.title = element_text(size=9))) 
))
p5




gene.of.interest <- c("MTRNR2L1", "MAGEA6", "HOXA5", "HOXB8", "HOXB5", "MTRNR2L8", "MAGEF1", "APOE", "LY6K", "XIST", "LY6E", "APOA1")

gene.of.interest <- current.comparison.subset$gene.name[current.comparison.subset$DEGs %in% c("Higher in vivo", "Higher in xrTestis")]

gene.of.interest <- c("MTRNR2L1", "MAGEB2", "ID2", "PLCG2", 
                      "ID3", "ID1", "HIST1H4C", 
                      "XIST", "APOE", "HSPA5", "KIT", "LY6E", "TMED9")

# Create the ggplot object
p <- ggplot(data = current.comparison.subset, aes(x = HS.counts, y = NCG.counts, label = gene.name, fill = DEGs)) +
  geom_point(aes(size = ifelse(DEGs == 'Not significant', 0.7, 1.7), color = ifelse(DEGs == 'Not significant', '#aaa9aa', 'black')), shape = 21) +
  geom_text(data = base::subset(current.comparison.subset, gene.name %in% gene.of.interest), hjust=-0.1, vjust=-0.1) +
  scale_fill_manual(values = c("#334796", "#eb2e30", '#00000000')) +
  scale_color_manual(values = c("#aaa9aa", "black", 'black')) +
  scale_size_identity() +
  theme_minimal() +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 15)) +
  labs(title = paste(list.of.cell.types[number], ", q =", p.val.cutoff, ", log2fc =", log2fc.cutoff, "r^2 =", round(r.squared.value, 2),
                     "n up =", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]),
                     "n down =", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])))

# Set custom breaks for x and y axes, and remove minor breaks
p <- p + scale_x_continuous(breaks = seq(0, 15, by = 1), minor_breaks = NULL)
p <- p + scale_y_continuous(breaks = seq(0, 15, by = 1), minor_breaks = NULL)

# Add a black border around the plot
p <- p + theme(
  panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
  panel.grid.major = element_line(color = "gray", linetype = "dotted", linewidth = 0.5),
  axis.ticks = element_line(color = "black"),
  axis.text = element_text(margin = margin(0.1, 0.1, 0.1, 0.1))
)

# Add diagonal line for y = x (average value)
p <- p + geom_abline(intercept = 0, slope = 1, linetype = "solid", color = "black")

# Add diagonal line for y = log2fc.cutoff * x
p <- p + geom_abline(intercept = -2, slope = 1, linetype = "dotted", color = "black")
p <- p + geom_abline(intercept = 2, slope = 1, linetype = "dotted", color = "black")
p <- p + geom_abline(intercept = -1, slope = 1, linetype = "dotted", color = "black")
p <- p + geom_abline(intercept = 1, slope = 1, linetype = "dotted", color = "black")

# Print the plot
print(p)



#repeat, this time for a specific gene

gene.of.interest <- c("IGFBP2",	"PLCG2",	"PLCG2",	"PLCG2","MT-ND2",	"FTL",	"MDK",	"RPS4X",	"MDK")
gene.of.interest <- c("CHCHD2",	"RPS26",	"RPS26",	"HBG2",	"RPS26",	"MORC1",	"HIST1H2AA",	"PRM1",	"PRM2")

gene.of.interest <- "STRA8"
gene.of.interest <- c("MTRNR2L1", "MAGEA6", "HOXA5", "HOXB8", "HOXB5", "MTRNR2L8", "MAGEF1", "APOE", "LY6K", "XIST", "LY6E", "APOA1")


for(number in 1:(length(list.of.cell.types))){
  # number <- 5
  current.comparison <- list.of.comparisons[[number]]
  # current.comparison <- marker.genes.ssc
  current.counts <- av.counts.hs.df[,c(number, (length(list.of.cell.types))+1)]
  
  current.counts2 <- av.counts.ncg.df[,c(number, (length(list.of.cell.types))+1)]
  
  
  current.comparison$gene <- rownames(current.comparison)
  # current.comparison.filtered <- filter(current.comparison, p_val_adj < p.val.cutoff)
  current.comparison.filtered.up <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC > log2fc.cutoff)
  current.comparison.filtered.down <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC < -log2fc.cutoff)
  current.comparison.dataframe <- data.frame(current.counts[,1], current.counts2[,1])
  rownames(current.comparison.dataframe) <- rownames(current.counts)
  current.comparison.subset <- current.comparison.dataframe[rowSums(current.comparison.dataframe[])>0,]
  # current.comparison.sig.up <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.up$gene)
  # current.comparison.sig.down <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.down$gene)
  current.comparison.subset <- na.omit(current.comparison.subset)
  current.comparison.subset$gene.name <- rownames(current.comparison.subset)
  current.comparison.subset$DEGs <- "Not significant"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene] <- "Up"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene] <- "Down"
  current.comparison.subset2 <- current.comparison.subset
  current.comparison.subset2$DEGs[current.comparison.subset2$gene.name %in% gene.of.interest] <- "GeneOfInterest"
  table(current.comparison.subset2$DEGs)
  print(paste0(number, " ", list.of.cell.types[number], " HS vs NCG", 
               " Up: ", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]), 
               " Down:", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])
  ))
  for(gene.number in 1:length(gene.of.interest)){
    # gene.number <- 1
    current.gene <-(gene.of.interest[gene.number])
    print(current.gene)
    print(current.comparison.subset[current.gene, "DEGs"])
  }
  
  colnames(current.comparison.subset) <- c("HS.counts", "NCG.counts", "gene.name", "DEGs")
  assign(paste0("r", number), ggplot(data = current.comparison.subset, aes(x = HS.counts, y = NCG.counts, 
                                                                           label = gene.name,
                                                                           color = DEGs)) + geom_point() + geom_text(data = base::subset(current.comparison.subset, gene.name %in% gene.of.interest), hjust=-0.1, vjust=-0.1) +
           scale_color_manual(values=c('#E69F00', "red4",'#999999', '#56B4E9')) + 
           theme_classic() + labs(title=paste(list.of.cell.types[number],", q", p.val.cutoff, ", log2fc", log2fc.cutoff, 
                                              " Up: ", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]), 
                                              " Down:", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])) )
         # geom_label(data = base::subset(current.comparison.subset, DEGs == "Up")))
  )
  
}

#T1 = r5

r5

plot_grid(ncol = 3, r1, r2, r3, r4, 
          r5, r6, r7, r8, 
          r9)



#Methylation related enzymes

# enzymes that deposit H3K9me2
H3K9methylases <- c("EHMT2", "SUV39H1", "SUV39H2")
H3K9demethylases <- c("KDM1A", "KDM3A", "KDM3B", "KDM4A")
H3K27methylases <- c("ED", "EZH2", "SUZ12","EZH1")
DNAdemethylases <- c("TET1", "TET2", "TET3", "DNMT1", "DNMT3A", "UHRF1", "DNMT3L", "DNMT3B")

p1 <- FeaturePlot(HS_NCG_only_HS, features = DNAdemethylases, order = T, min.cutoff = 0, ncol = length(DNAdemethylases)) & scale_color_viridis(option="plasma") & NoLegend() & DarkTheme()
p2 <- FeaturePlot(HS_NCG_only_NCG, features = DNAdemethylases, order = T, min.cutoff = 0, ncol = length(DNAdemethylases)) & scale_color_viridis(option="plasma") & NoLegend() & DarkTheme()
plot_grid(p1, p2, nrow = 2)

p1 <- VlnPlot(HS_NCG, features = "DNMT1", group.by = "cell.type", split.by = "experiment")
p1

Idents(HS_NCG) <- "experiment"
HS_NCG_noPGCLCs <- subset(HS_NCG, idents = "hPGCLC", invert = T)
VlnPlot(HS_NCG_noPGCLCs, features = DNAdemethylases, ncol = 3, split.plot = T, group.by = "cell.type", split.by = "experiment")

VlnPlot(HS_NCG_noPGCLCs, features = DNAdemethylases, ncol = 3, split.plot = TRUE, group.by = "cell.type", split.by = "experiment") & theme(axis.text.x = element_blank())

# Modify the plot to remove x-axis labels for all but the bottom row
vln_plot + theme(axis.text.x = element_blank()) + facet_wrap(~ experiment, ncol = 2, scales = "free_x")



FeaturePlot(HS_NCG_noPGCLCs, features = c("REC8", "TEX101"), split.by = "experiment") & scale_colour_viridis() & DarkTheme()

genes.1 <- c("DDX4", "DAZL", "SYCP3", "SYCP1", "TEX101", "REC8") 
genes.2 <- c("MAGEA3", "MAGEA4", "MEIOB", "MEIOC","SCML1", "SOHLH1")
genes.3 <- c("MAGEA3", "MAGEA4", "MEIOB", "MEIOC","SCML1", "SOHLH1")
genes.4 <- c("POU5F1", "UTF1", "MORC1", "PIWIL4","MEIOC", "SYCP1")

genes.5 <- c("DDX4", "MEIOC", "SYCP1", "MEIOB")

Ryo.genes <- c("TFAP2C", "NANOG", "PIWIL4","POU5F1", "DDX4", 
               "DAZL", "SOX17", "MAGEC2", "MAGEA3", "UTF1", "GFRA1", "MKI67", "KIT", "H2AFX", 
               "REC8", "SYCP1", "SYCP3")

FeaturePlot(HS_NCG, features = Ryo.genes, order = T, ncol = 6) & scale_colour_viridis(option = "plasma") & DarkTheme()
DimPlot(HS_NCG, label = T)

VlnPlot(HS_NCG_noPGCLCs, features = genes.5, 
        # ncol = 3, 
        split.plot = T, pt.size = 0,
        group.by = "cell.type", split.by = "experiment")

VlnPlot(HS_NCG_noPGCLCs, features = genes.4, ncol = 3, split.plot = T, pt.size = 0,
        group.by = "cell.type", split.by = "experiment")



VlnPlot(HS_NCG_only_HS, features = genes.1)
VlnPlot(HS_NCG_only_HS, features = genes.2)
VlnPlot(HS_NCG_only_NCG, features = genes.1)
VlnPlot(HS_NCG_only_NCG, features = genes.2)

Idents(HS_NCG_noPGCLCs) <- "cell.type"
VlnPlot(HS_NCG_noPGCLCs, features = c( "MEIOC"), split.by = "experiment", pt.size = 0)

Stacked_VlnPlot(seurat_object = HS_NCG_noPGCLCs, split.by = "experiment", features = "MEIOC", x_lab_rotate = TRUE)

FeaturePlot(HS_NCG, features = c("TFAP2C", "PIWIL4", "DDX4"), ncol = 3, order = T) & scale_colour_viridis(option = "magma") & DarkTheme()

DimPlot(HS_NCG_only_NCG)
table(Idents(HS_NCG_only_NCG))


FeaturePlot(HS_NCG, features = "XIST", split.by = "experiment", min.cutoff = 0, order = T) & DarkTheme() & scale_color_viridis()

DimPlot(HS_NCG_only_HS, split.by = "stage")
Idents(HS_NCG_only_HS) <- "stage"
HS_NCG_only_HS_12yo <- subset(HS_NCG_only_HS, idents = "12 yo")
HS_NCG_only_HS_5yo <- subset(HS_NCG_only_HS, idents = "5 yo")
HS_NCG_only_HS_fixed <- HS_NCG_only_HS[, !is.na(HS_NCG_only_HS$stage)]
HS_NCG_only_HS_prenatal <- subset(HS_NCG_only_HS_fixed, idents = c("12 yo", "5 yo"), invert = T)
Idents(HS_NCG_only_HS_12yo) <- HS_NCG_only_HS_12yo[["cell.type"]]
Idents(HS_NCG_only_HS_5yo) <- HS_NCG_only_HS_5yo[["cell.type"]]
Idents(HS_NCG_only_HS_prenatal) <- HS_NCG_only_HS_prenatal[["cell.type"]]

table(Idents(HS_NCG_only_HS_12yo))
table(Idents(HS_NCG_only_HS_5yo))
table(Idents(HS_NCG_only_HS_prenatal))


DimPlot(HS_NCG_noPGCLCs)
HS_NCG_noPGCLCs_nopostmeiotic <- subset(HS_NCG_noPGCLCs, idents = c("Leptotene/Zygotene", 
                                                                    "Pachytene spermatocytes",
                                                                    "Diplotene/2ndary spermatocytes", 
                                                                    "Spermatids"), invert = T)
VlnPlot(HS_NCG_noPGCLCs_nopostmeiotic, features = "MEIOC", split.by = "experiment", split.plot = T)


### Extended Data Fig 6A ###


p2 <- FeaturePlot(HS_NCG, features = c("DDX4"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "red", "orangered", "yellow", "white"))
p1 <- FeaturePlot(HS_NCG, features = c("TFAP2C"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "darkgreen", "green", "yellow", "white"))
p3 <- FeaturePlot(HS_NCG, features = c("PIWIL4"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "cyan4", "cyan", "white"))

plot_grid(p1, p2, p3, ncol = 3)

p4 <- FeaturePlot(HS_NCG, features = c("MAGEB2"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "blue2", "deepskyblue", "skyblue"))
p5 <- FeaturePlot(HS_NCG, features = c("MAGEA4"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "magenta4", "hotpink", "white"))
p6 <- FeaturePlot(HS_NCG, features = c("SYCP1"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "brown", "burlywood"))


p7 <- FeaturePlot(HS_NCG, features = c("MEIOC"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "orange4", "orange", "yellow"))
p8 <- FeaturePlot(HS_NCG, features = c("ACRV1"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "plum4", "plum", "white"))
p9 <- FeaturePlot(HS_NCG, features = c("PRM1"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "olivedrab", "olivedrab1", "yellow"))




plot_grid(p1, p2, p3,p4,p5,p6,p7,p8,p9,
          ncol = 3)

### EXTENDED DATA FIG 7 ###


Idents(HS_NCG) <- "cell.type"
DimPlot(HS_NCG)
spermatogonia_comparison <- subset(HS_NCG, idents = c("T1 prospermatogonia", "Spermatogonia", "Diff. spermatogonia (early)"))
spermatogonia_comparison <- FindVariableFeatures(object = spermatogonia_comparison, nfeatures = 2000, verbose = FALSE)
# DefaultAssay(object = current.integrated) <- "integrated"
spermatogonia_comparison <- ScaleData(object = spermatogonia_comparison, vars.to.regress = c("nCount_RNA", "percent.mito"), verbose = FALSE)
spermatogonia_comparison <- RunPCA(object = spermatogonia_comparison,
                                   npcs = 30,
                                   verbose = TRUE)
spermatogonia_comparison <- FindNeighbors(object = spermatogonia_comparison)
spermatogonia_comparison <- FindClusters(object = spermatogonia_comparison,
                                         resolution = 1.1,
                                         #algorithm = 1,
                                         verbose = TRUE
)
# current.integrated <- RunTSNE(object = current.integrated, dims = 1:20)
spermatogonia_comparison <- RunUMAP(object = spermatogonia_comparison, reduction = "pca",
                                    dims = 1:30)

DimPlot(spermatogonia_comparison, group.by = "experiment")

FeaturePlot(spermatogonia_comparison, features = c("GAGE12H", "KIT", "SOHLH2","PIWIL4", "ID4", "ETV5", "SOHLH1", "EOMES", "TSPAN8"), order = T) & DarkTheme() & scale_color_viridis(option = "C")



# spermatogonia_comparison$PhaseExperiment <- paste0(spermatogonia_comparison$Phase, spermatogonia_comparison$experiment)
DimPlot(spermatogonia_comparison, group.by = "experiment")
current.list <- Seurat::SplitObject(spermatogonia_comparison, split.by = c("experiment"))
table(spermatogonia_comparison)

# merged.sample <- merge(x = current.list[[2]], y = current.list[[3]])
# current.list2 <- c(merged.sample, current.list[[1]], current.list[[4]], current.list[[5]], current.list[[6]])

for (i in 1:length(x = current.list)) {
  current.list[[i]] <- NormalizeData(object = current.list[[i]], verbose = FALSE)
  current.list[[i]] <- FindVariableFeatures(object = current.list[[i]], nfeatures = 2000, verbose = FALSE)
}

current.anchors <- FindIntegrationAnchors(object.list = current.list, dims = 1:30)
# save(current.anchors, file="current.anchors.Robj")
current.integrated <- IntegrateData(anchorset = current.anchors, dims = 1:30)
# save(current.integrated, file="current.integrated.Robj")

DefaultAssay(object = current.integrated) <- "integrated"
current.integrated <- ScaleData(object = current.integrated, vars.to.regress = c("nCount_RNA", "percent.mito"), 
                                verbose = FALSE)
current.integrated <- RunPCA(object = current.integrated,
                             npcs = 30,
                             verbose = TRUE)
current.integrated <- FindNeighbors(object = current.integrated)
current.integrated <- FindClusters(object = current.integrated,
                                   resolution = 1.1,
                                   #algorithm = 1,
                                   verbose = TRUE
)
# current.integrated <- RunTSNE(object = current.integrated, dims = 1:20)
current.integrated <- RunUMAP(object = current.integrated, reduction = "pca",
                              dims = 1:30)



DefaultAssay(current.integrated) <- "RNA"
DimPlot(current.integrated, group.by = "Phase")




spermatogonia_comparison2 <- current.integrated
DefaultAssay(spermatogonia_comparison2) <- "RNA"

# c("GAGE12H", "KIT", "SOHLH2","PIWIL4", "ID4", "ETV5", "SOHLH1", "EOMES", "TSPAN8")

FeaturePlot(spermatogonia_comparison2, features = c("GAS5", "TOMM7", "APOE", "SERPING1", "TPM2", "MEG3", "NANOS2", 
                                                    "GAGE12H", "RET", "PIWIL4", "EGR4", "ID4","UTF1",  "SOHLH1", "EOMES", "TSPAN33","KIT", "SOHLH2",
                                                    "GFRA1", "L1TD1", "ETV5", "MSL3"), order = T) & DarkTheme() & scale_color_viridis(option = "C")

T1markers <- c("GAS5", "TOMM7", "SERPING1", "TPM2", "MEG3", "NANOS2")
S0markers <- c("EGR4", "PIWIL4", "SLC25A22", "ICA1L", "PPP1R36", "EGR4", "MSL3", "TSPAN33")
S1markers <- c("L1TD1", "ETV5", "GFRA1")
Diffmarkers <- c("KIT", "SOHLH2") 

spermatogonia_comparison2 <- AddModuleScore(spermatogonia_comparison2,
                                            features = list(T1markers),
                                            name="T1markers")


spermatogonia_comparison2 <- AddModuleScore(spermatogonia_comparison2,
                                            features = list(S0markers),
                                            name="S0markers")


spermatogonia_comparison2 <- AddModuleScore(spermatogonia_comparison2,
                                            features = list(S1markers),
                                            name="S1markers")


spermatogonia_comparison2 <- AddModuleScore(spermatogonia_comparison2,
                                            features = list(Diffmarkers),
                                            name="Diffmarkers")


FeaturePlot(spermatogonia_comparison2, split.by = "experiment",  features = c("T1markers1", "S0markers1", "S1markers1", "Diffmarkers1"), order = T) & DarkTheme() & scale_color_viridis(option = "C")

VlnPlot(spermatogonia_comparison2, ncol = 4, features = c("NANOS2", "TSPAN33", "L1TD1", "KIT"), pt.size = 0, split.plot = T, split.by = "experiment", group.by = "cell.type.spermatogonia")

# setwd("/Users/ewhelan/Documents/scRNAseq/Rscripts/savefiles")
# save(spermatogonia_comparison2, file = "spermatogonia_comparison2.Robj")

# # # 



Idents(spermatogonia_comparison2) <- "seurat_clusters"
DimPlot(spermatogonia_comparison2, label = T)

new.cluster.ids <- character(length(levels(Idents(spermatogonia_comparison2))))

T1.ids <- c(5,9) +1
SSC.ids <- c(8, 4, 3) +1
progenitor.ids <- c(12,0,7,10,2,11) +1
differentiating.ids <- c(6,1)+ 1

new.cluster.ids[T1.ids] <- "T1.prospg"
new.cluster.ids[SSC.ids] <- "S0.spg"
new.cluster.ids[progenitor.ids] <- "S1.spg"
new.cluster.ids[differentiating.ids] <- "Differentiating"



new.cluster.ids

names(x = new.cluster.ids) <- levels(x = spermatogonia_comparison2)
spermatogonia_comparison2 <- RenameIdents(object = spermatogonia_comparison2, new.cluster.ids)
table(Idents(spermatogonia_comparison2))


DimPlot(spermatogonia_comparison2)
spermatogonia_comparison2$cell.type.spermatogonia <- Idents(spermatogonia_comparison2)

spermatogonia_comparison2$cell.type.spermatogonia <- factor(x = spermatogonia_comparison2$cell.type.spermatogonia, levels = c("T1.prospg", 
                                                                                                                              "S0.spg",
                                                                                                                              "S1.spg", "Differentiating"))



Idents(spermatogonia_comparison2) <- spermatogonia_comparison2$cell.type.spermatogonia
DimPlot(spermatogonia_comparison2, group.by = "cell.type.spermatogonia", label = T)

setwd("~/Documents/scRNAseq/Robjects")
save(spermatogonia_comparison2, file = "spermatogonia_comparison2.Robj")
load("spermatogonia_comparison2.Robj")

table(spermatogonia_comparison2$cell.type.spermatogonia, spermatogonia_comparison2$experiment)


cds <- as.cell_data_set(spermatogonia_comparison2)


cds <- cluster_cells(cds, reduction_method = c("UMAP"), #cluster_method = "louvain",
                     resolution = 9e-5, 
                     k = 8,
                     partition_qval = 0.05
)

plot_cells(cds, show_trajectory_graph = FALSE)
cds <- learn_graph(cds, 
                   close_loop = F,
                   # learn_graph_control=list(ncenter=300, minimal_branch_len=7), #250, 18 is pretty good
                   use_partition = F)
print(plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)) #+ ggtitle(a))

# Figure 2C #######################################
cds <- order_cells(cds)#, root_cells = max.avp)
plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE, 
           label_branch_points = FALSE)+ NoLegend()



pseudotime.values <- pseudotime(cds, reduction_method = "UMAP")


Idents(spermatogonia_comparison2) <- "experiment"
spermatogonia_comparison_xrTestis <- subset(spermatogonia_comparison2, idents = "transplant")
cell.names <- Cells(spermatogonia_comparison_xrTestis)


NCG_germ4_newcells <- NCG_germ4
new.cell.names <- (Cells(NCG_germ4_newcells))
old.cell.names <- (Cells(NCG_germ4))

target.cell.names <- gsub("_\\d+$", "", cell.names)
table.of.cell.names <- data.frame(old.cell.names, new.cell.names)



# Create a vector of replacement values
replacement_vector <- table.of.cell.names$new.cell.names[match(table.of.cell.names$old.cell.names, target.cell.names)]
replacement.cell.names <- replacement_vector[!is.na(replacement_vector)]
# Print the replacement vector
print(replacement.cell.names)
length(replacement.cell.names)


# Find matching indices
matching_indices <- match(target.cell.names, table.of.cell.names$old.cell.names)

# Replace values in target.cell.names with values from new.cell.names
replacement_vector <- table.of.cell.names$new.cell.names[matching_indices]

head(target.cell.names)
head(replacement_vector)


spermatogonia_comparison_xrTestis2 <- RenameCells(spermatogonia_comparison_xrTestis, new.names = replacement_vector)
DimPlot(spermatogonia_comparison_xrTestis2)
head(Cells(spermatogonia_comparison_xrTestis2))

cells.spg <- Cells(spermatogonia_comparison_xrTestis2)

cells.spg.fixed <- gsub("-", "_", cells.spg)

spermatogonia_comparison_xrTestis3 <- RenameCells(object = spermatogonia_comparison_xrTestis2, new.names = cells.spg)



setwd("~/Documents/scRNAseq/Rscripts/savefiles")
save(spermatogonia_comparison_xrTestis2, file = "spermatogonia_comparison_xrTestis2.Robj")
save(spermatogonia_comparison_xrTestis3, file = "spermatogonia_comparison_xrTestis3.Robj")

cols = c(#"#a6cee2", 
  "#2179b4",  "#b4d88b", "#36a047", "#f6999a", "#f15a29", "#fdbf6f", "#f7941d", "#cab3d6", "#6b3f98")

cols2 = c(#"#a6cee2", 
  "#f15a29", "#fdbf6f", "darkgoldenrod3", "#f7941d")



Idents(spermatogonia_comparison2) <- "experiment"
DimPlot(spermatogonia_comparison2)

### Extended Data Fig 7 A - C ###

DimPlot(spermatogonia_comparison2, #split.by = "experiment", 
        group.by = "cell.type.spermatogonia", cols = cols2)

spermatogonia_comparison_inVivo <- subset(spermatogonia_comparison2, idents = "human in vivo")
spermatogonia_comparison_xrTestis <- subset(spermatogonia_comparison2, idents = "transplant")

Idents(spermatogonia_comparison_inVivo) <- "cell.type.spermatogonia"
Idents(spermatogonia_comparison_xrTestis) <- "cell.type.spermatogonia"

table(Idents(spermatogonia_comparison_inVivo))
table(Idents(spermatogonia_comparison_xrTestis))

DimPlot(spermatogonia_comparison_inVivo)

#I am here

all.genes <- c("MAGEA4", "ID4", "FGFR3", 
               "TCF3",  "UTF1", "GFRA1", 
               "ETV5", "L1TD1", "KIT", 
               "MKI67", "DMRT1", "STRA8",
               "PIWIL4", "PHGDH", "SLC25A22", 
               "ICA1L", "PPP1R36", "MAGEB1", 
               "EGR4", "MSL3", "TSPAN33")

T1markers <- c("GAS5", "TOMM7", "SERPING1", "TPM2", "MEG3", "NANOS2")
S0markers <- c("EGR4", "PIWIL4", "SLC25A22", "ICA1L", "PPP1R36", "EGR4", "MSL3", "TSPAN33")
S1markers <- c("L1TD1", "ETV5", "FAM25G")
Diffmarkers <- c("KIT", "SOHLH2") 
FeaturePlot(spermatogonia_comparison2, order = T, features = S0markers, min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(spermatogonia_comparison2, order = T, features = c("ASB9", "CYP26A1", "SAT1", "CITED2", "FAM25G", "LINC01419", "HOXA1", "PLBD1"), min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")

FeaturePlot(spermatogonia_comparison2, features = c("NANOS2", "EGR4", "PIWIL4", "MSL3", "L1TD1", "ETV5", "SOHLH2", "KIT"), min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")

Idents(spermatogonia_comparison2) <- "cell.type.spermatogonia"
DimPlot(spermatogonia_comparison2)
S0.vs.S1 <- FindMarkers(spermatogonia_comparison2, ident.1 = "S0.spg", ident.2 = "S1.spg", logfc.threshold = log2fc.cutoff, min.pct = 0)
S0.vs.S1$difference <- S0.vs.S1$pct.1 - S0.vs.S1$pct.2

# Order the dataframe by the difference column
S0.vs.S1.ordered <- S0.vs.S1[order(S0.vs.S1$difference, decreasing = F), ]
list.of.genes <- rownames(S0.vs.S1.ordered[49:72,])

list.of.genes <- c("NANOS2", "EGR4", "PIWIL4", "MSL3", "L1TD1", "ETV5", "DUSP6", "KIT") # primary figure
list.of.genes2 <- c("TSPAN33", "KLC3", "RGS14", "C19orf84", "MAPK8IP1",
                    "SPRY1", "FAM25G", "RARB", "ASB9", "LINC01091"
)

list.of.genes.all <- c(list.of.genes, list.of.genes2)

### Extended Data Fig 7D ###


FeaturePlot(ncol = 8, spermatogonia_comparison, order = T, features = list.of.genes, min.cutoff = 0) &
  DarkTheme() &
  scale_color_viridis(option = "A")

Idents(spermatogonia_comparison2)
DimPlot(spermatogonia_comparison2)


p1 <- FeaturePlot(ncol = 8, pt.size =2, spermatogonia_comparison_inVivo, order = T, features = list.of.genes, min.cutoff = 0) &
  DarkTheme() &
  scale_color_viridis(option = "A")

p2 <- FeaturePlot(ncol = 8, pt.size =2, spermatogonia_comparison_xrTestis, order = T, features = list.of.genes, min.cutoff = 0) &
  DarkTheme() &
  scale_color_viridis(option = "A")

p3 <- FeaturePlot(ncol = 10, spermatogonia_comparison_inVivo, order = T, features = list.of.genes2, min.cutoff = 0) &
  DarkTheme() &
  scale_color_viridis(option = "A")

p4 <- FeaturePlot(ncol = 10, spermatogonia_comparison_xrTestis, order = T, features = list.of.genes2, min.cutoff = 0) &
  DarkTheme() &
  scale_color_viridis(option = "A")


plot_grid(p1, p2, p3, p4, ncol = 1)
plot_grid(p1, p2, ncol = 1)
plot_grid(p3, p4, ncol = 1)


T1.markers <- FindMarkers(spermatogonia_comparison2, ident.1 = "T1.prospg", min.pct = 0, logfc.threshold = 0)
write.csv(T1.markers, file = "T1.markers.csv")

genes1 <- c("POU5F1", "NANOG", "TCL1A", "PDPN", "PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1",
            "UTF1", "KIT", "ASB9", "NANOS2", "PIWIL4", "MORC1", "MAGEB2", "TDRD1", "RHOXF1", "MAGEC2", "ZBTB16", "EGR4", "RET", "MAGEA4",
            "SOHLH1", "DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8", "DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9", "ZCWPW1", "RAD51AP2", "SYCP3",
            "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1")

genes2 <- c("GAS5",     "TOMM7",    "APOE",     "SERPING1", "TPM2",     "MEG3",     "NANOS2",   "GAGE12H",  "RET" ,     "PIWIL4",   "EGR4",    "ID4",     "UTF1",     "SOHLH1" ,  "EOMES" ,  
            "TSPAN33" , "KIT" ,    "SOHLH2"  , "GFRA1"  ,  "L1TD1"  ,  "ETV5" ,    "MSL3" ,    "SLC25A22", "ICA1L"  ,  "PPP1R36" , "MAGEA4" ,  "FGFR3"  ,  "TCF3" ,    "MKI67" ,   "DMRT1"   ,
            "STRA8",    "PHGDH" ,   "MAGEB1")

genes3 <- c("NANOS2", "TCF3", "UTF1", "MAGEC2", "SIX1", "DDX4",
            
            "SOHLH1", "PIWIL4", "EGR4", "TSPAN33", 
            
            "L1TD1", "ETV5", "MAGEA4", "FGFR3", "DMRT1",
            
            "KIT", "STRA8", "REC8", "SOHLH2")

genes4 <- c("DUSP6", "ETV5", "L1TD1", 
            "SPRY1", "STOX1", "SAT1", 
            "FAM25G", "RARB", "CLEC4C", 
            "ASB9", "LINC01091", 
            "TSPAN33", "KLC3", "RGS14", 
            "C19orf84", "MAPK8IP1",
            "HOXA1", "PLBD1")

genes <- unique(c(genes1, genes2, genes3, genes4))
# genes <- c("NANOG", "SOHLH1")

genes <- unique(c("NANOS2", "BCAM", "TCF3", "UTF1", "MAGEC2", "SIX1", "DDX4",
                  
                  "SOHLH1", "PIWIL4", "EGR4", "TSPAN33", "MORC1", "MAGEB2", "PHGDH", "KLC3", "RGS14", "C19orf84", "MAPK8IP1",
                  
                  "L1TD1", "ETV5", "MAGEA4", "FGFR3", "DMRT1", "TDRD1", "MAGEA4", "SPRY1", "FAM25G",
                  
                  "KIT", "STRA8", "REC8", "SOHLH2", "NANOS3", "PDPN", "DMRT1", "DUSP6", "SAT1", "RARB", "ASB9", "LINC01091", "PLBD1"))

genes <- unique(c("NANOS2", 
                  # "TCF3", 
                  # "UTF1", 
                  "BCAM", 
                  "MAGEC2", 
                  
                  # "DDX4",
                  "PIWIL4", 
                  "EGR4", 
                  "TSPAN33", 
                  # "MORC1", 
                  # "MAGEB2", 
                  "PHGDH", 
                  "KLC3", 
                  "RGS14", 
                  "C19orf84", 
                  "MAPK8IP1",
                  "SOHLH1",
                  "SIX1", 
                  
                  # "L1TD1", 
                  
                  "MAGEA4", 
                  # "FGFR3", 
                  "TDRD1",
                  # "MAGEA4", 
                  "SPRY1", 
                  "FAM25G",
                  "ETV5", 
                  
                  # "DMRT1", 
                  "KIT", 
                  # "STRA8",
                  # "REC8",
                  # "SOHLH2", 
                  # "NANOS3", 
                  "PDPN", 
                  # "DMRT1", 
                  # "DUSP6", 
                  "SAT1",
                  # "RARB",
                  "ASB9", 
                  # "LINC01091", 
                  "PLBD1"))



# genes <- genes4

DefaultAssay(spermatogonia_comparison_inVivo) <- "RNA"
DefaultAssay(spermatogonia_comparison_xrTestis) <- "RNA"
Idents(spermatogonia_comparison_inVivo) <- "cell.type.spermatogonia"
Idents(spermatogonia_comparison_xrTestis) <- "cell.type.spermatogonia"


# DimPlot(spermatogonia_comparison_inVivo)
# all.markers <- FindAllMarkers(spermatogonia_comparison_inVivo, min.pct = 0, logfc.threshold = 0, features = genes)

T1.markers <- FindMarkers(spermatogonia_comparison_inVivo, ident.1 = "T1.prospg", min.pct = 0, logfc.threshold = 0, features = genes)
S0.markers <- FindMarkers(spermatogonia_comparison_inVivo, ident.1 = "S0.spg", min.pct = 0, logfc.threshold = 0, features = genes)
S1.markers <- FindMarkers(spermatogonia_comparison_inVivo, ident.1 = "S1.spg", min.pct = 0, logfc.threshold = 0, features = genes)
Diff.markers <- FindMarkers(spermatogonia_comparison_inVivo, ident.1 = "Differentiating", min.pct = 0, logfc.threshold = 0, features = genes)

T1.markers$gene <- rownames(T1.markers)
S0.markers$gene <- rownames(S0.markers)
S1.markers$gene <- rownames(S1.markers)
Diff.markers$gene <- rownames(Diff.markers)

T1.markers$cluster <- "T1.prospg"
S0.markers$cluster <- "S0.spg"
S1.markers$cluster <- "S1.spg"
Diff.markers$cluster <- "Differentiating"

all.markers <- rbind(T1.markers, S0.markers, S1.markers, Diff.markers)


order.of.genes <- as.vector(genes)
# graph.data$Pathway2 <- factor(graph.data$Pathway, levels = rev(unique(graph.data$Pathway)))
all.markers$gene <- factor(all.markers$gene, levels = rev(genes))
all.markers$cluster <- factor(all.markers$cluster, levels = c("T1.prospg", 
                                                              "S0.spg", 
                                                              "S1.spg",
                                                              "Differentiating"))

### DOT PLOT ### 



dotplot1 <- ggplot(data = all.markers, na.rm = T,
                   aes(x = cluster, y = (gene),
                       color = avg_log2FC, size = pct.1))+
  geom_point() +
  # scale_color_gradient(high = "red" , low = "gray") +
  scale_color_viridis(option = "inferno")+
  scale_size(range = c(-7, 7)) +
  theme_bw() +
  ylab("Gene") +
  xlab("Cell Type") +
  ggtitle("In vivo") +
  scale_y_discrete(position = "right")+
  theme(axis.text.x = element_text( #face = "bold",
    color = "black", size = 14, angle = -60, hjust = 0
  ),
  axis.text.y = element_text( #face = "bold",
    color = "black",  hjust = 1)
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "left")




# DimPlot(spermatogonia_comparison_xrTestis)
# all.markers <- FindAllMarkers(spermatogonia_comparison_xrTestis, min.pct = 0, logfc.threshold = 0, features = genes)

T1.markers <- FindMarkers(spermatogonia_comparison_xrTestis, ident.1 = "T1.prospg", min.pct = 0, logfc.threshold = 0, features = genes)
S0.markers <- FindMarkers(spermatogonia_comparison_xrTestis, ident.1 = "S0.spg", min.pct = 0, logfc.threshold = 0, features = genes)
S1.markers <- FindMarkers(spermatogonia_comparison_xrTestis, ident.1 = "S1.spg", min.pct = 0, logfc.threshold = 0, features = genes)
Diff.markers <- FindMarkers(spermatogonia_comparison_xrTestis, ident.1 = "Differentiating", min.pct = 0, logfc.threshold = 0, features = genes)

T1.markers$gene <- rownames(T1.markers)
S0.markers$gene <- rownames(S0.markers)
S1.markers$gene <- rownames(S1.markers)
Diff.markers$gene <- rownames(Diff.markers)

T1.markers$cluster <- "T1.prospg"
S0.markers$cluster <- "S0.spg"
S1.markers$cluster <- "S1.spg"
Diff.markers$cluster <- "Differentiating"



#now repeat for all markers

T1.markers.xrTestis <- FindMarkers(spermatogonia_comparison_xrTestis, ident.1 = "T1.prospg", min.pct = 0, logfc.threshold = 0)
S0.markers.xrTestis <- FindMarkers(spermatogonia_comparison_xrTestis, ident.1 = "S0.spg", min.pct = 0, logfc.threshold = 0)
S1.markers.xrTestis <- FindMarkers(spermatogonia_comparison_xrTestis, ident.1 = "S1.spg", min.pct = 0, logfc.threshold = 0)
Diff.markers.xrTestis <- FindMarkers(spermatogonia_comparison_xrTestis, ident.1 = "Differentiating", min.pct = 0, logfc.threshold = 0)

T1.markers.xrTestis$gene <- rownames(T1.markers.xrTestis)
S0.markers.xrTestis$gene <- rownames(S0.markers.xrTestis)
S1.markers.xrTestis$gene <- rownames(S1.markers.xrTestis)
Diff.markers.xrTestis$gene <- rownames(Diff.markers.xrTestis)

T1.markers.xrTestis$cluster <- "T1.prospg"
S0.markers.xrTestis$cluster <- "S0.spg"
S1.markers.xrTestis$cluster <- "S1.spg"
Diff.markers.xrTestis$cluster <- "Differentiating"



setwd("~/Desktop/S0S1/")
write.csv(T1.markers.xrTestis, row.names = FALSE, col.names = FALSE, file="T1.markers_xrTestis.csv")
write.csv(S0.markers.xrTestis, row.names = FALSE, col.names = FALSE, file="S0.markers_xrTestis.csv")
write.csv(S1.markers.xrTestis, row.names = FALSE, col.names = FALSE, file="S1.markers_xrTestis.csv")
write.csv(Diff.markers.xrTestis, row.names = FALSE, col.names = FALSE, file="Diff.markers_xrTestis.csv")

#Repeat for just in vivo and also the merged object

T1.markers.invivo <- FindMarkers(spermatogonia_comparison_inVivo, ident.1 = "T1.prospg", min.pct = 0, logfc.threshold = 0)
S0.markers.invivo <- FindMarkers(spermatogonia_comparison_inVivo, ident.1 = "S0.spg", min.pct = 0, logfc.threshold = 0)
S1.markers.invivo <- FindMarkers(spermatogonia_comparison_inVivo, ident.1 = "S1.spg", min.pct = 0, logfc.threshold = 0)
Diff.markers.invivo <- FindMarkers(spermatogonia_comparison_inVivo, ident.1 = "Differentiating", min.pct = 0, logfc.threshold = 0)

T1.markers.invivo$gene <- rownames(T1.markers.invivo)
S0.markers.invivo$gene <- rownames(S0.markers.invivo)
S1.markers.invivo$gene <- rownames(S1.markers.invivo)
Diff.markers.invivo$gene <- rownames(Diff.markers.invivo)

T1.markers.invivo$cluster <- "T1.prospg"
S0.markers.invivo$cluster <- "S0.spg"
S1.markers.invivo$cluster <- "S1.spg"
Diff.markers.invivo$cluster <- "Differentiating"

setwd("~/Desktop/S0S1/")
write.csv(T1.markers.invivo, row.names = FALSE, col.names = FALSE, file="T1.markers_inVivo.csv")
write.csv(S0.markers.invivo, row.names = FALSE, col.names = FALSE, file="S0.markers_inVivo.csv")
write.csv(S1.markers.invivo, row.names = FALSE, col.names = FALSE, file="S1.markers_inVivo.csv")
write.csv(Diff.markers.invivo, row.names = FALSE, col.names = FALSE, file="Diff.markers_inVivo.csv")


T1.markers.merged <- FindMarkers(spermatogonia_comparison2, ident.1 = "T1.prospg", min.pct = 0, logfc.threshold = 0)
S0.markers.merged <- FindMarkers(spermatogonia_comparison2, ident.1 = "S0.spg", min.pct = 0, logfc.threshold = 0)
S1.markers.merged <- FindMarkers(spermatogonia_comparison2, ident.1 = "S1.spg", min.pct = 0, logfc.threshold = 0)
Diff.markers.merged <- FindMarkers(spermatogonia_comparison2, ident.1 = "Differentiating", min.pct = 0, logfc.threshold = 0)

T1.markers.merged$gene <- rownames(T1.markers.merged)
S0.markers.merged$gene <- rownames(S0.markers.merged)
S1.markers.merged$gene <- rownames(S1.markers.merged)
Diff.markers.merged$gene <- rownames(Diff.markers.merged)

T1.markers.merged$cluster <- "T1.prospg"
S0.markers.merged$cluster <- "S0.spg"
S1.markers.merged$cluster <- "S1.spg"
Diff.markers.merged$cluster <- "Differentiating"

setwd("~/Desktop/S0S1/")
write.csv(T1.markers.merged, row.names = FALSE, col.names = FALSE, file="T1.markers_merged.csv")
write.csv(S0.markers.merged, row.names = FALSE, col.names = FALSE, file="S0.markers_merged.csv")
write.csv(S1.markers.merged, row.names = FALSE, col.names = FALSE, file="S1.markers_merged.csv")
write.csv(Diff.markers.merged, row.names = FALSE, col.names = FALSE, file="Diff.markers_merged.csv")










all.markers <- rbind(T1.markers, S0.markers, S1.markers, Diff.markers)


order.of.genes <- as.vector(genes)
# graph.data$Pathway2 <- factor(graph.data$Pathway, levels = rev(unique(graph.data$Pathway)))
all.markers$gene <- factor(all.markers$gene, levels = rev(genes))
all.markers$cluster <- factor(all.markers$cluster, levels = c("T1.prospg", 
                                                              "S0.spg", 
                                                              "S1.spg",
                                                              "Differentiating"))

### DOT PLOT ### 



dotplot2 <- ggplot(data = all.markers, na.rm = T,
                   aes(x = cluster, y = (gene), 
                       color = avg_log2FC, size = pct.1))+ 
  geom_point() +
  # scale_color_gradient(high = "red" , low = "gray") +
  scale_color_viridis(option = "inferno")+
  scale_size(range = c(-7, 7)) +
  theme_bw() + 
  ylab("Gene") + 
  xlab("Cell Type") + 
  ggtitle("xrTestis") +
  scale_y_discrete(position = "right")+
  theme(axis.text.x = element_text( #face = "bold",
    color = "black", size = 14, angle = -60, hjust = 0
  ),
  axis.text.y = element_text( #face = "bold",
    color = "black",  hjust = 1)
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "left")


plot_grid(dotplot1, dotplot2)



genes.of.interest <- c("ALDOA", "TPI1", "PFKL", "PFKFB3", "ENO3", "PFKM", "GAPDH" ,
                       "QRSL1", "MRPS17", "GFM2", "MRPS18C", "MRPS15", "MRPS18A", "HARS2", "HARS", "MRPS21", "MRPL47", "MTERF4", "LARS2", "NOA1" )



#read in GO data

setwd("~/Desktop/S0S1/")

S0upInVivo <- read.csv("S0_up_inVivo.csv")
S1upInVivo<- read.csv("S1_up_inVivo.csv")
S0upXrTestis <- read.csv("S0_up_xrTestis.csv")
S1upXrTestis <- read.csv("S1_up_xrTestis.csv")

genes.of.interest1 <- S1upXrTestis[7,6]
genes.of.interest2 <- S1upInVivo[36,6]


genes.of.interest1 <- unlist(strsplit(genes.of.interest1, ", "))
genes.of.interest2 <- unlist(strsplit(genes.of.interest2, ", "))

genes.of.interest <- unique(c(genes.of.interest1, genes.of.interest2))


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 

### SCATTER PLOT IN VIVO ### Extended Data Fig 7F

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### 


p.val.cutoff <- 0.05
log2fc.cutoff <- 0.20

Idents(spermatogonia_comparison_inVivo) <- "cell.type.spermatogonia"
# DimPlot(spermatogonia_comparison_inVivo)
DefaultAssay(spermatogonia_comparison_inVivo) <- "RNA"

data.T1.prospg <- as.matrix(GetAssayData(spermatogonia_comparison_inVivo, slot = "data")[, WhichCells(spermatogonia_comparison_inVivo, ident = "T1.prospg")])
data.S0.spg <- as.matrix(GetAssayData(spermatogonia_comparison_inVivo, slot = "data")[, WhichCells(spermatogonia_comparison_inVivo, ident = "S0.spg")])
data.S1.spg <- as.matrix(GetAssayData(spermatogonia_comparison_inVivo, slot = "data")[, WhichCells(spermatogonia_comparison_inVivo, ident = "S1.spg")])
data.differentiating <- as.matrix(GetAssayData(spermatogonia_comparison_inVivo, slot = "data")[, WhichCells(spermatogonia_comparison_inVivo, ident = "Differentiating")])

av.counts.T1.prospg <- apply(data.T1.prospg, 1, mean)
av.counts.S0.spg <- apply(data.S0.spg, 1, mean)
av.counts.S1.spg <- apply(data.S1.spg, 1, mean)
av.counts.differentiating <- apply(data.differentiating, 1, mean)

av.counts.HS_NCG_only_NCG <- cbind(av.counts.T1.prospg, av.counts.S0.spg, av.counts.S1.spg, av.counts.differentiating)


av.counts.HS_NCG_only_NCG.df <- as.data.frame(av.counts.HS_NCG_only_NCG)

# gene.name.list <- rownames(av.counts.HS_NCG_only_NCG.df)
gene.name.list <- rownames(av.counts.HS_NCG_only_NCG.df)

av.counts.HS_NCG_only_NCG.df$gene.name <- gene.name.list

# genes.of.interest <- Keren.gene.list[c(1:48)] #this is just Keren's list of genes
# genes.of.interest <- (all.markers.subset$gene) #this is a log-fold cutoff
# genes.of.interest <- rownames(topNmarkers) #this is just taking the top N markers
# genes.of.interest <- topNmarkers$gene #this is just taking the top N markers
# genes.of.interest <- c("GAS5", "TOMM7", "APOE", "SERPING1", "TPM2", "MEG3", "NANOS2", 
#                        "GAGE12H", "RET", "PIWIL4", "EGR4", "ID4","UTF1",  "SOHLH1", "EOMES", "TSPAN33","KIT", "SOHLH2",
#                        "GFRA1", "L1TD1", "ETV5", "MSL3")
# genes.of.interest <- list.of.genes.all

av.counts.genes.of.interest <- av.counts.HS_NCG_only_NCG.df %>%
  dplyr::filter(gene.name %in% genes.of.interest)


av.counts.genes.of.interest.ordered <- av.counts.genes.of.interest[match(genes.of.interest, av.counts.genes.of.interest$gene.name),]

integer.col <- ncol(av.counts.genes.of.interest.ordered)

av.counts.genes.of.interest.ordered <- av.counts.genes.of.interest.ordered[,-integer.col]

av.counts.genes.of.interest.matrix1 <- as.matrix(av.counts.genes.of.interest.ordered)
colnames(av.counts.genes.of.interest.matrix1) <- c("T1.prospg", "S0.spg", "S1.spg", "Differentiating")

av.counts.genes.of.interest.matrix.no.na <- na.omit(av.counts.genes.of.interest.matrix1)

head(av.counts.genes.of.interest.matrix.no.na)





# dev.off()
l.mat <- rbind(c(0,4), c(0, 3), c(2,1), c(0,0))
l.mat


av.counts.genes.of.interest.matrix.no.na.t <- as.data.frame(t(av.counts.genes.of.interest.matrix.no.na))
av.counts.genes.of.interest.matrix.no.na.t$cell.type <- rownames(av.counts.genes.of.interest.matrix.no.na.t)

av.counts.genes.of.interest.matrix.no.na.t$cell.type <- factor(av.counts.genes.of.interest.matrix.no.na.t$cell.type, levels = c(
  "T1.prospg"  ,     "S0.spg"   ,       "S1.spg"      ,    "Differentiating"
))

av.counts <- av.counts.genes.of.interest.matrix.no.na.t

# Reorder the levels of cell.type according to the specified order
av.counts$cell.type <- factor(av.counts$cell.type, levels = c(
  "T1.prospg",       "S0.spg",          "S1.spg",          "Differentiating"
))



# Create a heatmap with row scaling
heatmap.2(av.counts.genes.of.interest.matrix.no.na, Rowv = NA, Colv = NA, 
          # col = cm.colors(256), 
          # scale_color_viridis(256),
          # col = c("blue", "white", "red"),
          # col = colorRampPalette(c("blue", "white", "red"))(n = 299),
          col = magma,
          trace = "none",
          scale = "row",
          margins = c(5,10), main = "In vivo")




T1.vs.S0 <- FindMarkers(spermatogonia_comparison_inVivo, ident.1 = "T1.prospg", ident.2 = "S0.spg", logfc.threshold = log2fc.cutoff, min.pct = 0)
S0.vs.S1 <- FindMarkers(spermatogonia_comparison_inVivo, ident.1 = "S0.spg", ident.2 = "S1.spg", logfc.threshold = log2fc.cutoff, min.pct = 0)
S1.vs.DiffSpg <- FindMarkers(spermatogonia_comparison_inVivo, ident.1 = "S1.spg", ident.2 = "Differentiating", logfc.threshold = log2fc.cutoff, min.pct = 0)



list.of.comparisons <- list(T1.vs.S0, S0.vs.S1, S1.vs.DiffSpg)
list.of.counts <- list(av.counts.T1.prospg,
                       av.counts.S0.spg,
                       av.counts.S1.spg,
                       av.counts.differentiating)

list.of.cell.types <- c(  "T1.prospg",       "S0.spg",          "S1.spg",          "Differentiating")

list.of.cell.types_abbr <- c(  "T1.prospg",       "S0.spg",          "S1.spg",          "Differentiating")




#CELL TYPE COMPARISON CHARTS
#starting with PGC-E
cell.type.colors <- c("#a6cee2", "#2179b4", "#b4d88b", "#36a047", "#f6999a", "#f15a29", "#fdbf6f", "#f7941d", "#cab3d6", "#6b3f98")
cols2 = c(#"#a6cee2", 
  "#f15a29", "#fdbf6f", "darkgoldenrod3", "#f7941d")

for(number in 1:(length(list.of.cell.types)-1)){
  # number <- 1
  print(number)
  # gene.of.interest <- c(list.of.lists.of.genes[[number]], list.of.lists.of.genes[[number+1]])
  # gene.of.interest <- list.of.lists.of.genes[[number]]
  gene.of.interest <- genes.of.interest
  current.comparison <- list.of.comparisons[[number]]
  current.counts <- list.of.counts[[number]]
  current.counts2 <- list.of.counts[[number+1]]
  
  current.comparison$gene <- rownames(current.comparison)
  # current.comparison.filtered <- filter(current.comparison, p_val_adj < p.val.cutoff)
  current.comparison.filtered.up <- filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC > 0)
  current.comparison.filtered.down <- filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC < 0)
  current.comparison.filtered.significant <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff)
  
  current.comparison.dataframe <- data.frame(current.counts, current.counts2)
  current.comparison.dataframe$log2fc <- current.counts - current.counts2
  current.comparison.dataframe$regulation <- "Nonsig"
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc > log2fc.cutoff] <- "Up"
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc < -log2fc.cutoff] <- "Down"
  print(table(current.comparison.dataframe$regulation))
  
  current.comparison.dataframe$regulation[!(current.comparison.dataframe$gene.name %in% current.comparison.filtered.significant$gene)] <- "Nonsig"
  print(table(current.comparison.dataframe$regulation))
  
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  # current.comparison.subset <- filter(current.comparison.dataframe, av.counts.diffL > 0 | av.counts.Spermatocytes>0)
  current.comparison.subset <- dplyr::filter(current.comparison.dataframe, current.counts > 0 | current.counts2 >0)
  # current.comparison.sig.up <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.up$gene)
  # current.comparison.sig.down <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.down$gene)
  current.comparison.subset <- na.omit(current.comparison.subset)
  current.comparison.subset$DEGs <- "Not significant"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene] <- "Up"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene] <- "Down"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% gene.of.interest] <- "GeneOfInterest"
  current.comparison.subset$regulation[current.comparison.subset$gene.name %in% gene.of.interest] <- "GeneOfInterest"
  print(paste0(list.of.cell.types[number], " vs ", list.of.cell.types[number+1], 
               " Up: ", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]), 
               " Down:", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])))
  
  colnames(current.comparison.subset) <- c("HS.counts", "NCG.counts", "log2FC_manual", "regulation", "gene.name", "DEGs")
  print(table(current.comparison.subset$regulation))
  
  degs_order <- c( "Down", "GeneOfInterest", "Not significant", "Up")
  
  
  current.comparison.subset_ordered <- current.comparison.subset %>%
    mutate(DEGs = factor(DEGs, levels = degs_order)) 
  
  setwd("~/Documents/scRNAseq/Rscripts/savefiles/xrTestis.cell.type.comparisons")
  current.comparison.subset_ordered_noDEG <- current.comparison.subset_ordered[,-6]
  write.csv(current.comparison.subset_ordered_noDEG, row.names = F, file=paste0(list.of.cell.types_abbr[number], "(up)_vs_", list.of.cell.types_abbr[number+1], "_(down)_p.val", p.val.cutoff, "log2fc", log2fc.cutoff, "_inVivo.csv"))
  
  
  
  assign(paste0("r", number), 
         
         ggplot(data = current.comparison.subset_ordered, aes(x = HS.counts, y = NCG.counts,
                                                              label = gene.name,
                                                              color = regulation)) +
           rasterise(geom_point(data = filter(current.comparison.subset_ordered, regulation == "Nonsig"), alpha = 0.5), dpi = 900) + # Plot "Not Significant" points first
           geom_point(data = filter(current.comparison.subset_ordered, regulation != "Nonsig")) + # Plot other points
           # geom_text(data = filter(current.comparison.subset_ordered, regulation != "Nonsig")) +
           geom_point(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest")) +
           geom_text(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest"), hjust = -0.1, vjust = -0.1) +
           
           
           scale_color_manual(values = c(cols2[number+1], "black", '#999999', cols2[number])) +
           theme_classic() +
           xlim(0, 2.5) + ylim(0, 2.5) +
           labs(title = paste(list.of.cell.types[number], " vs. ", list.of.cell.types[number + 1], ", q", p.val.cutoff, ", log2fc", log2fc.cutoff),
                x = list.of.cell.types[number], y = list.of.cell.types[number + 1]) +
           theme(plot.title = element_text(size = 9),
                 panel.grid.major = element_line(color = "gray", linetype = "dotted"),
                 panel.grid.minor = element_line(color = "lightgray", linetype = "dotted")) +
           geom_abline(intercept = 0, slope = 1, linetype = "dotted")
  )
}

plot_grid(ncol = 3, r1, r2, r3) 

### VIOLIN PLOTS ###
### ### ### ### ### 
DefaultAssay(spermatogonia_comparison_inVivo) <- "RNA"
p1 <- Stacked_VlnPlot(spermatogonia_comparison_inVivo, 
                      colors_use = cols2,x_lab_rotate = T,
                      features = c("NANOS2", "EGR4", "PIWIL4", 
                                   "MSL3", "L1TD1", "ETV5", 
                                   "DUSP6", "KIT"))

DefaultAssay(spermatogonia_comparison_xrTestis) <- "RNA"
p2 <- Stacked_VlnPlot(spermatogonia_comparison_xrTestis, 
                      colors_use = cols2,x_lab_rotate = T,
                      features = c("NANOS2", "EGR4", "PIWIL4", 
                                   "MSL3", "L1TD1", "ETV5", 
                                   "DUSP6", "KIT"))
plot_grid(p1, p2)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

### SCATTER PLOT XRTESTIS ### Extended Data Fig 7F (right)

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


Idents(spermatogonia_comparison_xrTestis3) <- "cell.type.spermatogonia"
# DimPlot(spermatogonia_comparison_xrTestis3)
DefaultAssay(spermatogonia_comparison_xrTestis3) <- "RNA"

data.T1.prospg <- as.matrix(GetAssayData(spermatogonia_comparison_xrTestis3, slot = "data")[, WhichCells(spermatogonia_comparison_xrTestis3, ident = "T1.prospg")])
data.S0.spg <- as.matrix(GetAssayData(spermatogonia_comparison_xrTestis3, slot = "data")[, WhichCells(spermatogonia_comparison_xrTestis3, ident = "S0.spg")])
data.S1.spg <- as.matrix(GetAssayData(spermatogonia_comparison_xrTestis3, slot = "data")[, WhichCells(spermatogonia_comparison_xrTestis3, ident = "S1.spg")])
data.differentiating <- as.matrix(GetAssayData(spermatogonia_comparison_xrTestis3, slot = "data")[, WhichCells(spermatogonia_comparison_xrTestis3, ident = "Differentiating")])

av.counts.T1.prospg <- apply(data.T1.prospg, 1, mean)
av.counts.S0.spg <- apply(data.S0.spg, 1, mean)
av.counts.S1.spg <- apply(data.S1.spg, 1, mean)
av.counts.differentiating <- apply(data.differentiating, 1, mean)

av.counts.HS_NCG_only_NCG <- cbind(av.counts.T1.prospg, av.counts.S0.spg, av.counts.S1.spg, av.counts.differentiating)


av.counts.HS_NCG_only_NCG.df <- as.data.frame(av.counts.HS_NCG_only_NCG)

# gene.name.list <- rownames(av.counts.HS_NCG_only_NCG.df)
gene.name.list <- rownames(av.counts.HS_NCG_only_NCG.df)

av.counts.HS_NCG_only_NCG.df$gene.name <- gene.name.list

# genes.of.interest <- Keren.gene.list[c(1:48)] #this is just Keren's list of genes
# genes.of.interest <- (all.markers.subset$gene) #this is a log-fold cutoff
# genes.of.interest <- rownames(topNmarkers) #this is just taking the top N markers
# genes.of.interest <- topNmarkers$gene #this is just taking the top N markers
# genes.of.interest <- c("GAS5", "TOMM7", "APOE", "SERPING1", "TPM2", "MEG3", "NANOS2", 
#                        "GAGE12H", "RET", "PIWIL4", "EGR4", "ID4","UTF1",  "SOHLH1", "EOMES", "TSPAN33","KIT", "SOHLH2",
#                        "GFRA1", "L1TD1", "ETV5", "MSL3")

av.counts.genes.of.interest <- av.counts.HS_NCG_only_NCG.df %>%
  dplyr::filter(gene.name %in% genes.of.interest)


av.counts.genes.of.interest.ordered <- av.counts.genes.of.interest[match(genes.of.interest, av.counts.genes.of.interest$gene.name),]

integer.col <- ncol(av.counts.genes.of.interest.ordered)

av.counts.genes.of.interest.ordered <- av.counts.genes.of.interest.ordered[,-integer.col]

av.counts.genes.of.interest.matrix1 <- as.matrix(av.counts.genes.of.interest.ordered)
colnames(av.counts.genes.of.interest.matrix1) <- c("T1.prospg", "S0.spg", "S1.spg", "Differentiating")

av.counts.genes.of.interest.matrix.no.na <- na.omit(av.counts.genes.of.interest.matrix1)

head(av.counts.genes.of.interest.matrix.no.na)




dev.off()
l.mat <- rbind(c(0,4), c(0, 3), c(2,1), c(0,0))
l.mat


av.counts.genes.of.interest.matrix.no.na.t <- as.data.frame(t(av.counts.genes.of.interest.matrix.no.na))
av.counts.genes.of.interest.matrix.no.na.t$cell.type <- rownames(av.counts.genes.of.interest.matrix.no.na.t)

av.counts.genes.of.interest.matrix.no.na.t$cell.type <- factor(av.counts.genes.of.interest.matrix.no.na.t$cell.type, levels = c(
  "T1.prospg"  ,     "S0.spg"   ,       "S1.spg"      ,    "Differentiating"
))

av.counts <- av.counts.genes.of.interest.matrix.no.na.t

# Reorder the levels of cell.type according to the specified order
av.counts$cell.type <- factor(av.counts$cell.type, levels = c(
  "T1.prospg",       "S0.spg",          "S1.spg",          "Differentiating"
))




scaled_data <- scale(av.counts.genes.of.interest.matrix.no.na, center = TRUE, scale = TRUE)

# Create a heatmap with row scaling
heatmap(scaled_data, Rowv = NA, Colv = NA, 
        # col = cm.colors(256), 
        scale = "row",
        margins = c(5,10), main = "xrTestis")


heatmap.2(av.counts.genes.of.interest.matrix.no.na, Rowv = NA, Colv = NA, 
          # col = cm.colors(256), 
          # scale_color_viridis(256),
          # col = c("blue", "white", "red"),
          # col = colorRampPalette(c("blue", "white", "red"))(n = 299),
          col = magma,
          trace = "none",
          scale = "row",
          margins = c(5,10), main = "xrTestis")





DimPlot(spermatogonia_comparison_xrTestis3)

T1.vs.S0 <- FindMarkers(spermatogonia_comparison_xrTestis3, ident.1 = "T1.prospg", ident.2 = "S0.spg", logfc.threshold = log2fc.cutoff, min.pct = 0)
S0.vs.S1 <- FindMarkers(spermatogonia_comparison_xrTestis3, ident.1 = "S0.spg", ident.2 = "S1.spg", logfc.threshold = log2fc.cutoff, min.pct = 0)
S1.vs.DiffSpg <- FindMarkers(spermatogonia_comparison_xrTestis3, ident.1 = "S1.spg", ident.2 = "Differentiating", logfc.threshold = log2fc.cutoff, min.pct = 0)



list.of.comparisons <- list(T1.vs.S0, S0.vs.S1, S1.vs.DiffSpg)
list.of.counts <- list(av.counts.T1.prospg,
                       av.counts.S0.spg,
                       av.counts.S1.spg,
                       av.counts.differentiating)

list.of.cell.types <- c(  "T1.prospg",       "S0.spg",          "S1.spg",          "Differentiating")

list.of.cell.types_abbr <- c(  "T1.prospg",       "S0.spg",          "S1.spg",          "Differentiating")



guo_gene_list <- read.csv("~/Desktop/guo_gene_list.csv")

state0rows <- guo_gene_list[guo_gene_list$Column2 == 2, ]
state1rows <- guo_gene_list[guo_gene_list$Column2 == 3, ]

for(number in 1:(length(list.of.cell.types)-1)){
  # number <- 1
  print(number)
  # gene.of.interest <- c(list.of.lists.of.genes[[number]], list.of.lists.of.genes[[number+1]])
  # gene.of.interest <- list.of.lists.of.genes[[number]]
  gene.of.interest <- genes.of.interest
  current.comparison <- list.of.comparisons[[number]]
  current.counts <- list.of.counts[[number]]
  current.counts2 <- list.of.counts[[number+1]]
  
  current.comparison$gene <- rownames(current.comparison)
  # current.comparison.filtered <- filter(current.comparison, p_val_adj < p.val.cutoff)
  current.comparison.filtered.up <- filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC > 0)
  current.comparison.filtered.down <- filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC < 0)
  current.comparison.filtered.significant <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff)
  
  
  
  
  current.comparison.dataframe <- data.frame(current.counts, current.counts2)
  current.comparison.dataframe$log2fc <- current.counts - current.counts2
  current.comparison.dataframe$regulation <- "Nonsig"
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc > log2fc.cutoff] <- "Up"
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc < -log2fc.cutoff] <- "Down"
  
  current.comparison.dataframe$regulation[!(current.comparison.dataframe$gene.name %in% current.comparison.filtered.significant$gene)] <- "Nonsig"
  
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  # current.comparison.subset <- filter(current.comparison.dataframe, av.counts.diffL > 0 | av.counts.Spermatocytes>0)
  current.comparison.subset <- dplyr::filter(current.comparison.dataframe, current.counts > 0 | current.counts2 >0)
  # current.comparison.sig.up <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.up$gene)
  # current.comparison.sig.down <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.down$gene)
  current.comparison.subset <- na.omit(current.comparison.subset)
  current.comparison.subset$DEGs <- "Not significant"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene] <- "Up"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene] <- "Down"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% gene.of.interest] <- "GeneOfInterest"
  current.comparison.subset$regulation[current.comparison.subset$gene.name %in% gene.of.interest] <- "GeneOfInterest"
  print(paste0(list.of.cell.types[number], " vs ", list.of.cell.types[number+1], 
               " Up: ", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]), 
               " Down:", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])))
  
  colnames(current.comparison.subset) <- c("HS.counts", "NCG.counts", "log2FC_manual", "regulation", "gene.name", "DEGs")
  print(table(current.comparison.subset$regulation))
  
  degs_order <- c( "Down", "GeneOfInterest", "Not significant", "Up")
  
  
  current.comparison.subset_ordered <- current.comparison.subset %>%
    mutate(DEGs = factor(DEGs, levels = degs_order)) 
  
  setwd("~/Documents/scRNAseq/Rscripts/savefiles/xrTestis.cell.type.comparisons")
  current.comparison.subset_ordered_noDEG <- current.comparison.subset_ordered[,-6]
  write.csv(current.comparison.subset_ordered_noDEG, row.names = F, file=paste0(list.of.cell.types_abbr[number], "(up)_vs_", list.of.cell.types_abbr[number+1], "_(down)_p.val", p.val.cutoff, "log2fc", log2fc.cutoff, "_xrTestis.csv"))
  
  
  
  current.comparison.subset_ordered.inVivo <- read.csv(paste0(list.of.cell.types_abbr[number], "(up)_vs_", list.of.cell.types_abbr[number+1], "_(down)_p.val", p.val.cutoff, "log2fc", log2fc.cutoff, "_inVivo.csv"))
  
  current.comparison.subset_ordered.inVivo$regulation2 <- current.comparison.subset_ordered.inVivo$regulation
  
  # current.comparison.subset_ordered$regulationInVivo <- current.comparison.subset_ordered.inVivo$regulation
  
  
  
  # Task 1: Filter rows where gene.name is found in the other data frame
  current.comparison.subset_ordered <- semi_join(current.comparison.subset_ordered, current.comparison.subset_ordered.inVivo, by = "gene.name")
  current.comparison.subset_ordered.inVivo <- semi_join(current.comparison.subset_ordered.inVivo, current.comparison.subset_ordered, by = "gene.name")
  
  # Task 2: Order rows in current.comparison.subset_ordered.inVivo by the other data frame
  current.comparison.subset_ordered.inVivo <- arrange(current.comparison.subset_ordered.inVivo, match(gene.name, current.comparison.subset_ordered$gene.name))
  
  # View the resulting data frames
  head(current.comparison.subset_ordered)
  head(current.comparison.subset_ordered.inVivo)
  
  
  
  current.comparison.subset_ordered$regulation2 <- current.comparison.subset_ordered.inVivo$regulation2
  
  current.comparison.subset_ordered$GUOgenes <- "NotFound"
  head(current.comparison.subset_ordered)
  
  # Assuming state_1_rows contains the genes in the first column
  state_0_genes <- state0rows$Column1
  state_1_genes <- state1rows$Column1
  
  # Update GUOgenes column in current.comparison.subset_ordered
  current.comparison.subset_ordered$GUOgenes[current.comparison.subset_ordered$gene.name %in% state_0_genes] <- "state0gene"
  current.comparison.subset_ordered$GUOgenes[current.comparison.subset_ordered$gene.name %in% state_1_genes] <- "state1gene"
  
  assign(paste0("p", number), #this is just xrTestis as normal
         
         ggplot(data = current.comparison.subset_ordered, aes(x = HS.counts, y = NCG.counts,
                                                              label = gene.name,
                                                              color = regulation)) +
           geom_point(data = filter(current.comparison.subset_ordered, regulation == "Nonsig"), alpha = 0.5) + # Plot "Not Significant" points first
           geom_point(data = filter(current.comparison.subset_ordered, regulation != "Nonsig")) + # Plot other points
           # geom_text(data = filter(current.comparison.subset_ordered, regulation != "Nonsig")) +
           geom_point(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest")) +
           geom_text(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest"), hjust = -0.1, vjust = -0.1) +
           
           
           scale_color_manual(values = c(cols2[number+1], "black", '#999999', cols2[number])) +
           theme_classic() +
           xlim(0, 2.5) + ylim(0, 2.5) +
           labs(title = paste(list.of.cell.types[number], " vs. ", list.of.cell.types[number + 1], ", q", p.val.cutoff, ", log2fc", log2fc.cutoff),
                x = list.of.cell.types[number], y = list.of.cell.types[number + 1]) +
           theme(plot.title = element_text(size = 9),
                 panel.grid.major = element_line(color = "gray", linetype = "dotted"),
                 panel.grid.minor = element_line(color = "lightgray", linetype = "dotted")) +
           geom_abline(intercept = 0, slope = 1, linetype = "dotted")
  )
  
  assign(paste0("s", number), #this adds in the significant genes from human
         
         ggplot(data = current.comparison.subset_ordered, aes(x = HS.counts, y = NCG.counts,
                                                              label = gene.name,
                                                              color = regulation2)) +
           geom_point(data = filter(current.comparison.subset_ordered, regulation2 == "Nonsig"), alpha = 0.5) + # Plot "Not Significant" points first
           geom_point(data = filter(current.comparison.subset_ordered, regulation2 != "Nonsig")) + # Plot other points
           # geom_text(data = filter(current.comparison.subset_ordered, regulation2 != "Nonsig")) +
           geom_point(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest")) +
           geom_text(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest"), hjust = -0.1, vjust = -0.1) +
           
           
           scale_color_manual(values = c(cols2[number+1], "black", '#999999', cols2[number])) +
           theme_classic() +
           xlim(0, 2.5) + ylim(0, 2.5) +
           labs(title = paste(list.of.cell.types[number], " vs. ", list.of.cell.types[number + 1], ", q", p.val.cutoff, ", log2fc", log2fc.cutoff),
                x = list.of.cell.types[number], y = list.of.cell.types[number + 1]) +
           theme(plot.title = element_text(size = 9),
                 panel.grid.major = element_line(color = "gray", linetype = "dotted"),
                 panel.grid.minor = element_line(color = "lightgray", linetype = "dotted")) +
           geom_abline(intercept = 0, slope = 1, linetype = "dotted")
  )
  
  assign(paste0("t", number), #this paints the S0 and S1 markers from Guo et al 
         
         ggplot(data = current.comparison.subset_ordered, aes(x = HS.counts, y = NCG.counts,
                                                              label = gene.name,
                                                              color = GUOgenes)) +
           geom_point(data = filter(current.comparison.subset_ordered, GUOgenes == "NotFound"), alpha = 0.5) + # Plot "Not Significant" points first
           geom_point(data = filter(current.comparison.subset_ordered, GUOgenes != "NotFound")) + # Plot other points
           # geom_text(data = filter(current.comparison.subset_ordered, GUOgenes != "NotFound")) 
           # geom_point(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest")) +
           # geom_text(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest"), hjust = -0.1, vjust = -0.1) +
           
           
           scale_color_manual(values = c( '#999999', cols2[number+1], cols2[number])) +
           theme_classic() +
           xlim(0, 2.5) + ylim(0, 2.5) +
           labs(title = paste(list.of.cell.types[number], " vs. ", list.of.cell.types[number + 1], ", q", p.val.cutoff, ", log2fc", log2fc.cutoff),
                x = list.of.cell.types[number], y = list.of.cell.types[number + 1]) +
           theme(plot.title = element_text(size = 9),
                 panel.grid.major = element_line(color = "gray", linetype = "dotted"),
                 panel.grid.minor = element_line(color = "lightgray", linetype = "dotted")) +
           geom_abline(intercept = 0, slope = 1, linetype = "dotted")
  )
  
  
}


### Extended Data Fig 7 G ###

# plot_grid(ncol = 3, r1, r2, r3) #in vivo
# plot_grid(ncol = 3, p1, p2, p3) #xrTestis
# plot_grid(ncol = 3, s1, s2, s3) #in vivo significant genes on xrTestis
# plot_grid(ncol = 3, t1, t2, t3) #Guo genes on xrTestis

plot_grid(ncol = 3, r1, r2, r3, 
          p1, p2, p3,
          s1, s2, s3,
          t1, t2, t3
)

plot_grid(ncol = 1,  r2,
          p2, 
          s2
          # t2
)




### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

# compare spermatogonia to spermatogonia LC 




Idents(spermatogonia_comparison2) <- "experiment"
DimPlot(spermatogonia_comparison2)
DimPlot(spermatogonia_comparison2, group.by = "cell.type.spermatogonia", split.by = "experiment")


DimPlot(spermatogonia_comparison_xrTestis3)
DimPlot(spermatogonia_comparison_inVivo)



Idents(spermatogonia_comparison2) <- "cell.type.spermatogonia"

spermatogonia_comparison.T1 <- subset(spermatogonia_comparison2, idents = "T1.prospg")
spermatogonia_comparison.S0 <- subset(spermatogonia_comparison2, idents = "S0.spg")
spermatogonia_comparison.S1 <- subset(spermatogonia_comparison2, idents = "S1.spg")
spermatogonia_comparison.Diff <- subset(spermatogonia_comparison2, idents = "Differentiating")

DimPlot(spermatogonia_comparison2, group.by = "experiment")
Idents(spermatogonia_comparison.T1) <- "experiment"
Idents(spermatogonia_comparison.S0) <- "experiment"
Idents(spermatogonia_comparison.S1) <- "experiment"
Idents(spermatogonia_comparison.Diff) <- "experiment"

T1.markers.hs.vs.ncg <- FindMarkers(spermatogonia_comparison.T1, ident.1 = "human in vivo", ident.2 = "transplant")
S0.markers.hs.vs.ncg <- FindMarkers(spermatogonia_comparison.S0, ident.1 = "human in vivo", ident.2 = "transplant")
S1.markers.hs.vs.ncg <- FindMarkers(spermatogonia_comparison.S1, ident.1 = "human in vivo", ident.2 = "transplant")
Diff.markers.hs.vs.ncg <- FindMarkers(spermatogonia_comparison.Diff, ident.1 = "human in vivo", ident.2 = "transplant")


# setwd("/Users/ewhelan/Desktop/NCG/NCG_marker_genes2")
# write.csv(PGCs.E.markers.hs.vs.ncg, file="HS_NCG_PGCs.E.markers.hs.vs.ncg.csv")
# write.csv(PGCs.L.markers.hs.vs.ncg, file="HS_NCG_PGCs.L.markers.hs.vs.ncg.csv")
# write.csv(M.ProSpermatogonia.markers.hs.vs.ncg, file="HS_NCG_M.ProSpermatogonia.markers.hs.vs.ncg.csv")
# write.csv(M2T1.markers.hs.vs.ncg, file="HS_NCG_M2T1.markers.hs.vs.ncg.csv")



DimPlot(spermatogonia_comparison_xrTestis3)
DimPlot(spermatogonia_comparison_inVivo)

Idents(spermatogonia_comparison_xrTestis3) <- "cell.type.spermatogonia"
Idents(spermatogonia_comparison_inVivo) <- "cell.type.spermatogonia"

data.T1.hs <- as.matrix(GetAssayData(spermatogonia_comparison_inVivo, slot = "data")[, WhichCells(spermatogonia_comparison_inVivo, ident = "T1.prospg")])
data.S0.hs <- as.matrix(GetAssayData(spermatogonia_comparison_inVivo, slot = "data")[, WhichCells(spermatogonia_comparison_inVivo, ident = "S0.spg")])
data.S1.hs <- as.matrix(GetAssayData(spermatogonia_comparison_inVivo, slot = "data")[, WhichCells(spermatogonia_comparison_inVivo, ident = "S1.spg")])
data.Diff.hs <- as.matrix(GetAssayData(spermatogonia_comparison_inVivo, slot = "data")[, WhichCells(spermatogonia_comparison_inVivo, ident = "Differentiating")])

data.T1.ncg <- as.matrix(GetAssayData(spermatogonia_comparison_xrTestis3, slot = "data")[, WhichCells(spermatogonia_comparison_xrTestis3, ident = "T1.prospg")])
data.S0.ncg <- as.matrix(GetAssayData(spermatogonia_comparison_xrTestis3, slot = "data")[, WhichCells(spermatogonia_comparison_xrTestis3, ident = "S0.spg")])
data.S1.ncg <- as.matrix(GetAssayData(spermatogonia_comparison_xrTestis3, slot = "data")[, WhichCells(spermatogonia_comparison_xrTestis3, ident = "S1.spg")])
data.Diff.ncg <- as.matrix(GetAssayData(spermatogonia_comparison_xrTestis3, slot = "data")[, WhichCells(spermatogonia_comparison_xrTestis3, ident = "Differentiating")])

av.counts.T1.hs <- apply(data.T1.hs, 1, mean)
av.counts.S0.hs <- apply(data.S0.hs, 1, mean)
av.counts.S1.hs <- apply(data.S1.hs, 1, mean)
av.counts.Diff.hs <- apply(data.Diff.hs, 1, mean)

av.counts.T1.ncg <- apply(data.T1.ncg, 1, mean)
av.counts.S0.ncg <- apply(data.S0.ncg, 1, mean)
av.counts.S1.ncg <- apply(data.S1.ncg, 1, mean)
av.counts.Diff.ncg <- apply(data.Diff.ncg, 1, mean)


av.counts.hs <- cbind(av.counts.T1.hs, av.counts.S0.hs, av.counts.S1.hs, av.counts.Diff.hs)

av.counts.ncg <- cbind(av.counts.T1.ncg,av.counts.S0.ncg,av.counts.S1.ncg,av.counts.Diff.ncg)


av.counts.hs.df <- as.data.frame(av.counts.hs)
av.counts.ncg.df <- as.data.frame(av.counts.ncg)

gene.name.list.hs <- rownames(av.counts.hs.df)
gene.name.list.ncg <- rownames(av.counts.ncg.df)

length(gene.name.list.hs)
length(gene.name.list.ncg)
length(intersect(gene.name.list.hs, gene.name.list.ncg))

###


###

av.counts.hs.df$gene <- gene.name.list.hs
av.counts.ncg.df$gene <- gene.name.list.ncg


list.of.comparisons <- list(T1.markers.hs.vs.ncg,
                            S0.markers.hs.vs.ncg,
                            S1.markers.hs.vs.ncg,
                            Diff.markers.hs.vs.ncg)

# setwd("~/Desktop/NCG/NCG_vs_HS")
# write.csv(PGCs.E.markers.hs.vs.ncg, file="PGCs.E.markers.hs.vs.ncg.csv")
# write.csv(PGCs.L.markers.hs.vs.ncg, file="PGCs.L.markers.hs.vs.ncg.csv")
# write.csv(M.ProSpermatogonia.markers.hs.vs.ncg, file="M.ProSpermatogonia.markers.hs.vs.ncg.csv")
# write.csv(M2T1.markers.hs.vs.ncg, file="M2T1.markers.hs.vs.ncg.csv")
# write.csv(T.ProSpermatogonia.markers.hs.vs.ncg, file="T.ProSpermatogonia.markers.hs.vs.ncg.csv")
# write.csv(Spg.markers.hs.vs.ncg, file="Spg.markers.hs.vs.ncg.csv")
# write.csv(Diff.Spg.E.markers.hs.vs.ncg, file="Diff.Spg.E.markers.hs.vs.ncg.csv")
# write.csv(Diff.Spg.L.markers.hs.vs.ncg, file="Diff.Spg.L.markers.hs.vs.ncg.csv")
# write.csv(Pre.Leptotene.markers.hs.vs.ncg, file="Pre.Leptotene.markers.hs.vs.ncg.csv")


# Select the columns you want to compare
hs_column <- av.counts.hs.df[, -ncol(av.counts.hs.df)]  # Exclude the "gene" column
ncg_column <- av.counts.ncg.df[, -ncol(av.counts.ncg.df)]  # Exclude the "gene" column

# Compute correlation between corresponding columns
correlation_values <- cor(hs_column, ncg_column, use = "pairwise.complete.obs")

# Show correlation values
print(correlation_values)
rsquared <- correlation_values^2
print(rsquared)

list.of.cell.types <- c("T1 prospermatogonia",
                        "S0 spermatogonia",
                        "S1 spermatogonia",
                        "Differentiating spermatogonia")

list.of.cell.types.short <- c("T1.prospg",
                              "S0.spg",
                              "S1.spg",
                              "Differentiating")


inVivoData <- av.counts.hs.df[,-5]
inVitroData <- av.counts.ncg.df[,-5]
correlation_matrix <- cor(inVivoData, inVitroData)
print(correlation_matrix)
correltation_matrix_r <- correlation_matrix
correlation_matrix <- correlation_matrix^2


### Extended Data Fig 7 ###

# Convert the correlation matrix to a tibble
correlation_df <- as_tibble(as.data.frame(rsquared), rownames = "Cell_Type_1")

# Factor the levels to maintain the original order
correlation_df$Cell_Type_1 <- factor(correlation_df$Cell_Type_1, levels = rownames(correlation_matrix))
correlation_df <- gather(correlation_df, key = "Cell_Type_2", value = "Correlation", -Cell_Type_1)

# Factor the levels for Cell_Type_2
correlation_df$Cell_Type_2 <- factor(correlation_df$Cell_Type_2, levels = colnames(correlation_matrix))

# Create the heatmap using ggplot2
heatmap_plot <- ggplot(correlation_df, aes(x = Cell_Type_1, y = Cell_Type_2, fill = Correlation)) +
  geom_tile(color = "white", size = 0.5) +
  geom_text(aes(label = round(Correlation, 2)), vjust = 1) +  # Add correlation values
  # scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.8, na.value = NA) +
  # scale_fill_gradient2(low = "#451d5b", mid = "white", high = "red", midpoint = 0.8, na.value = NA) +
  scale_fill_viridis(begin = 0.5) +
  theme_minimal() +
  labs(title = "Correlation Heatmap",
       x = "Cell Type (In Vivo)",
       y = "Cell Type (In Vitro)") +
  theme(axis.text.x = element_text(angle = -22, hjust = 0))  # Rotate x-axis labels at a 45-degree angle

# Print the heatmap
print(heatmap_plot)



p.val.cutoff <- 0.05

log2fc.cutoff <- 0.3

viridis.colours <- inferno(5)

genes.of.interest <- c("SMC1B", "VCX", "VCX2", "SMC3", "TEX101", "SYCE",
                       "IL13RA2", "MDK", "TUBB", "RPS5", "ID4", "RB1")

gene.of.interest.PGCE <- c("TFAP2C", "SOX17", "YBX1", "MIF", "POU5F1", 
                           "NANOS3", "TOP2A", "KHDC3L", "MKI67", "DND1", "DAZL")
gene.of.interest.PGCL1 <- c("DDX4", "DAZL", "MKI67",  "TOP2A", "ID4")
gene.of.interest.PGCL2 <- c("ACTG2",  "POU5F1", "S100A13", "CCNB1", "PCSK1N" )
gene.of.interest.M1 <- c("DAZL", "ID4", "DDX43", "DDX4", "E2F8")
gene.of.interest.M2 <- c("NANOG", "ACTG2", "MKI67", "NANOS3", "POU5F1", "STMN1", "CENPA")
gene.of.interest.transitional1 <- c("MDK", "ASB9", "CITED2", "MAP4K4", "TDRD1", "ID3")
gene.of.interest.transitional2 <- c("MDK", "ASB9", "NANOS3", "POU5F1", "TSPAN6", "APOE", "STMN1", "MEG3")
gene.of.interest.T1 <- c("NANOS2", "RHOXF1", "TDRD1", "MAGEB2", "TEX15", "DCAF4L1", "MAGEC2", "PAGE2B",
                         "SIX1", "DDX4", "FOXP1", "TCF3", "VCX", 'EGR4', "VCX3B", "DPPA2")
gene.of.interest.T1_2 <- c("MEG3", "NANOS2", "TPM2")
gene.of.interest.spg1 <- c("PIWIL4", "SIX1", "HES5", "MAGEB2",  "TCF3", "UTF1", "EGR4")
gene.of.interest.spg2 <- c("DPPA2", "SIX1", "FGFR3", "PIWIL4", "EGR4", "ID4", "SOHLH1", "UTF1")
gene.of.interest.earlydiff1 <- c("ASB9", "STRA8","PRSS21", "SYCP3", "MLEC", "PARP4", "REC8", "LARP1", "DNMT1", "PIWIL1", "TOP2B")
gene.of.interest.earlydiff2 <- c("KIT", "PIWIL1", "ASB9", "ETV7", "L1TD1", "SCML1", "MKI67", "TEX15", "HORMAD1", "MEIOC", "SYCP3")
gene.of.interest.latediff1 <- c("SYCP3", "SMC1B", "MEIOC", "MAGEA4", "UBE2C", "ESX1", "MKI67")
gene.of.interest.latediff2 <- c("CENPF", "KIT", "ITGA6", "UBE2C", "ELAVL2",  "RHOXF1")
gene.of.interest.prelep1 <- c("TEX101", "MEIOB", "HORMAD1", "SYCP2", "SYCP1", "ZCWPW1","SCML1", "TOP2A","SYCP3", "SPO11", "PRDM9")


genes.of.interest.PGCE <- c("TFAP2C", "SOX17", "YBX1", "POU5F1", "PRDM1", 
                            "KIT",   "KRT8", "KRT18",  
                            "PCSK1N", #up 
                            "CENPF", "SCG5", 
                            "TOP2A", "DPPA4", "NANOS3", "TOP2A", "KHDC3L",  "DND1") #down
genes.of.interest.PGCL <- c("TFAP2C", "SOX17", "GSTP1", "MIF", "POU5F1", "DDX4", "DAZL", "MKI67",  "TOP2A",
                            "KRT18",  
                            "VIM", "IFI16", "PDPN", "TMED9", "CHEK1", "CHCHD2", 
                            "NANOS3", "TOP2A", "KHDC3L", "MKI67", "DND1", "DAZL") #down

# genes.of.interest.PGCL <- c("TFAP2C", "SOX17", "YBX1", "MIF", "POU5F1", "DDX4", "DAZL", "MKI67",  "TOP2A", "ID4")
genes.of.interest.M <- c("POU5F1", "SOX17", "TFAP2C", "PLCG2", "CRABP1",
                         "DDX4", "DAZL", "MKI67",  "TOP2A", "ID4", "NANOG", "ACTG2", "MKI67", "NANOS3", "PRDX5", 
                         "SIVA1", "MT-ND1", "MT-ND2", "XIST", "PHLDA3", "TMED9", "PDPN")
genes.of.interest.intermediate <- c("PRXL2A", "TLE5", "POU5F1", "SOX17", "TFAP2C", 
                                    "RHOXF1", "PAGEB2", "NANOG", "ID1", "DCAF4L1", "MAGEC2", 
                                    "MKI67", "THY1", "CHEK1", "HEBP2", "SKAP2", 
                                    "CITED2", "MAP4K4", "TDRD1")
genes.of.interest.T1 <- c("NANOS2", "DCAF4L1", "MAGEC2", "PAGE2B", "MT-ND4L", "PRAME", 
                          "DDX4",  "TCF3", "VCX", 'EGR4', "PLPPR3",   "GAGE12H", "APOE", "RPL23", 
                          "XIST",  "ID2", "ID3")
genes.of.interest.spermatogonia.MAIN <- c("FGFR3", "SIX1", "MAGEB2", 
                                          "ID2", "ID3", "ETV7", "ACTB", "ID4", "DDX4", "MORC1", "PIWIL4")

genes.of.interest.spermatogonia <- c("LY6K", "GAGE12H", "JUND", "PRDX5", "ZNF518A",  "DHRS2",
                                     # "DPPA2", , "CTAG2"
                                     # "SIX1", "FGFR3", 
                                     "EGR4",  "SOHLH1", "UTF1", 
                                     "MT-ND5", 
                                     # "MAGEB2", 
                                     "TEX15")

genes.of.interest.earlydiff <- c("KIT","DAZL", "PAGE5", "PAGE2", "ID1", "JUND", "CYP26A1", 
                                 # "CYP26B1", 
                                 "RPL41", "RPL10A", "ZNF428", 
                                 "ASB9", "SYCP3", "SOHLH2", "L1TD1", "MKI67"
)
genes.of.interest.earlydiff.MAIN <- c("STRA8", "PIWIL1", "REC8", "TEX11", "ACTB", 
                                      "SOHLH1", "UTF1", "ID4", "MAGEA4", "DDX4")


genes.of.interest.latediff <- c("SYCP3", "SMC1B", "MEIOC", "MAGEA4", "UBE2C", "MKI67", "TOP2A", 
                                "VCX", "RPL10", "LDHB", "RACK1", "TUBB", "MDK",  "PCNA", 
                                "CENPF", "KIT", "UBE2C", "ELAVL2")
genes.of.interest.prelep <- c("TEX101", "MEIOB", "HORMAD1", "SYCP2", "SYCP1", "ZCWPW1","SCML1", "TOP2A","SYCP3", "SPO11", "PRDM9", 
                              "MEIOC", "TEX19", "RAD51AP2",  
                              "RPS12", "RPL10", "RPS5", "MDK", "TUBB")

gene.list.premigratory <- c("PCDH1", "ID1", "MSX2", "KAL1", "SMAD7", "KRT19",
                            "CYFIP2", "SLC16A2", "EPHB1", "ASS1", "PCSK1N", 
                            "CABLES1", "ADAM19", "WFDC2", "SH3BGRL3", "C1orf85",
                            "LGALS3", "TRO", "FAM54A", "CDT1", "POLA1", "GDF3", "PMEL", "GREM1", 
                            "RIN2", "PPP1R14C", "ENPP2", "EPAS1", "WSCD1", "MYH2", "HAPLN3", 
                            "TPH2", "LONRF3", "GLDC", "SCG5", "SLC7A3", "ATP6V1H", 
                            "TFRAF3IP2", "CTSL2", "T", "SC5DL", "GSN")

# genes.of.interest.T1 <- c("GAGE12H", "APOE", "RPL23", "RPS26")
# genes.of.interest.spermatogonia1 <- c("MT-ND5", "ID3", "ID1", "ID2", "MAGEB2", "TEX15", "CTAG2")
# genes.of.interest.spermatogonia2 <- c("ID2", "ID3", "ETV7",  "ID4", 
#                        "MORC1", "PIWIL4", "ACTB", "DDX4")
# gene.of.interest.earlydiff <- c("REC8","TEX11", "ACTB", "PIWIL1", "STRA8",
#                                 "SOHLH1", "UTF1", "ID4", "MAGEA4", "DDX4")
list.of.lists.of.genes <- list(genes.of.interest.PGCE,
                               genes.of.interest.PGCL,
                               genes.of.interest.M,
                               genes.of.interest.intermediate,
                               genes.of.interest.T1,
                               genes.of.interest.spermatogonia,
                               genes.of.interest.earlydiff,
                               genes.of.interest.latediff,
                               genes.of.interest.prelep
)

length(list.of.comparisons)

for(number in 1:(length(list.of.comparisons))){
  # number <- 5
  print(number)
  # gene.of.interest <- c(list.of.lists.of.genes[[number]])
  # gene.of.interest <- c("FGFR2", "FGFR3")
  gene.of.interest <- list.of.lists.of.genes[[number]]
  # gene.of.interest <- gene.list.premigratory
  current.comparison <- list.of.comparisons[[number]]
  
  # current.counts <- av.counts.hs.df[[number]]
  # current.counts2 <- av.counts.ncg.df[[number]]
  
  current.counts <- av.counts.hs.df[,number]
  current.counts2 <- av.counts.ncg.df[,number]
  
  current.comparison$gene <- rownames(current.comparison)
  # current.comparison.filtered <- filter(current.comparison, p_val_adj < p.val.cutoff)
  current.comparison.filtered.up <- filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC > 0)
  current.comparison.filtered.down <- filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC < 0)
  
  all.significant.genes <- filter(current.comparison, p_val_adj < p.val.cutoff)
  all.significant.genes <- all.significant.genes$gene
  
  current.comparison.dataframe <- data.frame(current.counts, current.counts2)
  current.comparison.dataframe$log2fc <- current.counts - current.counts2
  current.comparison.dataframe$regulation <- "notDifferent"
  # current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  current.comparison.dataframe$gene.name <- rownames(av.counts.hs.df)
  rownames(current.comparison.dataframe) <- rownames(av.counts.hs.df)
  
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc > log2fc.cutoff] <- "higherInVivo"
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc < -log2fc.cutoff] <- "higherInVitro"
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  
  all.genes <- rownames(current.comparison)
  nonsignificant.genes <- all.genes[!(all.genes %in% all.significant.genes)]
  
  
  # current.comparison.dataframe$regulation[!(current.comparison.dataframe$gene.name %in% all.significant.genes)] <- "notDifferent"
  
  
  # current.comparison.subset <- filter(current.comparison.dataframe, av.counts.diffL > 0 | av.counts.Spermatocytes>0)
  current.comparison.subset <- dplyr::filter(current.comparison.dataframe, current.counts > 0 | current.counts2 >0)
  # current.comparison.sig.up <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.up$gene)
  # current.comparison.sig.down <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.down$gene)
  current.comparison.subset <- na.omit(current.comparison.subset)
  current.comparison.subset$DEGs <- "Not significant"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene] <- "higherInVivo"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene] <- "higherInVitro"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% gene.of.interest] <- "GeneOfInterest"
  current.comparison.subset$regulation[current.comparison.subset$gene.name %in% gene.of.interest] <- "GeneOfInterest"
  print(paste0(list.of.cell.types[number], ": xrTestis vs in Vivo", 
               " Up: ", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]), 
               " Down:", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])))
  
  colnames(current.comparison.subset) <- c("HS.counts", "NCG.counts", "log2FC_manual", "regulation", "gene.name", "DEGs")
  
  print(list.of.cell.types[number])
  print(table(current.comparison.subset$regulation))
  
  
  degs_order <- c( "higherInVitro", "GeneOfInterest", "Not significant", "higherInVivo")
  
  
  current.comparison.subset_ordered <- current.comparison.subset %>%
    mutate(DEGs = factor(DEGs, levels = degs_order)) 
  
  # setwd("~/Desktop/xrTestisFigures/inVivo.vs.xrTestis")
  write.csv(current.comparison.subset_ordered, file = paste0(list.of.cell.types.short[number], "_comparison_UPinVivo_DOWNxrTestis_",
                                                             log2fc.cutoff, "_log2FC.csv"))
  table(current.comparison.subset_ordered$regulation)
  
  assign(paste0("r", number), 
         ggplot(data = current.comparison.subset_ordered, aes(x = HS.counts, y = NCG.counts,
                                                              label = gene.name,
                                                              color = regulation)) +
           geom_point(data = filter(current.comparison.subset_ordered, regulation == "notDifferent"), alpha = 0.5) + # Plot "Not Significant" points first
           geom_point(data = filter(current.comparison.subset_ordered, regulation != "notDifferent")) + # Plot other points
           # geom_text(data = filter(current.comparison.subset_ordered, regulation != "notDifferent"),hjust=-0.1, vjust=-0.1) +
           geom_point(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest")) +
           geom_text(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest"), hjust = -0.1, vjust = -0.1) +
           
           
           # scale_color_manual(values = c(cell.type.colors[number+1], "red4", '#999999', cell.type.colors[number])) +
           scale_color_manual(values = c("black","blue",   "red3",'#999999')) +
           theme_classic() +
           xlim(0, 3.5) + ylim(0, 3.5) + 
           labs(title = paste(list.of.cell.types[number], " xrTestis vs in Vivo", ", q", p.val.cutoff, ", log2fc", log2fc.cutoff),
                x = paste0(list.of.cell.types[number], " in Vivo"), y = paste0(list.of.cell.types[number], " xrTestis")) +
           theme(plot.title = element_text(size = 9),
                 panel.grid.major = element_line(color = "gray", linetype = "dotted"),
                 panel.grid.minor = element_line(color = "lightgray", linetype = "dotted")) +
           geom_abline(intercept = 0, slope = 1, linetype = "dotted")
  )
  

  
}

plot_grid(ncol = 4, r1, r2, r3, r4)


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### PART 3          ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### NANOS3 DND1     ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###


###############################################################################################################################
#Read in data


setwd("~/Documents/PROJECTS/ksasaki/NANOS3KO")

load("NANOS3_project/NANOS3_part1_20210809.rdata")
nanos3_raw <- readRDS("NANOS3_project/big_NANOS3_raw.seurat.object.rds")
load("NANOS3_bulk/NANOS3_bulk.rdata")

load("DND1_project/DND1_part1_20210809.rdata")
dnd1_raw <- readRDS("DND1_project/big_DND1_raw.seurat.object.rds")


DimPlot(big_DND1)
DimPlot(big_NANOS3)
FeaturePlot(big_DND1, features = "TFAP2C")
FeaturePlot(big_NANOS3, features = "NANOS3", split.by = "sample")
VlnPlot(big_NANOS3, features = "NANOS3")
d81_xrTestis <- Read10X(data.dir = "NatComm2020_xrTestis/d124_xrTestis/")
d124_xrTestis <- Read10X(data.dir = "NatComm2020_xrTestis/d81_xrTestis/")

###############################################################################################################################
# Merge datasets

names.data <- c("C1_191121", "C3_191203", "C1_191213", "d81_xrTestis", "d124_xrTestis")
human.list <- c(C1_191121_Hu, C3_191203_hu, C1_191213_hu, d81_xrTestis, d124_xrTestis)


for(a in 1:length(human.list)){
  print(length(Cells(human.list[[a]])))
}

inflection.point1 <- c(1:length(human.list))
log_lib_size_at_inflection1 <- c(1:length(human.list))

par(mfrow = c(3,3))

for(y in 1:length(human.list)){
  #y=1
  print(names.data[y])
  expression.df <- as.data.frame(human.list[y])
  
  umi_per_barcode <- colSums(expression.df)
  barcode_rank <- rank(-umi_per_barcode)
  # plot(barcode_rank, umi_per_barcode,
  #      #ylim=c(1,2000),
  #      xlim=c(1,10000))
  
  log_lib_size <- log10(umi_per_barcode)
  #plot(barcode_rank, log_lib_size, xlim=c(1,10000))
  
  o <- order(barcode_rank)
  log_lib_size <- log_lib_size[o]
  barcode_rank <- barcode_rank[o]
  
  rawdiff <- diff(log_lib_size)/diff(barcode_rank)
  inflection <- which(rawdiff == min(rawdiff[1:4000], na.rm=TRUE))
  
  
  
  plot(barcode_rank, log_lib_size, xlim=c(1,8000),
       pch = 20,
       main = paste(y, names.data[y], "cells", inflection, sep="_", "UMIs", round(10^log_lib_size[inflection]))
  )
  
  abline(v=inflection, col="blue", lwd=2)
  abline(h=log_lib_size[inflection], col="green", lwd=2)
  
  
  inflection.point1[y] <- inflection
  log_lib_size_at_inflection1[y] <- log_lib_size[inflection]
}

min.UMIs1 <- 10^log_lib_size_at_inflection1 #something is wrong with this, just make it manually

min.UMIs1 <- vector("numeric" , length(human.list)+2)



min.UMIs1[1] <- 10^4
min.UMIs1[2] <- 10^3.8
min.UMIs1[3] <- 10^4
min.UMIs1[4] <- 10^4
min.UMIs1[5] <- 10^3.5



min.features <- 100
max.features <- 3000
max.mito <- 20



C1_191121 <- CreateSeuratObject(C1_191121_Hu, project = "C1_191121")
C1_191121[["orig.ident"]] <- "C1_191121"
C3_191203 <- CreateSeuratObject(C3_191203_hu, project = "C3_191203")
C3_191203[["orig.ident"]] <- "C3_191203"
C1_191213 <- CreateSeuratObject(C1_191213_hu, project = "C1_191213")
C1_191213[["orig.ident"]] <- "C1_191213"

d81_xrTestis.seurat <- CreateSeuratObject(d81_xrTestis, project = "d81_xrTestis")
d81_xrTestis.seurat[["orig.ident"]] <- "d81_xrTestis"
d124_xrTestis.seurat <- CreateSeuratObject(d124_xrTestis, project = "d124_xrTestis")
d124_xrTestis.seurat[["orig.ident"]] <- "d124_xrTestis"

DefaultAssay(big_DND1) <- "SCT"
DefaultAssay(big_NANOS3) <- "SCT"

C1_191121 <- SCTransform(C1_191121, verbose = FALSE)
C3_191203 <- SCTransform(C3_191203, verbose = FALSE)
C1_191213 <- SCTransform(C1_191213,verbose = FALSE)
d81_xrTestis.seurat <- SCTransform(d81_xrTestis.seurat, verbose = FALSE)
d124_xrTestis.seurat <- SCTransform(d124_xrTestis.seurat, verbose = FALSE)

big_DND1_WT <- subset(big_DND1, idents = "DND1_WT_d14")
big_NANOS3_WT <- subset(big_NANOS3, idents = "NANOS3_WT_d9_xrT")

list.of.seurats <- c(
  C1_191121, C3_191203, C1_191213,
  
  
  big_NANOS3_WT,
  big_DND1_WT,
  d81_xrTestis.seurat, d124_xrTestis.seurat
)



experiment <- c("hiPSCs", "hPGCLC", "iMELCs",
                "d9", "d14",  "d81_xrTestis", "d124_xrTestis"
)


names.data <- c(
  "C1_191121", "C3_191203", "C1_191213", 
  
  "NANOS3_KO_WT",
  "DND1_KO_WT",
  "d81_xrTestis", "d124_xrTestis"
)

for(current.sample in 1:length(list.of.seurats)){
  # current.sample <- 1
  print(current.sample)
  current.seurat <- list.of.seurats[[current.sample]]
  
  current.seurat$replicate <- names.data[current.sample]
  current.seurat$experiment <- experiment[current.sample]
  current.seurat[['percent.mito']] <- PercentageFeatureSet(current.seurat, pattern = "^MT-")
  list.of.seurats[[current.sample]] <- current.seurat
}

for(current.sample in 1:(length(list.of.seurats))){
  # current.sample <- 1
  current.seurat <- list.of.seurats[[current.sample]]
  current.seurat <- subset(x = current.seurat,
                           # subset = nFeature_RNA > min.features & nFeature_RNA < max.features & percent.mito < max.mito & nCount_RNA > min.UMIs[current.sample])
                           subset = percent.mito < max.mito & nCount_RNA > min.UMIs1[current.sample]) #
  # current.seurat <- subset(x = current.seurat, subset = species == "human")
  list.of.seurats[[current.sample]] <- current.seurat
}



number.of.samples <- length(list.of.seurats)
number.of.cells.per.sample <-  data.frame(
  SampleName = c(1:number.of.samples),
  CellNumber = c(1:number.of.samples),
  GeneNumber = c(1:number.of.samples),
  SampleNumber = c(1:number.of.samples),
  # SampleOrigin = c(1:number.of.samples),
  medianUMI = c(1:number.of.samples),
  medianGenes = c(1:number.of.samples),
  medianMito =c(1:number.of.samples),
  cutoff = c(1:number.of.samples)
)

for(a in 1:length(list.of.seurats)){
  number.of.cells.per.sample[a,4]<-a
  number.of.cells.per.sample[a,1]<-names.data[a]
  number.of.cells.per.sample[a,2]<-length(Cells(list.of.seurats[[a]]))
  number.of.cells.per.sample[a,3]<-length(rownames(list.of.seurats[[a]]))
  # number.of.cells.per.sample[a,5]<-stage.data[a]
  number.of.cells.per.sample[a,5]<-median(list.of.seurats[[a]]$nCount_RNA)
  number.of.cells.per.sample[a,6]<-median(list.of.seurats[[a]]$nFeature_RNA)
  number.of.cells.per.sample[a,7]<-median(list.of.seurats[[a]]$percent.mito)
  number.of.cells.per.sample[a,8]<-min.UMIs1[a]
  #print(number.of.cells.per.sample[a,])
}

print(number.of.cells.per.sample)

current.list <- list.of.seurats

for (i in 1:length(x = current.list)) {
  current.list[[i]] <- NormalizeData(object = current.list[[i]], verbose = FALSE)
  
  current.list[[i]] <- FindVariableFeatures(object = current.list[[i]], nfeatures = 2000, verbose = FALSE)
}

current.anchors <- FindIntegrationAnchors(object.list = current.list, dims = 1:30)
# save(current.anchors, file="current.anchors.Robj")
current.integrated <- IntegrateData(anchorset = current.anchors,#normalization.method = c("LogNormalize"), 
                                    dims = 1:30)
# save(current.integrated, file="current.integrated.Robj")

DefaultAssay(object = current.integrated) <- "integrated"
current.integrated <- ScaleData(object = current.integrated, vars.to.regress = c("nCount_RNA", "percent.mito"), verbose = FALSE)
current.integrated <- RunPCA(object = current.integrated,
                             npcs = 30,
                             verbose = TRUE)
current.integrated <- FindNeighbors(object = current.integrated)
current.integrated <- FindClusters(object = current.integrated,
                                   resolution = 1.1,
                                   #algorithm = 1,
                                   verbose = TRUE
)
# current.integrated <- RunTSNE(object = current.integrated, dims = 1:20)
current.integrated <- RunUMAP(object = current.integrated, reduction = "pca",
                              dims = 1:30)

DimPlot(object = current.integrated, reduction = "umap", group.by = "experiment", label = T)
DefaultAssay(current.integrated) <- "SCT"
FeaturePlot(current.integrated, features = c("TFAP2C", "DDX4"))
FeaturePlot(current.integrated, split.by = "experiment", 
            features = c("TCL1A"))

current.integrated$experiment <- factor(x = current.integrated$experiment, levels = c("hiPSCs",
                                                                                      "iMELCs",
                                                                                      "hPGCLC",
                                                                                      "d9",
                                                                                      "d14",
                                                                                      "d81_xrTestis",
                                                                                      "d124_xrTestis"
))

gene.of.interest <- "DND1"
DefaultAssay(current.integrated) <- "SCT"
p1 <- VlnPlot(current.integrated, features = gene.of.interest, group.by = "experiment")
p2 <- FeaturePlot(current.integrated, split.by = "experiment", features = gene.of.interest)

plot_grid(p1, p2, ncol = 1)

current.integrated -> ALI
DimPlot(current.integrated)
Idents(current.integrated) <- "experiment"
table(Idents(current.integrated))
ALI_subset <- subset(current.integrated, idents = c("d81_xrTestis", "d124_xrTestis"), invert = T)

p1 <- VlnPlot(ALI_subset, features = "DND1", group.by = "experiment")
p2 <- VlnPlot(ALI_subset, features = "NANOS3", group.by = "experiment")
plot_grid(p1, p2, ncol = 1)

p1 <- VlnPlot(ALI_subset, features = "DND1", group.by = "experiment", pt.size = 0)
p2 <- VlnPlot(ALI_subset, features = "NANOS3", group.by = "experiment", pt.size = 0)
plot_grid(p1, p2, ncol = 1)



###############################################################################################################################
### Scatter plot comparing hiPSCs to d9 or d14




DimPlot(ALI)
Idents(ALI) <- ALI$experiment


data.pgclc <- as.matrix(GetAssayData(ALI, slot = "data")[, WhichCells(ALI, ident = "hPGCLC")])
data.d9 <- as.matrix(GetAssayData(ALI, slot = "data")[, WhichCells(ALI, ident = "d9")])
data.d14 <- as.matrix(GetAssayData(ALI, slot = "data")[, WhichCells(ALI, ident = "d14")])
data.d124 <- as.matrix(GetAssayData(ALI, slot = "data")[, WhichCells(ALI, ident = "d124_xrTestis")])

av.counts.pgclc <- apply(data.pgclc, 1, mean)
av.counts.d9 <- apply(data.d9, 1, mean)
av.counts.d14 <- apply(data.d14, 1, mean)
av.counts.d124 <- apply(data.d124, 1, mean)

# av.counts.HS_NCG_only_NCG <- cbind(av.counts.pgclc, av.counts.d9, av.counts.d14)
# 
# 
# av.counts.HS_NCG_only_NCG.df <- as.data.frame(av.counts.HS_NCG_only_NCG)
# 
# gene.name.list <- rownames(av.counts.HS_NCG_only_NCG.df)
# 
# av.counts.HS_NCG_only_NCG.df$gene.name <- gene.name.list
# 


###
# ALI <- NCG_germ4
p.val.cutoff <- 0.05
log2fc.cutoff <- 0.25

###
#banana
# ALI_original -> ALI

# ALI <- ALI

DimPlot(ALI)

DefaultAssay(ALI) <- "SCT"
slot(ALI$SCT@SCTModel.list[[1]], 'median_umi') = median(ALI$SCT@SCTModel.list[[1]]@cell.attributes$umi)

ALI <- SCTransform(ALI)
ALI <- PrepSCTFindMarkers(ALI)
PGCLC.vs.d9 <- FindMarkers(ALI, ident.1 = "hPGCLC", ident.2 = "d9", logfc.threshold = log2fc.cutoff, min.pct = 0)
PGCLC.vs.d14 <- FindMarkers(ALI, ident.1 = "hPGCLC", ident.2 = "d14", logfc.threshold = log2fc.cutoff, min.pct = 0)
PGCLC.vs.d124 <- FindMarkers(ALI, ident.1 = "hPGCLC", ident.2 = "d124_xrTestis", logfc.threshold = log2fc.cutoff, min.pct = 0)


list.of.comparisons <- list(PGCLC.vs.d9, PGCLC.vs.d14, PGCLC.vs.d124)
list.of.counts <- list(av.counts.pgclc, av.counts.d9, av.counts.d14, av.counts.d124)

list.of.cell.types <- c("hPGCLC", 
                        "d9", "d14", "d124")



for(number in 1:(length(list.of.cell.types)-1)){
  # number <- 1
  
  current.comparison <- list.of.comparisons[[number]]
  
  current.counts <- list.of.counts[[1]] #note this is always comparing hPGCLCs
  # current.counts <- (log2(exp((current.counts))))
  current.counts2 <- list.of.counts[[number+1]]
  # current.counts2 <- (log2(exp((current.counts2))))
  current.comparison$gene <- rownames(current.comparison)
  # current.comparison.filtered <- filter(current.comparison, p_val_adj < p.val.cutoff)
  
  
  current.comparison.filtered.up <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC > log2fc.cutoff)
  current.comparison.filtered.down <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC < -log2fc.cutoff)
  current.comparison.filtered.significant <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff)
  
  current.comparison.dataframe <- data.frame(current.counts, current.counts2)
  # current.counts.df <- data.frame(current.counts, current.counts2)
  current.comparison.dataframe$log2fc <- current.counts - current.counts2
  current.comparison.dataframe$regulation <- "not"
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc > log2fc.cutoff] <- "up"
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc < -log2fc.cutoff] <- "down"
  current.comparison.dataframe$regulation[!(current.comparison.dataframe$gene.name %in% current.comparison.filtered.significant$gene)] <- "not"
  
  
  # #simply deduct count values
  # current.comparison.dataframe$difference <- current.comparison.dataframe$current.counts - current.comparison.dataframe$current.counts2
  # 
  
  
  
  current.comparison.subset <- dplyr::filter(current.comparison.dataframe, current.counts > 0 | current.counts2 >0)
  # current.comparison.sig.up <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.up$gene)
  # current.comparison.sig.down <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.down$gene)
  current.comparison.subset <- na.omit(current.comparison.subset)
  current.comparison.subset$DEGs <- "Not significant"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene] <- "Up"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene] <- "Down"
  print(paste0(list.of.cell.types[number], " vs ", list.of.cell.types[number+1], 
               " Up: ", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]), 
               " Down:", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])))
  
  print(table(current.comparison.subset$regulation))
  
  
  assign(paste0("p", number),
         (ggplot() + geom_point(data = current.comparison.subset, aes(x = current.counts, y = current.counts2,  color = regulation)) +
            scale_color_manual(values=c('#E69F00','#999999', '#56B4E9')) +
            theme_classic() + labs(title=paste(list.of.cell.types[number], " vs. ", list.of.cell.types[number+1], ", q", p.val.cutoff, ", log2fc", log2fc.cutoff),
                                   x=list.of.cell.types[number], y = list.of.cell.types[number+1]) + theme(plot.title = element_text(size=9)))
  )
  
}

p1
p2
p3


cell.type.colors <- c("#2179b4", "#b4d88b", "#36a047", "#f6999a", "#f15a29", "#fdbf6f", "#f7941d", "#cab3d6", "#6b3f98")


gene.of.interest.PGCE <- c("TFAP2C", "SOX17", "YBX1", "MIF", "POU5F1")
gene.of.interest.PGCL1 <- c("DDX4", "DAZL", "MKI67",  "TOP2A", "ID4")
gene.of.interest.PGCL2 <- c("ACTG2",  "POU5F1", "S100A13", "CCNB1", "PCSK1N" )
gene.of.interest.M1 <- c("DAZL", "ID4", "DDX43", "DDX4", "E2F8")
gene.of.interest.M2 <- c("NANOG", "ACTG2", "MKI67", "NANOS3", "POU5F1", "STMN1", "CENPA")
gene.of.interest.transitional1 <- c("MDK", "ASB9", "CITED2", "MAP4K4", "TDRD1", "ID3")
gene.of.interest.transitional2 <- c("MDK", "ASB9", "NANOS3", "POU5F1", "TSPAN6", "APOE", "STMN1", "MEG3")
gene.of.interest.T1 <- c("NANOS2", "RHOXF1", "TDRD1", "MAGEB2", "TEX15", "DCAF4L1", "MAGEC2", "PAGE2B",
                         "SIX1", "DDX4", "FOXP1", "TCF3", "VCX", 'EGR4', "VCX3B", "DPPA2")
gene.of.interest.T1_2 <- c("MEG3", "NANOS2", "TPM2")
gene.of.interest.spg1 <- c("PIWIL4", "SIX1", "HES5", "MAGEB2",  "TCF3", "UTF1", "EGR4")
gene.of.interest.spg2 <- c("DPPA2", "SIX1", "FGFR3", "PIWIL4", "EGR4", "ID4", "SOHLH1", "UTF1")
gene.of.interest.earlydiff1 <- c("ASB9", "STRA8","PRSS21", "SYCP3", "MLEC", "PARP4", "REC8", "LARP1", "DNMT1", "PIWIL1", "TOP2B")
gene.of.interest.earlydiff2 <- c("KIT", "PIWIL1", "ASB9", "ETV7", "L1TD1", "SCML1", "MKI67", "TEX15", "HORMAD1", "MEIOC", "SYCP3")
gene.of.interest.latediff1 <- c("SYCP3", "SMC1B", "MEIOC", "MAGEA4", "UBE2C", "ESX1", "MKI67")
gene.of.interest.latediff2 <- c("CENPF", "KIT", "ITGA6", "UBE2C", "ELAVL2",  "RHOXF1")
gene.of.interest.prelep1 <- c("TEX101", "MEIOB", "HORMAD1", "SYCP2", "SYCP1", "ZCWPW1","SCML1", "TOP2A","SYCP3", "SPO11", "PRDM9")

gene.list.r1 <- c(gene.of.interest.PGCE, gene.of.interest.PGCL1)
gene.list.r2 <- c(gene.of.interest.PGCL2, gene.of.interest.M1)
gene.list.r3 <- c(gene.of.interest.M2, gene.of.interest.transitional1)
gene.list.r4 <- c(gene.of.interest.transitional2, gene.of.interest.T1)
gene.list.r5 <- c(gene.of.interest.T1_2, gene.of.interest.spg1)
gene.list.r6 <- c(gene.of.interest.spg2, gene.of.interest.earlydiff1)
gene.list.r7 <- c(gene.of.interest.earlydiff2, gene.of.interest.latediff1)
gene.list.r8 <- c(gene.of.interest.latediff2, gene.of.interest.prelep1)




list.of.lists.of.genes <- list(gene.list.r1,
                               gene.list.r2,
                               gene.list.r3,
                               gene.list.r4,
                               gene.list.r5,
                               gene.list.r6,
                               gene.list.r7,
                               gene.list.r8)



# log2fc.cutoff <- 0.25

for(number in 1:(length(list.of.cell.types)-1)){
  # number <- 1
  print(number)
  # gene.of.interest <- c(list.of.lists.of.genes[[number]], list.of.lists.of.genes[[number+1]])
  # gene.of.interest <- list.of.lists.of.genes[[number]]
  gene.of.interest <- c("NANOS3", "DND1")
  current.comparison <- list.of.comparisons[[number]]
  current.counts <- list.of.counts[[1]] #NOTE THIS IS ONLY COMPARING PGCLCs to everything
  current.counts2 <- list.of.counts[[number+1]]
  
  current.comparison$gene <- rownames(current.comparison)
  # current.comparison.filtered <- filter(current.comparison, p_val_adj < p.val.cutoff)
  current.comparison.filtered.up <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC > 0)
  current.comparison.filtered.down <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC < 0)
  current.comparison.filtered.significant <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff)
  
  
  
  
  current.comparison.dataframe <- data.frame(current.counts, current.counts2)
  current.comparison.dataframe$log2fc <- current.counts - current.counts2
  current.comparison.dataframe$regulation <- "not"
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc > log2fc.cutoff] <- "up"
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc < -log2fc.cutoff] <- "down"
  
  current.comparison.dataframe$regulation[!(current.comparison.dataframe$gene.name %in% current.comparison.filtered.significant$gene)] <- "not"
  
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  # current.comparison.subset <- filter(current.comparison.dataframe, av.counts.diffL > 0 | av.counts.Spermatocytes>0)
  current.comparison.subset <- dplyr::filter(current.comparison.dataframe, current.counts > 0 | current.counts2 >0)
  # current.comparison.sig.up <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.up$gene)
  # current.comparison.sig.down <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.down$gene)
  current.comparison.subset <- na.omit(current.comparison.subset)
  current.comparison.subset$DEGs <- "Not significant"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene] <- "Up"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene] <- "Down"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% gene.of.interest] <- "GeneOfInterest"
  current.comparison.subset$regulation[current.comparison.subset$gene.name %in% gene.of.interest] <- "GeneOfInterest"
  print(paste0(list.of.cell.types[number], " vs ", list.of.cell.types[number+1], 
               " Up: ", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]), 
               " Down:", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])))
  
  colnames(current.comparison.subset) <- c("HS.counts", "NCG.counts", "log2FC_manual", "regulation", "gene.name", "DEGs")
  print(table(current.comparison.subset$regulation))
  
  degs_order <- c( "Down", "GeneOfInterest", "Not significant", "Up")
  
  
  current.comparison.subset_ordered <- current.comparison.subset %>%
    mutate(DEGs = factor(DEGs, levels = degs_order)) 
  
  
  assign(paste0("r", number), 
         
         
         
         ggplot(data = current.comparison.subset_ordered, aes(x = HS.counts, y = NCG.counts,
                                                              label = gene.name,
                                                              color = regulation)) +
           geom_point(data = filter(current.comparison.subset_ordered, regulation == "not"), alpha = 0.5) + # Plot "Not Significant" points first
           geom_point(data = filter(current.comparison.subset_ordered, regulation != "not")) + # Plot other points
           # geom_text(data = filter(current.comparison.subset_ordered, regulation != "not")) +
           geom_point(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest")) +
           geom_text(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest"), hjust = -0.1, vjust = -0.1) +
           
           
           scale_color_manual(values = c(cell.type.colors[number+1], "black", '#999999', cell.type.colors[number])) +
           theme_classic() +
           xlim(0, 3) + ylim(0, 3) + 
           labs(title = paste(list.of.cell.types[1], " vs. ", list.of.cell.types[number + 1], ", q", p.val.cutoff, ", log2fc", log2fc.cutoff),
                x = list.of.cell.types[number], y = list.of.cell.types[number + 1]) +
           theme(plot.title = element_text(size = 9),
                 panel.grid.major = element_line(color = "gray", linetype = "dotted"),
                 panel.grid.minor = element_line(color = "lightgray", linetype = "dotted")) +
           geom_abline(intercept = 0, slope = 1, linetype = "dotted")
  )
}



r1
r2
r3


###############################################################################################################################
# Scatterplot of WT vs KO



table(Idents(big_DND1))
table(Idents(big_NANOS3))

data.dnd1KO <- as.matrix(GetAssayData(big_DND1, slot = "counts")[, WhichCells(big_DND1, ident = "DND1_KO_d14")])
data.dnd1WT <- as.matrix(GetAssayData(big_DND1, slot = "counts")[, WhichCells(big_DND1, ident = "DND1_WT_d14")])
data.nanos3KO <- as.matrix(GetAssayData(big_NANOS3, slot = "counts")[, WhichCells(big_NANOS3, ident = "NANOS3_KO_d9_xrT")])
data.nanos3WT <- as.matrix(GetAssayData(big_NANOS3, slot = "counts")[, WhichCells(big_NANOS3, ident = "NANOS3_WT_d9_xrT")])



av.counts.dnd1KO <- apply(data.dnd1KO, 1, mean)
av.counts.dnd1WT <- apply(data.dnd1WT, 1, mean)
av.counts.nanos3KO <- apply(data.nanos3KO, 1, mean)
av.counts.nanos3WT <- apply(data.nanos3WT, 1, mean)

sum.dnd1KO <- sum(av.counts.dnd1KO)
sum.dnd1WT <- sum(av.counts.dnd1WT)
sum.nanos3KO <- sum(av.counts.nanos3KO)
sum.nanos3WT <- sum(av.counts.nanos3WT)


norm.av.counts.dnd1KO <- log2(((av.counts.dnd1KO) / sum.dnd1KO)*10^6 +1)
norm.av.counts.dnd1WT <- log2(((av.counts.dnd1WT) / sum.dnd1WT)*10^6 +1)
norm.av.counts.nanos3KO <- log2(((av.counts.nanos3KO) / sum.nanos3KO)*10^6 +1)
norm.av.counts.nanos3WT <- log2(((av.counts.nanos3WT) / sum.nanos3WT)*10^6 +1)



###

p.val.cutoff <- 0.05
log2fc.cutoff <- 0.4

###


DND1.KO.v.WT <- FindMarkers(big_DND1, ident.1 = "DND1_KO_d14", ident.2 = "DND1_WT_d14", logfc.threshold = log2fc.cutoff)#, min.pct = 0.2)
NANOS3.KO.v.WT <- FindMarkers(big_NANOS3, ident.1 = "NANOS3_KO_d9_xrT", ident.2 = "NANOS3_WT_d9_xrT", logfc.threshold = log2fc.cutoff)#, min.pct = 0.2)


list.of.comparisons <- list(DND1.KO.v.WT, NANOS3.KO.v.WT)

# list.of.counts <- list(av.counts.dnd1KO, av.counts.dnd1WT, av.counts.nanos3KO, av.counts.nanos3WT)
list.of.counts <- list(norm.av.counts.dnd1KO, norm.av.counts.dnd1WT, norm.av.counts.nanos3KO, norm.av.counts.nanos3WT)


list.of.cell.types <- c("dnd1KO", 
                        "dnd1WT", "nanos3KO", "nanos3WT")



for(number in 1:2){
  
  # number <- 1
  
  current.comparison <- list.of.comparisons[[number]]
  
  if(number == 1){
    current.counts <- list.of.counts[[1]]
    current.counts2 <- list.of.counts[[2]]
  }
  if(number == 2){
    current.counts <- list.of.counts[[3]]
    current.counts2 <- list.of.counts[[4]]
  }
  current.comparison$gene <- rownames(current.comparison)
  # current.comparison.filtered <- filter(current.comparison, p_val_adj < p.val.cutoff)
  
  
  current.comparison.filtered.up <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC > log2fc.cutoff)
  current.comparison.filtered.down <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC < -log2fc.cutoff)
  current.comparison.filtered.significant <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff)
  
  current.comparison.dataframe <- data.frame(current.counts, current.counts2)
  # current.counts.df <- data.frame(current.counts, current.counts2)
  current.comparison.dataframe$log2fc <- current.counts - current.counts2
  current.comparison.dataframe$regulation <- "not"
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc > log2fc.cutoff] <- "up"
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc < -log2fc.cutoff] <- "down"
  current.comparison.dataframe$regulation[!(current.comparison.dataframe$gene.name %in% current.comparison.filtered.significant$gene)] <- "not"
  
  
  # #simply deduct count values
  # current.comparison.dataframe$difference <- current.comparison.dataframe$current.counts - current.comparison.dataframe$current.counts2
  # 
  
  
  
  current.comparison.subset <- dplyr::filter(current.comparison.dataframe, current.counts > 0 | current.counts2 >0)
  # current.comparison.sig.up <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.up$gene)
  # current.comparison.sig.down <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.down$gene)
  current.comparison.subset <- na.omit(current.comparison.subset)
  current.comparison.subset$DEGs <- "Not significant"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene] <- "Up"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene] <- "Down"
  print(paste0(list.of.cell.types[number+(number-1)], " vs ", list.of.cell.types[number+number], 
               " Up: ", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]), 
               " Down:", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])))
  
  print(table(current.comparison.subset$regulation))
  
  
  # assign(paste0("p", number),
  #        (ggplot() + geom_point(data = current.comparison.subset, aes(x = current.counts, y = current.counts2,  color = DEGs)) +
  #           scale_color_manual(values=c('#E69F00','#999999', '#56B4E9')) +     #      xlim(0, 3) + ylim(0, 3) + 
  #           theme_classic() + labs(title=paste(list.of.cell.types[number+(number-1)], " vs. ", list.of.cell.types[number+number], ", q", p.val.cutoff, ", log2fc", log2fc.cutoff),
  #                                  x=list.of.cell.types[number+(number-1)], y = list.of.cell.types[number+number]) + theme(plot.title = element_text(size=9)))
  # )
  
  assign(paste0("p", number),
         (ggplot() + geom_point(data = current.comparison.subset, aes(x = current.counts, y = current.counts2,  color = regulation)) +
            scale_color_manual(values=c('#E69F00','#999999', '#56B4E9')) +     #      xlim(0, 3) + ylim(0, 3) + 
            theme_classic() + labs(title=paste(list.of.cell.types[number+(number-1)], " vs. ", list.of.cell.types[number+number], ", q", p.val.cutoff, ", log2fc", log2fc.cutoff),
                                   x=list.of.cell.types[number+(number-1)], y = list.of.cell.types[number+number]) + theme(plot.title = element_text(size=9)))
  )
  
}

p1
p2
# p3


cell.type.colors <- c("#2179b4", "#b4d88b", "#36a047", "#f6999a", "#f15a29", "#fdbf6f", "#f7941d", "#cab3d6", "#6b3f98")




gene.list.r1 <- c("HMGB2","TCL1A", "TCL1B", "NANOS3", "ACTG2", "SYNE2",
                  "GAL", "NES", "EN2", "SOX11",  "RTN4", "RCN2", "FEZ1", "NECAB1", "NPTN", "SELENOW")
gene.list.r2 <- c("METTL9", "SYNE2", "TCL1B", "PDLIM1", 
                  "ASRGL1", "TUBB2B", "SELENOW", "PSAP", "TBCB", "VAMP8", "TUBGCP2")





list.of.lists.of.genes <- list(gene.list.r1,
                               gene.list.r2)



# log2fc.cutoff <- 0.25

for(number in 1:2){
  # number <- 1
  print(number)
  gene.of.interest <- list.of.lists.of.genes[[number]]
  # gene.of.interest <- list.of.lists.of.genes[[number]]
  
  current.comparison <- list.of.comparisons[[number]]
  current.counts <- list.of.counts[[2*number-1]] #NOTE THIS IS ONLY COMPARING PGCLCs to everything
  current.counts2 <- list.of.counts[[2*number]]
  
  current.comparison$gene <- rownames(current.comparison)
  # current.comparison.filtered <- filter(current.comparison, p_val_adj < p.val.cutoff)
  current.comparison.filtered.up <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC > 0)
  current.comparison.filtered.down <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC < 0)
  current.comparison.filtered.significant <- dplyr::filter(current.comparison, p_val_adj < p.val.cutoff)
  
  
  
  
  current.comparison.dataframe <- data.frame(current.counts, current.counts2)
  current.comparison.dataframe$log2fc_manual <- current.counts - current.counts2
  
  current.comparison.dataframe$regulation <- "not"
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc > log2fc.cutoff] <- "up"
  current.comparison.dataframe$regulation[current.comparison.dataframe$log2fc < -log2fc.cutoff] <- "down"
  
  current.comparison.dataframe$regulation[!(current.comparison.dataframe$gene.name %in% current.comparison.filtered.significant$gene)] <- "not"
  
  current.comparison$regulation <- "NotSignificant"
  current.comparison[current.comparison.filtered.up$gene,]$regulation <- "HigherInKO"
  current.comparison[current.comparison.filtered.down$gene,]$regulation <- "HigherInWT"
  
  
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  # current.comparison.subset <- filter(current.comparison.dataframe, av.counts.diffL > 0 | av.counts.Spermatocytes>0)
  current.comparison.subset <- dplyr::filter(current.comparison.dataframe, current.counts > 0 | current.counts2 >0)
  # current.comparison.sig.up <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.up$gene)
  # current.comparison.sig.down <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.down$gene)
  current.comparison.subset <- na.omit(current.comparison.subset)
  current.comparison.subset$DEGs <- "Not significant"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene] <- "Up"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene] <- "Down"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% gene.of.interest] <- "GeneOfInterest"
  current.comparison.subset$regulation[current.comparison.subset$gene.name %in% gene.of.interest] <- "GeneOfInterest"
  print(paste0(list.of.cell.types[number], " vs ", list.of.cell.types[number+1], 
               " Up: ", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]), 
               " Down:", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])))
  
  
  
  colnames(current.comparison.subset) <- c("HS.counts", "NCG.counts", "log2FC", "regulation", "gene.name", "DEGs")
  print(table(current.comparison.subset$regulation))
  
  degs_order <- c( "Down", "GeneOfInterest", "Not significant", "Up")
  
  
  
  current.comparison.subset_ordered <- current.comparison.subset %>%
    mutate(DEGs = factor(DEGs, levels = degs_order)) 
  setwd("~/Documents/PROJECTS/ksasaki/NANOS3KO")
  if(number == 1){write.csv(current.comparison, file = "DND1-KO_DEGs.csv")}
  if(number == 2){write.csv(current.comparison, file = "NANOS3-KO_DEGs.csv")}
  
  assign(paste0("r", number), 
         
         
         
         ggplot(data = current.comparison.subset_ordered, aes(x = HS.counts, y = NCG.counts,
                                                              label = gene.name,
                                                              color = DEGs)) +
           geom_point(data = filter(current.comparison.subset_ordered, DEGs == "Not significant"), alpha = 0.5) + # Plot "Not Significant" points first
           geom_point(data = filter(current.comparison.subset_ordered, regulation != "Not significant")) + # Plot other points
           # geom_text(data = filter(current.comparison.subset_ordered, regulation != "Not significant")) +
           geom_point(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest")) +
           geom_text(data = filter(current.comparison.subset_ordered, regulation == "GeneOfInterest"), hjust = -0.1, vjust = -0.1) +
           
           
           scale_color_manual(values = c(  '#999999',"darkred","black", "darkgreen")) +
           theme_classic() +
           xlim(0, 3.5) + ylim(0, 3.5) +
           labs(title = paste(list.of.cell.types[number+number-1], " vs. ", list.of.cell.types[2*number], ", q", p.val.cutoff, ", log2fc", log2fc.cutoff),
                x = list.of.cell.types[2*number-1], y = list.of.cell.types[2*number]) +
           theme(plot.title = element_text(size = 9),
                 panel.grid.major = element_line(color = "gray", linetype = "dotted"),
                 panel.grid.minor = element_line(color = "lightgray", linetype = "dotted")) +
           geom_abline(intercept = 0, slope = 1, linetype = "dotted")
  )
}


# Extended Data Fig 5 r and s

plot_grid(r1, r2)
r1
r2



### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### PART 4              ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### Macaque xrTestis    ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

library(Seurat)
library(SeuratWrappers)
library(SeuratObject)
library(phateR)
library(RColorBrewer)
library(scCustomize)
library(ggridges)
library(dplyr)
library(stringr)
library(patchwork)
library(viridis)
library(DoubletFinder)
library(monocle3)
library(gplots)
library(forcats)
library(ComplexHeatmap)
library(ggplot2)
library(circlize)
library(tidyverse)
library(cowplot)
library(lattice)
library(scales)

setwd("/Users/ewhelan/Documents/PROJECTS/ksasaki/Macaque/")

#Read in macaque 10X data. As with human, data is aligned against macaque reference and mouse separately in order to distinguish macaque from mouse later.

Count_output_1C2AGVTDay170_v_macaque <- Read10X(data.dir = "Count_output_1C2AGVTDay170_v_macaque/filtered_feature_bc_matrix")
Count_output_1C2AGVTDay170_v_mouse <- Read10X(data.dir = "Count_output_1C2AGVTDay170_v_mouse/raw_feature_bc_matrix")

Count_output_1C2AGVT9M_v_macaque <- Read10X(data.dir = "Count_output_1C2AGVT9M_v_macaque/filtered_feature_bc_matrix")
Count_output_1C2AGVT9M_v_mouse <- Read10X(data.dir = "Count_output_1C2AGVT9M_v_mouse/raw_feature_bc_matrix")

Count_output_1C2_4M_transplant_v_macaque <- Read10X(data.dir = "Count_output_1C2_4M_transplant_v_macaque/filtered_feature_bc_matrix")
Count_output_1C2_4M_transplant_v_mouse <- Read10X(data.dir = "Count_output_1C2_4M_transplant_v_mouse/raw_feature_bc_matrix")

Count_output_2M_Transplant_c38_v_macaque <- Read10X(data.dir = "Count_output_2M_Transplant_c38_v_macaque/filtered_feature_bc_matrix")
Count_output_2M_Transplant_c38_v_mouse <- Read10X(data.dir = "Count_output_2M_Transplant_c38_v_mouse/raw_feature_bc_matrix")

Count_output_1M_Transplant_c23_v_macaque <- Read10X(data.dir = "Count_output_1M_Transplant_c23_v_macaque/filtered_feature_bc_matrix")
Count_output_1M_Transplant_c23_v_mouse <- Read10X(data.dir = "Count_output_1M_Transplant_c23_v_mouse/raw_feature_bc_matrix")

Count_output_01_Yasu_10X_1c2_c34_5M_trans_v_macaque <- Read10X(data.dir = "Count_output_01_Yasu_10X_1c2_c34_5M_trans_v_macaque/filtered_feature_bc_matrix")
Count_output_01_Yasu_10X_1c2_c34_5M_trans_v_mouse <- Read10X(data.dir = "Count_output_01_Yasu_10X_1c2_c34_5M_trans_v_mouse/raw_feature_bc_matrix")

macaque.list <- c(Count_output_1C2AGVTDay170_v_macaque, Count_output_1C2AGVT9M_v_macaque, Count_output_1C2_4M_transplant_v_macaque,
                  Count_output_2M_Transplant_c38_v_macaque, Count_output_1M_Transplant_c23_v_macaque,
                  Count_output_01_Yasu_10X_1c2_c34_5M_trans_v_macaque)
mouse.list <- c(Count_output_1C2AGVTDay170_v_mouse, Count_output_1C2AGVT9M_v_mouse, Count_output_1C2_4M_transplant_v_mouse, 
                Count_output_2M_Transplant_c38_v_mouse, Count_output_1M_Transplant_c23_v_mouse, 
                Count_output_01_Yasu_10X_1c2_c34_5M_trans_v_mouse)

#names refere to time transplanted in months
names.data <- c("macaque6M", "macaque9M", "macaque4M", 
                "macaque2M", "macaque1M", 
                "macaque5M")

macaque.list.original <- macaque.list

mouse.or.macaque.list <- macaque.list


#this assignes macaque or mouse identity on a per-cell basis.

for(current.dataset in 1:length(macaque.list)){
  # current.dataset <- 6
  current.dataset.hu <- macaque.list[[current.dataset]]
  current.dataset.mu <- mouse.list[[current.dataset]]
  colnames.hu <- colnames(current.dataset.hu)
  length(colnames.hu)
  current.dataset.mu.filtered <- current.dataset.mu[,colnames.hu]
  colnames.mu <- colnames(current.dataset.mu.filtered)
  length(colnames.mu)
  number.of.cells <- length(intersect(colnames.hu, colnames.mu))
  if(number.of.cells != length(colnames.hu)){print("ERROR!!!")}
  counts.hu <- colSums(current.dataset.hu)
  counts.mu <- colSums(current.dataset.mu.filtered)
  macaque.or.mouse <- counts.hu
  for(current.cell in 1:(length(macaque.or.mouse))){
    # current.cell <- 1
    # if((counts.hu[current.cell])>6000){
    if((counts.hu[current.cell]/counts.mu[current.cell])>2){
      macaque.or.mouse[current.cell] <- "macaque"
    } else {
      if((counts.hu[current.cell]/counts.mu[current.cell])<0.5){
        macaque.or.mouse[current.cell] <- "mouse"
      } else {
        macaque.or.mouse[current.cell] <- "unknown"
      }
    }
    # }else{macaque.or.mouse[current.cell] <- "lowQ"}
  }
  print(names.data[current.dataset])
  print(table(macaque.or.mouse))
  mouse.or.macaque.list[[current.dataset]] <- macaque.or.mouse
  
}

#examine the UMI/cell:

for(a in 1:length(macaque.list)){
  print(length(Cells(macaque.list[[a]])))
}

inflection.point1 <- c(1:length(macaque.list))
log_lib_size_at_inflection1 <- c(1:length(macaque.list))

par(mfrow = c(3,2))

for(y in 1:length(macaque.list)){
  # y=2
  print(names.data[y])
  expression.df <- as.data.frame(macaque.list[y])
  
  umi_per_barcode <- colSums(expression.df)
  barcode_rank <- rank(-umi_per_barcode)
  # plot(barcode_rank, umi_per_barcode,
  #      #ylim=c(1,2000),
  #      xlim=c(1,10000))
  
  log_lib_size <- log10(umi_per_barcode)
  #plot(barcode_rank, log_lib_size, xlim=c(1,10000))
  
  o <- order(barcode_rank)
  log_lib_size <- log_lib_size[o]
  barcode_rank <- barcode_rank[o]
  
  rawdiff <- diff(log_lib_size)/diff(barcode_rank)
  # inflection <- which(rawdiff == min(rawdiff[500:4000], na.rm=TRUE))
  
  
  
  plot(barcode_rank, log_lib_size, xlim=c(1,8000),
       pch = 20,
       main = paste(y, names.data[y])
       # main = paste(y, names.data[y], "cells", inflection, sep="_", "UMIs", round(10^log_lib_size[inflection]))
  )
  
}


#here flat cutoffs were used:

min.features <- 1000
max.features <- 10000
max.mito <- 20
min.UMIs <- 10^4



all.gene.human <- rownames(NCG_germ4)
all.gene.human.df <- as.data.frame(all.gene.human)
mt_genes <- all.gene.human[grep("^MT-", all.gene.human)]
gene_names_cleaned <- sub("^MT-", "", mt_genes)

all.gene.macaque <- rownames(Count_output_1C2AGVTDay170_v_macaque)
all.gene.macaque.df <- as.data.frame(all.gene.macaque)
length(all.gene.macaque)
# Add "MT-" prefix to gene names in the larger list if they match the cleaned gene names
gene_names_with_prefix <- all.gene.macaque

for (i in seq_along(gene_names_with_prefix)) {
  if (gene_names_with_prefix[i] %in% gene_names_cleaned) {
    gene_names_with_prefix[i] <- paste0("MT-", gene_names_with_prefix[i])
  }
}

# Print the result
print(gene_names_with_prefix)
length(gene_names_with_prefix)

rownames(Count_output_1C2AGVTDay170_v_macaque) <- gene_names_with_prefix
rownames(Count_output_1C2AGVT9M_v_macaque) <- gene_names_with_prefix
rownames(Count_output_1C2_4M_transplant_v_macaque) <- gene_names_with_prefix
rownames(Count_output_2M_Transplant_c38_v_macaque) <- gene_names_with_prefix
rownames(Count_output_1M_Transplant_c23_v_macaque) <- gene_names_with_prefix
rownames(Count_output_01_Yasu_10X_1c2_c34_5M_trans_v_macaque) <- gene_names_with_prefix


#make seurats

macaque6m <- CreateSeuratObject(Count_output_1C2AGVTDay170_v_macaque, project = "macaque6m")
macaque9m <- CreateSeuratObject(Count_output_1C2AGVT9M_v_macaque, project = "macaque9m")
macaque4m <- CreateSeuratObject(Count_output_1C2_4M_transplant_v_macaque, project = "macaque4m")
macaque2m <- CreateSeuratObject(Count_output_2M_Transplant_c38_v_macaque, project = "macaque2m")
macaque1m <- CreateSeuratObject(Count_output_1M_Transplant_c23_v_macaque, project = "macaque1m")
macaque5m <- CreateSeuratObject(Count_output_01_Yasu_10X_1c2_c34_5M_trans_v_macaque, project = "macaque5m")


list.of.seurats <- c(macaque6m, macaque9m, 
                     macaque4m, macaque2m, macaque1m, macaque5m)

Cells.merged <- merge(
  x = list.of.seurats[[1]],
  y = list.of.seurats[2:length(list.of.seurats)])

VlnPlot(object = Cells.merged,
        #log = FALSE,
        pt.size = 0,
        group.by = "orig.ident",
        #y.max = 2000,
        #ncol = 3,
        #cols = c("blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue"),
        features = c("nFeature_RNA", "nCount_RNA", "percent.mito"))



for(current.sample in 1:length(list.of.seurats)){
  # current.sample <- 1
  current.seurat <- list.of.seurats[[current.sample]]
  current.mouse.or.macaque <- mouse.or.macaque.list[[current.sample]]
  current.seurat$species <- current.mouse.or.macaque
  
  current.seurat[['percent.mito']] <- PercentageFeatureSet(current.seurat, pattern = "^MT-")
  list.of.seurats[[current.sample]] <- current.seurat
}

for(current.sample in 1:(length(list.of.seurats))){
  # current.sample <- 2
  
  current.seurat <- list.of.seurats[[current.sample]]
  current.seurat <- subset(x = current.seurat,
                           subset = nFeature_RNA > min.features & 
                             nFeature_RNA < max.features & 
                             percent.mito < max.mito & nCount_RNA > min.UMIs)
  # subset = percent.mito < max.mito & nCount_RNA > min.UMIs) #
  # current.seurat <- subset(x = current.seurat, subset = species == "macaque")
  list.of.seurats[[current.sample]] <- current.seurat
}

#statistics:

number.of.samples <- length(list.of.seurats)
number.of.cells.per.sample <-  data.frame(
  SampleName = c(1:number.of.samples),
  CellNumber = c(1:number.of.samples),
  GeneNumber = c(1:number.of.samples),
  SampleNumber = c(1:number.of.samples),
  SampleOrigin = c(1:number.of.samples),
  medianUMI = c(1:number.of.samples),
  medianGenes = c(1:number.of.samples),
  medianMito =c(1:number.of.samples),
  cutoff = c(1:number.of.samples)
)

for(a in 1:length(list.of.seurats)){
  number.of.cells.per.sample[a,4]<-a
  number.of.cells.per.sample[a,1]<-names.data[a]
  number.of.cells.per.sample[a,2]<-length(Cells(list.of.seurats[[a]]))
  number.of.cells.per.sample[a,3]<-length(rownames(list.of.seurats[[a]]))
  number.of.cells.per.sample[a,5]<-"na"
  number.of.cells.per.sample[a,6]<-median(list.of.seurats[[a]]$nCount_RNA)
  number.of.cells.per.sample[a,7]<-median(list.of.seurats[[a]]$nFeature_RNA)
  number.of.cells.per.sample[a,8]<-median(list.of.seurats[[a]]$percent.mito)
  # number.of.cells.per.sample[a,9]<-min.UMIs1[a]
  number.of.cells.per.sample[a,9]<-min.UMIs
  #print(number.of.cells.per.sample[a,])
}

print(number.of.cells.per.sample)

Cells.merged <- merge(
  x = list.of.seurats[[1]],
  y = list.of.seurats[2:length(list.of.seurats)])

VlnPlot(object = Cells.merged,
        #log = FALSE,
        pt.size = 0,
        group.by = "orig.ident",
        #y.max = 2000,
        #ncol = 3,
        #cols = c("blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue", "blue"),
        features = c("nFeature_RNA", "nCount_RNA", "percent.mito")
)

# #examine seurats one at a time.
# 
# list.of.seurats1 <- list.of.seurats
# list.of.seurats2 <- list.of.seurats
# 
# 
# for(a in 1:length(list.of.seurats1)){
#   # a <- 1
#   current.seurat <- list.of.seurats1[[a]]
#   current.seurat <- NormalizeData(object = current.seurat, verbose = FALSE)
#   current.seurat  <- FindVariableFeatures(object = current.seurat, nfeatures = 2000, verbose = FALSE)
#   current.seurat <- ScaleData(object = current.seurat, vars.to.regress = c("nCount_RNA", "percent.mito"), verbose = FALSE)
#   current.seurat <- RunPCA(object = current.seurat,
#                              npcs = 30,
#                              verbose = TRUE)
#   current.seurat <- FindNeighbors(object = current.seurat)
#   current.seurat <- FindClusters(object = current.seurat,
#                                    resolution = 1.1,
#                                    #algorithm = 1,
#                                    verbose = TRUE
#   )
#   # current.seurat <- RunTSNE(object = current.seurat, dims = 1:20)
#   current.seurat <- RunUMAP(object = current.seurat, reduction = "pca",
#                                 dims = 1:30)
# 
#   p1 <- (DimPlot(object = current.seurat, reduction = "umap", group.by = "species",label = F) & ggtitle(names.data[a]))
#   p2 <- (FeaturePlot(object = current.seurat, reduction = "umap", features = "REC8",label = T) & ggtitle(names.data[a]) & DarkTheme() & scale_color_viridis())
#   
#   print(plot_grid(p1,p2))
#   
#   current.seurat -> list.of.seurats1[[a]]
#   
#   Idents(current.seurat) <- "species"
#   current.seurat <- subset(current.seurat, idents = "macaque")
#   
#   current.seurat <- NormalizeData(object = current.seurat, verbose = FALSE)
#   current.seurat  <- FindVariableFeatures(object = current.seurat, nfeatures = 2000, verbose = FALSE)
#   current.seurat <- ScaleData(object = current.seurat, vars.to.regress = c("nCount_RNA", "percent.mito"), verbose = FALSE)
#   current.seurat <- RunPCA(object = current.seurat,
#                            npcs = 30,
#                            verbose = TRUE)
#   current.seurat <- FindNeighbors(object = current.seurat)
#   current.seurat <- FindClusters(object = current.seurat,
#                                  resolution = 1.1,
#                                  #algorithm = 1,
#                                  verbose = TRUE
#   )
#   # current.seurat <- RunTSNE(object = current.seurat, dims = 1:20)
#   current.seurat <- RunUMAP(object = current.seurat, reduction = "pca",
#                             dims = 1:30)
#   
#   p1 <- (DimPlot(object = current.seurat, reduction = "umap", group.by = "species",label = F) & ggtitle(names.data[a]))
#   p2 <- (FeaturePlot(object = current.seurat, reduction = "umap", features = "REC8",label = T) & ggtitle(names.data[a]) & DarkTheme() & scale_color_viridis())
#   
#   print(plot_grid(p1,p2))
#   
#   
#   
#   current.seurat -> list.of.seurats2[[a]]
#   
# }
# 
# 
# gene.of.interest <- "SYCP3"
# p1 <- FeaturePlot(list.of.seurats2[[1]], features = gene.of.interest, order =T) & DarkTheme() & scale_color_viridis() & ggtitle(paste(gene.of.interest, names.data[1]))
# p2 <- FeaturePlot(list.of.seurats2[[2]], features = gene.of.interest, order =T) & DarkTheme() & scale_color_viridis()& ggtitle(paste(gene.of.interest, names.data[2]))
# p3 <- FeaturePlot(list.of.seurats2[[3]], features = gene.of.interest, order =T) & DarkTheme() & scale_color_viridis()& ggtitle(paste(gene.of.interest, names.data[3]))
# p4 <- FeaturePlot(list.of.seurats2[[4]], features = gene.of.interest, order =T) & DarkTheme() & scale_color_viridis()& ggtitle(paste(gene.of.interest, names.data[4]))
# p5 <- FeaturePlot(list.of.seurats2[[5]], features = gene.of.interest, order =T) & DarkTheme() & scale_color_viridis()& ggtitle(paste(gene.of.interest, names.data[5]))
# p6 <- FeaturePlot(list.of.seurats2[[6]], features = gene.of.interest, order =T) & DarkTheme() & scale_color_viridis()& ggtitle(paste(gene.of.interest, names.data[6]))
# plot_grid(p5, p4, p3, p6, p1, p2)


# current.list <- c(HS_SCT_subset, list.of.seurats)
current.list <- list.of.seurats


#integrate samples together:

for (i in 1:length(x = current.list)) {
  current.list[[i]] <- NormalizeData(object = current.list[[i]], verbose = FALSE)
  current.list[[i]] <- FindVariableFeatures(object = current.list[[i]], nfeatures = 3000, verbose = FALSE)
}

current.anchors <- FindIntegrationAnchors(object.list = current.list, dims = 1:30)
# save(current.anchors, file="current.anchors.Robj")
current.integrated <- IntegrateData(anchorset = current.anchors, dims = 1:30)
# save(current.integrated, file="current.integrated.Robj")

DefaultAssay(object = current.integrated) <- "integrated"
current.integrated <- ScaleData(object = current.integrated, vars.to.regress = c("nCount_RNA", "percent.mito"), verbose = FALSE)
current.integrated <- RunPCA(object = current.integrated,
                             npcs = 30,
                             verbose = TRUE)
current.integrated <- FindNeighbors(object = current.integrated)
current.integrated <- FindClusters(object = current.integrated,
                                   resolution = 1.1,
                                   #algorithm = 1,
                                   verbose = TRUE
)
# current.integrated <- RunTSNE(object = current.integrated, dims = 1:20)
current.integrated <- RunUMAP(object = current.integrated, reduction = "pca",
                              dims = 1:30)

DimPlot(object = current.integrated, reduction = "umap", group.by = "orig.ident", label = T)


NCG_germ_macaque <- current.integrated
p1 <- DimPlot(object = NCG_germ_macaque, reduction = "umap", group.by = "orig.ident", label = T)
p2 <- DimPlot(object = NCG_germ_macaque, reduction = "umap", group.by = "species", label = T)
print(p1 | p2)

#CELL CYCLE REGRESSION




#These are Seurat's default list of cell cycle genes

s.genes <- (cc.genes.updated.2019$s.genes)
g2m.genes <- (cc.genes.updated.2019$g2m.genes)

#Cell cycle scoring:
DefaultAssay(NCG_germ_macaque) <- "RNA"
NCG_germ_macaque <- CellCycleScoring(NCG_germ_macaque, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


fp1 <- FeaturePlot(NCG_germ_macaque, features = c("S.Score"), min.cutoff = 0) + labs(title = "S score")
fp2 <- FeaturePlot(NCG_germ_macaque, features = c("G2M.Score"), min.cutoff = 0) + labs(title = expression(bold(paste(G[2],"-M score",sep=""))))
wrap_plots(fp1, fp2)

DimPlot(NCG_germ_macaque, group.by = "Phase")

setwd("/Users/ewhelan/Documents/scRNAseq/Rscripts/savefiles")
save(NCG_germ_macaque, file = "NCG_germ_macaque.2023.12.21.Robj")

FeaturePlot(NCG_germ_macaque, features = "REC8", order = T, split.by = "orig.ident") & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(NCG_germ_macaque, features = "nCount_RNA", max.cutoff = 20000) & DarkTheme() & scale_color_viridis()


#Seurat's typical cell cycle regression did not do a good job here, so not used but code retained for comparison:
#Regress out the cell cycle genes
# 
# DefaultAssay(NCG_germ_macaque) <- "integrated"
# NCG_germ_macaque2 <- NCG_germ_macaque
# 
# NCG_germ_macaque2 <- ScaleData(NCG_germ_macaque2, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(NCG_germ_macaque))
# 
# 
# 
# #Germ.culture2 is the variation with cell cycle regression
# DefaultAssay(NCG_germ_macaque2) <- "integrated"
# # NCG_germ_macaque2 <- FindVariableFeatures(object = NCG_germ_macaque2, 
# #                                       verbose = FALSE)
# n.pcs <- 30
# NCG_germ_macaque2 <- RunPCA(NCG_germ_macaque2, npcs = n.pcs)#, features = VariableFeatures(Germ.culture),npcs = 30, nfeatures.print = 10)
# NCG_germ_macaque2 <- FindNeighbors(object = NCG_germ_macaque2,
#                            k.param = n.pcs)
# NCG_germ_macaque2 <- FindClusters(object = NCG_germ_macaque2, 
#                           # resolution = 0.34,
#                           algorithm = 1,
#                           verbose = TRUE
# )
# # Germ.culture2<- RunTSNE(object = Germ.culture2, dims = 1:n.pcs, reduction = "pca")
# NCG_germ_macaque2 <- RunUMAP(object = NCG_germ_macaque2, 
#                      dims = 1:n.pcs,
#                      reduction = "pca")
# 
# 
# DefaultAssay(NCG_germ_macaque2) <- "RNA"
# DimPlot(NCG_germ_macaque2, #split.by = "experiment", 
#         group.by = "Phase")
# 


# Instead cell cycle was regressed by splitting and re-integrating samples:
DefaultAssay(NCG_germ_macaque) <- "RNA"
Idents(NCG_germ_macaque) <- "species"
DimPlot(NCG_germ_macaque)
NCG_germ_macaque_nomouse <- subset(NCG_germ_macaque, idents = "macaque")


FeaturePlot(NCG_germ_macaque_nomouse, features = c("POU5F1", "DDX4", "PIWIL4", "STRA8", "REC8", "MEIOB"))


Idents(NCG_germ_macaque_nomouse) <- "Phase"

NCG_S <- subset(NCG_germ_macaque_nomouse, idents = "S")
NCG_G1 <- subset(NCG_germ_macaque_nomouse, idents = "G1")
NCG_G2M <- subset(NCG_germ_macaque_nomouse, idents = "G2M")

n.dims <- 25

current.list <- c(NCG_S, NCG_G1, NCG_G2M)
for (i in 1:length(x = current.list)) {
  # current.list[[i]] <- NormalizeData(object = current.list[[i]], verbose = FALSE)
  current.list[[i]] <- FindVariableFeatures(object = current.list[[i]], nfeatures = 3000, verbose = FALSE)
}

current.anchors <- FindIntegrationAnchors(object.list = current.list, dims = 1:n.dims)
# save(current.anchors, file="current.anchors.Robj")
current.integrated <- IntegrateData(anchorset = current.anchors, dims = 1:n.dims)
# save(current.integrated, file="current.integrated.Robj")


DefaultAssay(object = current.integrated) <- "integrated"
current.integrated <- ScaleData(object = current.integrated, verbose = FALSE)
current.integrated <- RunPCA(object = current.integrated,
                             npcs = n.dims,
                             verbose = TRUE)
current.integrated <- FindNeighbors(object = current.integrated)
current.integrated <- FindClusters(object = current.integrated,
                                   resolution = 1.1,
                                   #algorithm = 1,
                                   verbose = TRUE
)
# current.integrated <- RunTSNE(object = current.integrated, dims = 1:20)
current.integrated <- RunUMAP(object = current.integrated, reduction = "pca",
                              dims = 1:n.dims)

p4 <- DimPlot(object = current.integrated, reduction = "umap", group.by = "seurat_clusters", label = T)
p1 <- DimPlot(object = current.integrated, reduction = "umap", group.by = "Phase", label = T)
p2 <- DimPlot(object = current.integrated, reduction = "umap", group.by = "orig.ident")
p3 <- DimPlot(object = current.integrated, reduction = "umap", group.by = "species")
wrap_plots(p4, p1, p2, p3)

NCG_germ_macaque3 <- current.integrated

setwd("/Users/ewhelan/Documents/scRNAseq/Rscripts/savefiles")
save(NCG_germ_macaque3, file="NCG_germ_macaque3_2023.12.18")
# load("NCG_germ_macaque3_cells_renamed_2023.07.26")

DimPlot(NCG_germ_macaque3)


fp1 <- FeaturePlot(NCG_germ_macaque3, features = c("S.Score"), min.cutoff = 0) + labs(title = "S score")
fp2 <- FeaturePlot(NCG_germ_macaque3, features = c("G2M.Score"), min.cutoff = 0) + labs(title = expression(bold(paste(G[2],"-M score",sep=""))))
wrap_plots(fp1, fp2)

#pseudotime for NCG samples

# setwd("/Users/ewhelan/Documents/scRNAseq/Rscripts/savefiles")
# load("NCG_germ_macaque3_redone_2023.06.18")
# DimPlot(NCG_germ_macaque3, reduction = "PHATE", group.by = "seurat_clusters", label = T)
# DimPlot(NCG_germ_macaque3, reduction = "umap", group.by = "seurat_clusters", label = T)


Keren.gene.list <- c("POU5F1", "NANOG", "TCL1A", "PDPN", "PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1",
                     "UTF1", "KIT", "ASB9", "NANOS2", "PIWIL4", "MORC1", "MAGEB2", "TDRD1", "RHOXF1", "MAGEC2", "ZBTB16", "EGR4", "RET", "MAGEA4",
                     "SOHLH1", "DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8", "DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9", "ZCWPW1", "RAD51AP2", "SYCP3",
                     "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1", "MLH3", "PIWIL1", "MLH1", "SPATA8", "CCDC112", "OVOL1", "OVOL2",
                     "CCDC42", "ACR", "CDRT15", "TEKT1", "H1FOO", "TJP3", "FAM24A", "CATSPER3", "SPACA3", "SPACA1", "ACRV1", "SAXO1", "CAPZA3", "HEMGN", "CA2", "PRM3",
                     "TNP2", "TSSK6", "SMCP", "SPEM2", "LELP1", "MKI67", "UBE2C", "TOP2A")

DefaultAssay(NCG_germ_macaque3) <- "RNA"
FeaturePlot(NCG_germ_macaque3, order = T, features = c("POU5F1", "PRDM1", "KIT", "PIWIL4", "DAZL", "ZBTB16", "SOHLH1",
                                                       "STRA8", "REC8", "TEX101", "MEIOB", "MEIOC")) & DarkTheme() & scale_color_viridis(option = "C")

FeaturePlot(NCG_germ_macaque3, order = T, features = c(Keren.gene.list[46:56])) & DarkTheme() & scale_color_viridis(option = "C")



gene.list.PGC <- c("POU5F1", "NANOG", "TCL1A", "PDPN", "PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1",
                   "UTF1", "KIT", "ASB9", "NANOS2")
gene.list.spermatogonia <- c(                 
  "PIWIL4", "MORC1", "MAGEB2", "TDRD1", "RHOXF1", "MAGEC2", "ZBTB16", "EGR4", "RET")
gene.list.diff.spermatgonia <- c(
  
  "MAGEA4",
  "SOHLH1", "DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8")

gene.list.early.meiosis <- c(
  "DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9", "ZCWPW1", "RAD51AP2", "SYCP3",
  "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1", "MLH3")

gene.list.late.meiosis <- c(
  "PIWIL1", "MLH1", "SPATA8", "CCDC112", "OVOL1", "OVOL2",
  "CCDC42", "ACR", "CDRT15", "TEKT1", "H1FOO", "TJP3", "FAM24A", "CATSPER3", "SPACA3", "SPACA1", "ACRV1", "SAXO1")
gene.list.spermatids <- c("CAPZA3", "HEMGN", "CA2", "PRM3",
                          "TNP2", "TSSK6", "SMCP", "SPEM2", "LELP1", "MKI67", "UBE2C", "TOP2A")

FeaturePlot(NCG_germ_macaque3, order = T, features = gene.list.PGC) & DarkTheme() & scale_color_viridis(option = "C")


gene.list.sertoli <- c("CLU", "AMHR2", "SOX9", "CTSL")
genes.macrophage <- c("APOE", "DAB2", "CD74", "ADGRE1")
genes.endothelial <- c("VWF", "TIE1", "TEK")
genes.myoid <- c("ACTA2", "MYH11", "MYL6", "PDGFRB")
genes.leydig <- c("CYP17A1", "CYP11A1", "STAR", "HSD3B1")
genes.T <- c("CD8A", "CD4", "CD3E", "CD28")
genes.telocytes <- c("DCN", "GSN", "TCF21")

germ.cells <- c("POU5F1", "NANOG", "ZBTB16", "SOHLH1", "PRDM9", "ACRV1", "SPACA1", "PRM1", "TNP2")
somatic.cells <- c("SOX9", "APOE", "VWF", "ACTA2", "STAR", "DCN", "TCF21", "ACTA2", "DAB2")

gene.list <- c("POU5F1", "TCL1A", "KIT", "NANOS2", "MAGEB2", "TDRD1", "STRA8", "SOHLH2", "SYCP1" )






# apoptotic.genes <- c("CASP2", "CASP3", "CASP7", "CASP8", "CASP9", "BAX", "BAK1", "BAD")
apoptotic.genes <- c("ADD1","AIFM3","ANKH","ANXA1","APP","ATF3","AVPR1A","BAX","BCAP31","BCL10","BCL2L1",
                     "BCL2L10","BCL2L11","BCL2L2","BGN","BID","BIK","BIRC3","BMF","BMP2","BNIP3L","BRCA1",
                     "BTG2","BTG3","CASP1","CASP2","CASP3","CASP4","CASP6","CASP7","CASP8","CASP9","CAV1",
                     "CCNA1","CCND1","CCND2","CD14","CD2","CD38","CD44","CD69","CDC25B","CDK2","CDKN1A","CDKN1B",
                     "CFLAR","CLU","CREBBP","CTH","CTNNB1","CYLD","DAP","DAP3","DCN","DDIT3","DFFA","DIABLO","DNAJA1","DNAJC3","DNM1L","DPYD","EBP",
                     "EGR3","EMP1","ENO2","ERBB2","ERBB3","EREG","ETF1","F2","F2R","FAS","FASLG","FDXR","FEZ1","GADD45A","GADD45B","GCH1","GNA15",#"GPX1",
                     "GPX3","GPX4","GSN","GSR","GSTM1","GUCY2D",#"H1-0",
                     "HGF","HMGB2","HMOX1","HSPB1","IER3","IFITM3","IFNB1","IFNGR1","IGF2R","IGFBP6",
                     "IL18","IL1A","IL1B","IL6","IRF1","ISG20","JUN","KRT18","LEF1","LGALS3","LMNA","PLPPR4","LUM","MADD","MCL1","MGMT","MMP2","NEDD9",
                     "NEFH","PAK1","PDCD4","PDGFRB","PEA15","PLAT","PLCB2","PMAIP1","PPP2R5B","PPP3R1","PPT1","PRF1","PSEN1","PSEN2","PTK2","RARA",
                     "RELA","RETSAT","RHOB","RHOT2","RNASEL","ROCK1","SAT1","SATB1","SC5D","SLC20A1","SMAD7","SOD1","SOD2","SPTAN1","SQSTM1","TAP1",
                     "TGFB2","TGFBR3","TIMP1","TIMP2","TIMP3","TNF","TNFRSF12A","TNFSF10","TOP2A","TSPO","TXNIP","VDAC2","WEE1","XIAP")


gene.list.PGC.E <- c("POU5F1", "NANOG", "TCL1A", "PDPN" )
gene.list.PGC.L <- c("PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1", "UTF1")
gene.list.prospermatogonia <- c("KIT", 
                                # "ASB9", 
                                "NANOS2")
gene.list.spermatogonia <- c(                 
  "PIWIL4", "MORC1", "MAGEB2", "TDRD1", "RHOXF1", "MAGEC2", "ZBTB16", "EGR4", "RET")

gene.list.diff.spermatgonia.E <- c(
  "MAGEA4",
  "SOHLH1")

gene.list.diff.spermatgonia.L   <- c("DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8")

gene.list.preleptotene <- c("DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9")

gene.list.early.meiosis <- c("ZCWPW1", "RAD51AP2", "SYCP3",
                             "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1", "MLH3")

Idents(NCG_germ_macaque3) <- "Phase"
DimPlot(NCG_germ_macaque3)


DimPlot(NCG_germ_macaque3, group.by = "seurat_clusters",label = T)







NCG_germ_macaque3_temp <- NCG_germ_macaque3
DefaultAssay(NCG_germ_macaque3_temp) <- "RNA"

NCG_germ_macaque3_temp <- AddModuleScore(NCG_germ_macaque3_temp,
                                         features = list(apoptotic.genes),
                                         name="Apoptotic")

NCG_germ_macaque3_temp <- AddModuleScore(NCG_germ_macaque3_temp,
                                         features = list(gene.list.PGC.E),
                                         name="PGCE")


NCG_germ_macaque3_temp <- AddModuleScore(NCG_germ_macaque3_temp,
                                         features = list(gene.list.PGC.L),
                                         name="PGCL")

NCG_germ_macaque3_temp <- AddModuleScore(NCG_germ_macaque3_temp,
                                         features = list(gene.list.prospermatogonia),
                                         name="prospermatogonia")


NCG_germ_macaque3_temp <- AddModuleScore(NCG_germ_macaque3_temp,
                                         features = list(gene.list.spermatogonia),
                                         name="spg")


NCG_germ_macaque3_temp <- AddModuleScore(NCG_germ_macaque3_temp,
                                         features = list(gene.list.diff.spermatgonia.E),
                                         name="diffspgE")

NCG_germ_macaque3_temp <- AddModuleScore(NCG_germ_macaque3_temp,
                                         features = list(gene.list.diff.spermatgonia.L),
                                         name="diffspgL")

NCG_germ_macaque3_temp <- AddModuleScore(NCG_germ_macaque3_temp,
                                         features = list(gene.list.preleptotene),
                                         name="preleptotene")

NCG_germ_macaque3_temp <- AddModuleScore(NCG_germ_macaque3_temp,
                                         features = list(gene.list.early.meiosis),
                                         name="meiosis")


# scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

# Plot scores

plot_grid(
  
  FeaturePlot(NCG_germ_macaque3_temp,
              features = "Apoptotic1", label = F, repel = FALSE) +
    scale_color_viridis(option="plasma") + DarkTheme(),
  
  FeaturePlot(NCG_germ_macaque3_temp,
              features = "PGCE1", label = F, repel = FALSE) +
    scale_color_viridis(option="plasma") + DarkTheme(),
  
  FeaturePlot(NCG_germ_macaque3_temp,
              features = "PGCL1", label = F, repel = FALSE) +
    scale_color_viridis(option="plasma") + DarkTheme(),
  
  FeaturePlot(NCG_germ_macaque3_temp,
              features = "prospermatogonia1", label = F, repel = FALSE) +
    scale_color_viridis(option="plasma") + DarkTheme(),
  
  FeaturePlot(NCG_germ_macaque3_temp,
              features = "spg1", label = F, repel = FALSE) +
    scale_color_viridis(option="plasma") + DarkTheme(),
  
  FeaturePlot(NCG_germ_macaque3_temp,
              features = "diffspgE1", label = F, repel = FALSE) +
    scale_color_viridis(option="plasma") + DarkTheme(),
  
  FeaturePlot(NCG_germ_macaque3_temp,
              features = "diffspgL1", label = F, repel = FALSE) +
    scale_color_viridis(option="plasma") + DarkTheme(),
  
  FeaturePlot(NCG_germ_macaque3_temp,
              features = "preleptotene1", label = F, repel = FALSE) +
    scale_color_viridis(option="plasma") + DarkTheme(),
  
  FeaturePlot(NCG_germ_macaque3_temp,
              features = "meiosis1", label = F, repel = FALSE) +
    scale_color_viridis(option="plasma") + DarkTheme()
  
)


FeaturePlot(NCG_germ_macaque3_temp,
            features = c("ASB9", "NANOS2"), label = F, repel = FALSE) &
  scale_color_viridis(option="plasma") & DarkTheme()


FeaturePlot(NCG_germ_macaque3_temp,
            features = c("CITED2", "CYP11A1"), label = F, repel = FALSE) &
  scale_color_viridis(option="plasma") & DarkTheme()


DimPlot(NCG_germ_macaque3_temp, label = T)
#remove the apoptotic clusters
# NCG_germ_macaque3_temp <- subset(NCG_germ_macaque3_temp, idents = c(15, 16, 17), invert = T)


table(Idents(NCG_germ_macaque3_temp))
#orange
DefaultAssay(NCG_germ_macaque3_temp) <- "RNA"
suppressMessages(require(DoubletFinder))

NCG_germ_macaque3_temp <- NormalizeData(NCG_germ_macaque3_temp)
NCG_germ_macaque3_temp = FindVariableFeatures(NCG_germ_macaque3_temp, verbose = F)
NCG_germ_macaque3_temp = ScaleData(NCG_germ_macaque3_temp, vars.to.regress = c("nFeature_RNA", "percent.mito"),
                                   verbose = F)
NCG_germ_macaque3_temp = RunPCA(NCG_germ_macaque3_temp, verbose = F, npcs = 25)
NCG_germ_macaque3_temp <- FindNeighbors(object = NCG_germ_macaque3_temp)

NCG_germ_macaque3_temp = FindClusters(NCG_germ_macaque3_temp)

NCG_germ_macaque3_temp = RunUMAP(NCG_germ_macaque3_temp, dims = 1:25, verbose = F)

nExp <- round(ncol(NCG_germ_macaque3_temp) * 0.04)  # expect 4% doublets
# nExp <- 1000 #increase this to account for closer to what clusters 8 2 and 4
NCG_germ_macaque3_temp <- doubletFinder_v3(NCG_germ_macaque3_temp, pN = 0.25, pK = 0.09, nExp = nExp, PCs = 1:10)

# name of the DF prediction can change, so extract the correct column name.
DF.name = colnames(NCG_germ_macaque3_temp@meta.data)[grepl("DF.classification", colnames(NCG_germ_macaque3_temp@meta.data))]


cowplot::plot_grid(ncol = 3, DimPlot(NCG_germ_macaque3_temp, group.by = "orig.ident", label = T) + NoAxes(),
                   DimPlot(NCG_germ_macaque3_temp, group.by = DF.name) + NoAxes(), #FeaturePlot(NCG_germ_macaque3_temp, features = "MEIOB", min.cutoff = 0, order = T) + NoAxes())
                   DimPlot(NCG_germ_macaque3_temp, group.by = "Phase"))

FeaturePlot(NCG_germ_macaque3_temp, features = "KIT", split.by = "orig.ident", min.cutoff = 0, order = T) & NoAxes() & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(NCG_germ_macaque3_temp, features = "MEIOB", min.cutoff = 0, order = T)
# 
NCG_germ_macaque3_temp2 <- NCG_germ_macaque3_temp
# 
# NCG_germ_macaque3_temp2$doublets <- NCG_germ_macaque3_temp$DF.classifications_0.25_0.09_461
# 
# DimPlot(NCG_germ_macaque3_temp2, group.by = "doublets")

# Idents(NCG_germ_macaque3_temp2) <- "doublets"

# NCG_germ_macaque3_singlets <- subset(NCG_germ_macaque3_temp2, idents = "Singlet")





# Idents(NCG_germ_macaque3_singlets) <- "seurat_clusters"
# DimPlot(NCG_germ_macaque3_singlets, label = T)
# VlnPlot(NCG_germ_macaque3_singlets, features = c("nFeature_RNA", "nCount_RNA", "percent.mito"))
# 
# 
# NCG_germ_macaque3_singlets <- subset(NCG_germ_macaque3_singlets,
#                              subset = nCount_RNA < 70000)


#SUBSET to remove doublet-heavy clusters##


save(NCG_germ_macaque3, file = "NCG_germ_macaque3.2023.12.21")

Idents(NCG_germ_macaque3) <- "seurat_clusters"
DimPlot(NCG_germ_macaque3, label = T, group.by = "seurat_clusters")
NCG_germ_macaque3_subset <- subset(NCG_germ_macaque3, idents = c(21, 20, 7), invert = T) #7
DimPlot(NCG_germ_macaque3_subset, label = T)



Idents(NCG_germ_macaque3_subset) <- "Phase"

NCG_S <- subset(NCG_germ_macaque3_subset, idents = "S")
NCG_G1 <- subset(NCG_germ_macaque3_subset, idents = "G1")
NCG_G2M <- subset(NCG_germ_macaque3_subset, idents = "G2M")

n.dims <- 35

current.list <- c(NCG_S, NCG_G1, NCG_G2M)
for (i in 1:length(x = current.list)) {
  # current.list[[i]] <- NormalizeData(object = current.list[[i]], verbose = FALSE)
  current.list[[i]] <- FindVariableFeatures(object = current.list[[i]], nfeatures = 2000, verbose = FALSE)
}

current.anchors <- FindIntegrationAnchors(object.list = current.list, dims = 1:n.dims)
# save(current.anchors, file="current.anchors.Robj")
current.integrated <- IntegrateData(anchorset = current.anchors, dims = 1:n.dims)
# save(current.integrated, file="current.integrated.Robj")


DefaultAssay(object = current.integrated) <- "integrated"
current.integrated <- ScaleData(object = current.integrated, vars.to.regress = c("nFeature_RNA", "percent.mito"), verbose = FALSE)
current.integrated <- RunPCA(object = current.integrated,
                             npcs = n.dims,
                             verbose = TRUE)
current.integrated <- FindNeighbors(object = current.integrated)
current.integrated <- FindClusters(object = current.integrated,
                                   resolution = 1.1,
                                   #algorithm = 1,
                                   verbose = TRUE
)
# current.integrated <- RunTSNE(object = current.integrated, dims = 1:20)
current.integrated <- RunUMAP(object = current.integrated, reduction = "pca",
                              dims = 1:n.dims)

DimPlot(current.integrated, group.by = "Phase")


NCG_germ_macaque4 <- current.integrated

save(NCG_germ_macaque4, file="NCG_germ_macaque4.2023.12.21.Robj")
load("NCG_germ_macaque4.2023.12.21.Robj")

gene.list <- c("POU5F1", "KIT", "REC8", "STRA8", "SYCP1", "MEIOB")

Keren.gene.list <- c("POU5F1", "NANOG", "TCL1A", "PDPN", "PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1",
                     "UTF1", "KIT", "ASB9", "NANOS2", "PIWIL4", "MORC1", "MAGEB2", "TDRD1", "RHOXF1", "MAGEC2", "ZBTB16", "EGR4", "RET", "MAGEA4",
                     "SOHLH1", "DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8", "DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9", "ZCWPW1", "RAD51AP2", "SYCP3",
                     "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1", "MLH3", "PIWIL1", "MLH1", "SPATA8", "CCDC112", "OVOL1", "OVOL2",
                     "CCDC42", "ACR", "CDRT15", "TEKT1", "H1FOO", "TJP3", "FAM24A", "CATSPER3", "SPACA3", "SPACA1", "ACRV1", "SAXO1", "CAPZA3", "HEMGN", "CA2", "PRM3",
                     "TNP2", "TSSK6", "SMCP", "SPEM2", "LELP1", "MKI67", "UBE2C", "TOP2A")


FeaturePlot(NCG_germ_macaque4, features = gene.list,reduction = "umap", #cols = c("gray18", "blue","lightblue", "white"), 
            min.cutoff = 0, order = T) & DarkTheme() & scale_color_viridis(option = "C")
p0 <- FeaturePlot(NCG_germ_macaque4, features = c("POU5F1"), cols = c("gray20", "bisque4", "bisque1", "yellow", "white"), order = T, min.cutoff = 0) & DarkTheme() #& scale_color_continuous(type = "gradient")
p1 <- FeaturePlot(NCG_germ_macaque4, split.by = "orig.ident", features = c("TFAP2C"), cols = c("gray20", "green"), order = T, min.cutoff = 0) & DarkTheme() #& scale_color_continuous(type = "gradient")
# p2 <- FeaturePlot(NCG_germ_macaque4, features = c("PDPN"), cols = c("gray20", "red", "orange", "yellow"), order = T, min.cutoff = 0) & DarkTheme() #& scale_color_continuous(type = "gradient")
p2 <- FeaturePlot(NCG_germ_macaque4, split.by = "orig.ident", features = c("DDX4"), cols = c("gray20", "red"), order = T, min.cutoff = 0) & DarkTheme() #& scale_color_continuous(type = "gradient")
p3 <- FeaturePlot(NCG_germ_macaque4, split.by = "orig.ident", features = c("PIWIL4"), cols = c("gray20", "cyan"), order = T, min.cutoff = 0) & DarkTheme() #& scale_color_continuous(type = "gradient")

plot_grid(p1, p2, p3, ncol = 1)

p4 <- FeaturePlot(NCG_germ_macaque4, features = c("MEIOB"), cols = c("gray20", "blue4", "magenta4", "violet", "white"), order = T, min.cutoff = 0) & DarkTheme() #& scale_color_continuous(type = "gradient")
p5 <- FeaturePlot(NCG_germ_macaque4, features = c("SYCP1"), cols = c("gray20", "saddlebrown", "tan", "white"), order = T, min.cutoff = 0) & DarkTheme() #& scale_color_continuous(type = "gradient")
p6 <- FeaturePlot(NCG_germ_macaque4, features = c("MAGEC2"), cols = c("gray20", "deeppink4", "lightpink", "white"), order = T, min.cutoff = 0) & DarkTheme() #& scale_color_continuous(type = "gradient")


FeaturePlot(NCG_germ_macaque4, features = c("ID3"), cols = c("gray20", "red", "orange", "yellow"), order = T, min.cutoff = 0) & DarkTheme() #& scale_color_continuous(type = "gradient")


plot_grid(p1, p2, p3, 
          p4, p5, p6, ncol = 3)


FeaturePlot(NCG_germ_macaque4, features = Keren.gene.list[1:18], order = T, min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(NCG_germ_macaque4, features = Keren.gene.list[19:36],  order = T, min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(NCG_germ_macaque4, features = Keren.gene.list[37:54],  order = T, min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")
FeaturePlot(NCG_germ_macaque4, features = Keren.gene.list[55:76],  order = T, min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")

list.of.genes <- c("DNMT3L", 
                   "UTF1", 
                   "PIWIL4", 
                   "ID4", 
                   "SOHLH1", 
                   "EGR4",
                   "FGFR3", 
                   "NANOS3", 
                   "UCHL1", 
                   "L1TD1", 
                   "ASB9",
                   "DMRT1",
                   "DMRTB1", 
                   "GCNA", 
                   "STRA8",
                   "KIT",
                   "REC8",
                   "SYCP3",
                   "PRDM9",
                   "SPO11")

list.of.genes2 <- c("SOHLH1", "EGR4", "UCHL1", "KIT","REC8", "STRA8", "SPO11", "MEIOB")

list.of.genes3 <- c("ID4", "ITGA6", "ZBTB16", "ETV5", "GFRA1", "DPPA4", "SALL4",
                    "UCHL1", "FOXO1", "MKI67", "DMRT1", "STRA8", "SOHLH2", "REC8", "SYCP2", "HORMAD1", "TEX101", "LY6K")


FeaturePlot(NCG_germ_macaque4, features = list.of.genes2,  order = T, min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")
p1 <- DimPlot(NCG_germ_macaque4, group.by = "seurat_clusters", label = T) + NoLegend()
p2 <- FeaturePlot(NCG_germ_macaque4, features = list.of.genes2[1],  order = T, min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")
p3 <- FeaturePlot(NCG_germ_macaque4, features = list.of.genes2[2],  order = T, min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")
p4 <- FeaturePlot(NCG_germ_macaque4, features = list.of.genes2[3],  order = T, min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")
p5 <- FeaturePlot(NCG_germ_macaque4, features = list.of.genes2[4],  order = T, min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")
p6 <- FeaturePlot(NCG_germ_macaque4, features = list.of.genes2[5],  order = T, min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")
p7 <- FeaturePlot(NCG_germ_macaque4, features = list.of.genes2[6],  order = T, min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")
p8 <- FeaturePlot(NCG_germ_macaque4, features = list.of.genes2[7],  order = T, min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")
p9 <- FeaturePlot(NCG_germ_macaque4, features = list.of.genes2[8],  order = T, min.cutoff = 0) & DarkTheme() & scale_color_viridis(option = "C")

plot_grid(p1, p2, p3, 
          p4, p5, p6, p7, p8, p9)

NCG_germ_macaque4$original.seurat_clusters <- NCG_germ_macaque4$seurat_clusters


Idents(NCG_germ_macaque4) <- "seurat_clusters"
DimPlot(NCG_germ_macaque4, label = T)


NCG_germ_macaque4[["UMAP"]] <- NCG_germ_macaque4[["umap"]]


NCG_germ_macaque4subset <- subset(NCG_germ_macaque4, idents = "transplant")
Idents(NCG_germ_macaque4subset) <- "seurat_clusters"

DimPlot(NCG_germ_macaque4subset, reduction = "UMAP", label = T)
FeaturePlot(NCG_germ_macaque4subset, features = gene.list, reduction = "UMAP", order = T)

save(NCG_germ_macaque4, file = "NCG_germ_macaque4_newCellNames_noPseudotime.2023.07.27.Robj")
# save(NCG_germ_macaque4, file = "NCG_germ_macaque4_oldCellNames.2023.08.04.Robj")



###PSEUDOTIME STARTS HERE###


cds <- as.cell_data_set(NCG_germ_macaque4)


cds <- cluster_cells(cds, reduction_method = c("UMAP"), #cluster_method = "louvain",
                     resolution = 2e-4,
                     k = 8,
                     partition_qval = 0.05
)

plot_cells(cds, show_trajectory_graph = FALSE)
cds <- learn_graph(cds,
                   close_loop = F,
                   learn_graph_control=list(ncenter=400, minimal_branch_len=10), #250, 18 is pretty good
                   use_partition = F)
print(plot_cells(cds, label_groups_by_cluster = FALSE, label_leaves = FALSE, label_branch_points = FALSE)) #+ ggtitle(a))

# Figure 6B ###
cds <- order_cells(cds)#, root_cells = max.avp)
p1 <- plot_cells(cds, color_cells_by = "pseudotime", label_cell_groups = FALSE, label_leaves = FALSE,
                 label_branch_points = FALSE)+ NoLegend()

p2 <- DimPlot(NCG_germ_macaque4, group.by = "cell.type") + NoLegend()

plot_grid(p2, p1)

pseudotime.values <- pseudotime(cds, reduction_method = "UMAP")

NCG_germ_macaque4[["pseudotime"]] <- pseudotime.values



DefaultAssay(NCG_germ_macaque4) <- "RNA"
macaque.assay.data <- GetAssayData(NCG_germ_macaque4)f
data.frame.of.interest <- NCG_germ_macaque4[["pseudotime"]]
# gene.of.interest <- c("NANOS3", "EGR4")

# values.of.interest <- t(as.data.frame(macaque.assay.data[gene.of.interest,]))
values.of.interest <- t(as.data.frame(macaque.assay.data))
data.frame.of.interest <- cbind(data.frame.of.interest, values.of.interest)

ggplot(data.frame.of.interest, aes(x = pseudotime, y = ACRV1)) + geom_point()



#redo module scores


gene.list.PGC.E <- c("POU5F1", "NANOS3", "MIF", "ACTG2")
gene.list.PGC.L <- c("SOX15", "L1TD1", "PDPN", "TUBA1C")
gene.list.M <- c("KIT", "CCND2", "EEF1G", "CCPG1" )
gene.list.M2T1 <- c("MAP4K4", "ASB9", "PRXL2A", "TLE5")
gene.list.T1 <- c("DCAF4L1", "ID1", "ID3", "MORC1")
gene.list.spermatogonia <- c("UTF1", "PAGE2", "RHOXF1", "MAGEB2")
gene.list.diff.spermatgonia.E <- c("SOHLH1", "KCNQ2", "MAGEC2", "MAGEA4")
gene.list.diff.spermatgonia.L   <- c("DMRT1", "SOHLH2", "ESX1", "CTCFL")
gene.list.preleptotene <- c("RAD51AP2", "SCML1", "TEX19", "PRDM9")
gene.list.leptotene <- c("TEX101", "SYCP1", "DMRTC2", "TDRG1")
gene.list.zygotene <- c("CTAGE1", "EQTN", "SPACA1", "GOLGA6L2")


NCG_germ_macaque4 <- AddModuleScore(NCG_germ_macaque4,
                                    features = list(apoptotic.genes),
                                    name="Apoptotic")
NCG_germ_macaque4 <- AddModuleScore(NCG_germ_macaque4,
                                    features = list(gene.list.PGC.E),
                                    name="PGCE")
NCG_germ_macaque4 <- AddModuleScore(NCG_germ_macaque4,
                                    features = list(gene.list.PGC.L),
                                    name="PGCL")
NCG_germ_macaque4 <- AddModuleScore(NCG_germ_macaque4,
                                    features = list(gene.list.M),
                                    name="M")
NCG_germ_macaque4 <- AddModuleScore(NCG_germ_macaque4,
                                    features = list(gene.list.M2T1),
                                    name="M2T")
NCG_germ_macaque4 <- AddModuleScore(NCG_germ_macaque4,
                                    features = list(gene.list.T1),
                                    name="T")
NCG_germ_macaque4 <- AddModuleScore(NCG_germ_macaque4,
                                    features = list(gene.list.spermatogonia),
                                    name="spg")
NCG_germ_macaque4 <- AddModuleScore(NCG_germ_macaque4,
                                    features = list(gene.list.diff.spermatgonia.E),
                                    name="diffspgE")
NCG_germ_macaque4 <- AddModuleScore(NCG_germ_macaque4,
                                    features = list(gene.list.diff.spermatgonia.L),
                                    name="diffspgL")
NCG_germ_macaque4 <- AddModuleScore(NCG_germ_macaque4,
                                    features = list(gene.list.preleptotene),
                                    name="preleptotene")
NCG_germ_macaque4 <- AddModuleScore(NCG_germ_macaque4,
                                    features = list(gene.list.leptotene),
                                    name="leptotene")
NCG_germ_macaque4 <- AddModuleScore(NCG_germ_macaque4,
                                    features = list(gene.list.zygotene),
                                    name="zygotene")


# scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))

# Plot scores

plot_grid(
  
  FeaturePlot(NCG_germ_macaque4,
              features = "Apoptotic1",min.cutoff = 0,max.cutoff = 2, order = T, label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ_macaque4,
              features = "PGCE1",min.cutoff = 0, max.cutoff = 2,order = T, label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ_macaque4,
              features = "PGCL1",min.cutoff = 0, max.cutoff = 2,order = T, label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ_macaque4,
              features = "M1", min.cutoff = 0, max.cutoff = 2,order = T,label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  FeaturePlot(NCG_germ_macaque4,
              features = "M2T1", min.cutoff = 0, max.cutoff = 2,order = T,label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  FeaturePlot(NCG_germ_macaque4,
              features = "T1", min.cutoff = 0, max.cutoff = 2,order = T,label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ_macaque4,
              features = "spg1",min.cutoff = 0, max.cutoff = 2,order = T, label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ_macaque4,
              features = "diffspgE1",min.cutoff = 0, max.cutoff = 2,order = T, label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ_macaque4,
              features = "diffspgL1",min.cutoff = 0, max.cutoff = 2,order = T, label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ_macaque4,
              features = "preleptotene1", min.cutoff = 0,max.cutoff = 2, order = T,label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ_macaque4,
              features = "leptotene1", min.cutoff = 0,max.cutoff = 2,order = T,label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme(),
  
  FeaturePlot(NCG_germ_macaque4,
              features = "zygotene1", min.cutoff = 0, max.cutoff = 2,order = T,label = F, repel = FALSE) +
    scale_color_viridis(option="plasma", label = function(x) sprintf("%.2f", x)) + DarkTheme()
  
)



Idents(NCG_germ_macaque4) <- "seurat_clusters"
DimPlot(NCG_germ_macaque4, label = T)
new.cluster.ids <- character(length(levels(Idents(NCG_germ_macaque4))))


PGC.E.ids <- c(16, 12)+1
PGC.L.ids <- c(0,5,8)+1
M.ids <- c(9,4,11)+1
M2T1.ids <- c(7)+1
T1.ids <- c(14,2,6,13)+1
spermatogonia.ids <- c(1,3,10,15)+1


new.cluster.ids[PGC.E.ids] <- "PGCs (early)"
new.cluster.ids[PGC.L.ids] <- "PGCs (late)"
new.cluster.ids[M.ids] <- "M prospermatogonia"
new.cluster.ids[M2T1.ids] <- "M_T1 prospermatogonia"
new.cluster.ids[T1.ids] <- "T1 prospermatogonia"
new.cluster.ids[spermatogonia.ids] <- "Spermatogonia"


new.cluster.ids

names(x = new.cluster.ids) <- levels(x = NCG_germ_macaque4)
NCG_germ_macaque4 <- RenameIdents(object = NCG_germ_macaque4, new.cluster.ids)
NCG_germ_macaque4$cell.type <- Idents(NCG_germ_macaque4)
DimPlot(NCG_germ_macaque4)
NCG_germ_macaque4$cell.type <- factor(x = NCG_germ_macaque4$cell.type, levels = c(
  "PGCs (early)", 
  "PGCs (late)",
  "M prospermatogonia",
  "M_T1 prospermatogonia",
  "T1 prospermatogonia",
  "Spermatogonia"
  
))
Idents(NCG_germ_macaque4) <- NCG_germ_macaque4$cell.type
#Plot out designations with a feature of interest.
DimPlot(NCG_germ_macaque4, label = F, reduction = "umap", group.by = "cell.type",cols = DiscretePalette(9, palette = "glasbey"))

setwd("~/Documents/scRNAseq/Rscripts/savefiles")
save(NCG_germ_macaque4, file="NCG_germ_macaque4_2024.01.03.Robj")
# load("NCG_germ_macaque4_newCellNames_2023.07.27.Robj")



#HEATMAP
setwd("~/Documents/scRNAseq/Rscripts/savefiles")
load("NCG_germ_macaque4_2024.01.03.Robj")
NCG_germ_macaque4_subset <- NCG_germ_macaque4
Idents(NCG_germ_macaque4_subset) <- "cell.type"
DimPlot(NCG_germ_macaque4_subset)

# all.markers.subset <- filter(all.markers_NCG, avg_log2FC > 1)

PGCs.E.markers.all_macaque <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "PGCs (early)")
PGCs.L.markers.all_macaque <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "PGCs (late)")
M.ProSpermatogonia.markers.all_macaque <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "M prospermatogonia")
M2T1.markers.all_macaque <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "M_T1 prospermatogonia")
T.ProSpermatogonia.markers.all_macaque <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "T1 prospermatogonia")
Spg.markers.all_macaque <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "Spermatogonia")


# setwd("/Users/ewhelan/Desktop/NCG/NCG_marker_genes")
# write.csv(PGCs.E.markers.all, "macaque_PGCs.E.markers.all.csv")
# write.csv(PGCs.L.markers.all, "macaque_PGCs.L.markers.all.csv")
# write.csv(M.ProSpermatogonia.markers.all, "macaque_M.ProSpermatogonia.markers.all.csv")
# write.csv(M2T1.markers.all, "macaque_M2T1.markers.all.csv")
# write.csv(T.ProSpermatogonia.markers.all, "macaque_T.ProSpermatogonia.markers.all.csv")
# write.csv(Spg.markers.all, "macaque_Spg.markers.all.csv")
# # write.csv(Diff.Spg.E.markers.all, "Diff.Spg.E.markers.all.csv")
# # write.csv(Diff.Spg.L.markers.all, "Diff.Spg.L.markers.all.csv")
# # write.csv(Pre.Leptotene.markers.all, "Pre.Leptotene.markers.all.csv")

PGCs.E.markers.all_macaque$gene <- rownames(PGCs.E.markers.all_macaque)
PGCs.L.markers.all_macaque$gene <- rownames(PGCs.L.markers.all_macaque)
M.ProSpermatogonia.markers.all_macaque$gene <- rownames(M.ProSpermatogonia.markers.all_macaque)
M2T1.markers.all_macaque$gene <- rownames(M2T1.markers.all_macaque)
T.ProSpermatogonia.markers.all_macaque$gene <- rownames(T.ProSpermatogonia.markers.all_macaque)
Spg.markers.all_macaque$gene <- rownames(Spg.markers.all_macaque)



log2FC_cutoff <- 1

PGCs.E.markers.ordered <- PGCs.E.markers.all_macaque %>% filter(p_val_adj < 0.05 & avg_log2FC > log2FC_cutoff )
PGCs.L.markers.ordered <- PGCs.L.markers.all_macaque %>% filter(p_val_adj < 0.05 & avg_log2FC > log2FC_cutoff )
M.ProSpermatogonia.markers.ordered <- M.ProSpermatogonia.markers.all %>% filter(p_val_adj < 0.05 & avg_log2FC > log2FC_cutoff )
M2T1.markers.ordered <- M2T1.markers.all_macaque %>% filter(p_val_adj < 0.05 & avg_log2FC > log2FC_cutoff )
T.ProSpermatogonia.markers.ordered <- T.ProSpermatogonia.markers.all_macaque %>% filter(p_val_adj < 0.05 & avg_log2FC > log2FC_cutoff )
Spg.markers.ordered <- Spg.markers.all_macaque %>% filter(p_val_adj < 0.05 & avg_log2FC > log2FC_cutoff)


genes.per.cluster <- 3000
PGCs.E.markers.ordered <- PGCs.E.markers.ordered %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
PGCs.L.markers.ordered <- PGCs.L.markers.ordered %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
M.ProSpermatogonia.markers.ordered <- M.ProSpermatogonia.markers.ordered %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
M2T1.markers.ordered <- M2T1.markers.ordered %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
T.ProSpermatogonia.markers.ordered <- T.ProSpermatogonia.markers.ordered %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
Spg.markers.ordered <- Spg.markers.ordered %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)


PGCs.E.markers.ordered$cell.type <- "PGCs.E" 
PGCs.L.markers.ordered$cell.type <- "PGCs.L" 
M.ProSpermatogonia.markers.ordered$cell.type <- "M" 
M2T1.markers.ordered$cell.type <- "M2T1" 
T.ProSpermatogonia.markers.ordered$cell.type <- "T"
Spg.markers.ordered$cell.type <- "Spg"


topNmarkers <- rbind(PGCs.E.markers.ordered, PGCs.L.markers.ordered, M.ProSpermatogonia.markers.ordered, M2T1.markers.ordered,
                     T.ProSpermatogonia.markers.ordered, Spg.markers.ordered) #, #Diff.Spg.E.markers.ordered, Diff.Spg.L.markers.ordered, 
# Pre.Leptotene.markers.ordered)

# write.csv(topNmarkers, file=paste0("top", genes.per.cluster, "markers_macaque.csv"))
# write.csv(topNmarkers, file="top3000genes_macaque_with_cutoff_0.585.csv")

dim(topNmarkers)


# Heatmap of select genes by cluster ----

Idents(NCG_germ_macaque4_subset) <- NCG_germ_macaque4_subset$cell.type
DimPlot(NCG_germ_macaque4_subset)


#fix the typo in M/T1 propermaogonia
# levels(NCG_germ_macaque4_subset$cell.type)
# levels(NCG_germ_macaque4_subset$cell.type)[levels(NCG_germ_macaque4_subset$cell.type) == "M/T1 propermaogonia"] <- "M/T1 prospermatogonia"
# levels(NCG_germ_macaque4_subset$cell.type)
# Idents(NCG_germ_macaque4_subset) <- NCG_germ_macaque4_subset[["cell.type"]]


data.pgcE <- as.matrix(GetAssayData(NCG_germ_macaque4_subset, slot = "data")[, WhichCells(NCG_germ_macaque4_subset, ident = "PGCs (early)")])
data.pgcL <- as.matrix(GetAssayData(NCG_germ_macaque4_subset, slot = "data")[, WhichCells(NCG_germ_macaque4_subset, ident = "PGCs (late)")])
data.m <- as.matrix(GetAssayData(NCG_germ_macaque4_subset, slot = "data")[, WhichCells(NCG_germ_macaque4_subset, ident = "M prospermatogonia")])
data.mt1 <- as.matrix(GetAssayData(NCG_germ_macaque4_subset, slot = "data")[, WhichCells(NCG_germ_macaque4_subset, ident = "M_T1 prospermatogonia")])
data.t1 <- as.matrix(GetAssayData(NCG_germ_macaque4_subset, slot = "data")[, WhichCells(NCG_germ_macaque4_subset, ident = "T1 prospermatogonia")])
data.undiff <- as.matrix(GetAssayData(NCG_germ_macaque4_subset, slot = "data")[, WhichCells(NCG_germ_macaque4_subset, ident = "Spermatogonia")])


av.counts.pgcE <- apply(data.pgcE, 1, mean)
av.counts.pgcL <- apply(data.pgcL, 1, mean)
av.counts.m <- apply(data.m, 1, mean)
av.counts.mt1 <- apply(data.mt1, 1, mean)
av.counts.t1 <- apply(data.t1, 1, mean)
av.counts.undiff <- apply(data.undiff, 1, mean)


av.counts.HS_NCG_only_NCG <- cbind(av.counts.pgcE, av.counts.pgcL, av.counts.m, av.counts.mt1,av.counts.t1,
                                   av.counts.undiff )



av.counts.HS_NCG_only_NCG.df <- as.data.frame(av.counts.HS_NCG_only_NCG)

gene.name.list <- rownames(av.counts.HS_NCG_only_NCG.df)

av.counts.HS_NCG_only_NCG.df$gene.name <- gene.name.list

# genes.of.interest <- Keren.gene.list[c(1:48)] #this is just Keren's list of genes
# genes.of.interest <- (all.markers.subset$gene) #this is a log-fold cutoff

genes.of.interest <- rownames(topNmarkers) #this is just taking the top N markers
# genes.of.interest <- (topNmarkers) 

av.counts.genes.of.interest <- av.counts.HS_NCG_only_NCG.df %>%
  dplyr::filter(gene.name %in% genes.of.interest)


av.counts.genes.of.interest.ordered <- av.counts.genes.of.interest[match(genes.of.interest, av.counts.genes.of.interest$gene.name),]

integer.col <- ncol(av.counts.genes.of.interest.ordered)

av.counts.genes.of.interest.ordered <- av.counts.genes.of.interest.ordered[,-integer.col]

av.counts.genes.of.interest.matrix1 <- as.matrix(av.counts.genes.of.interest.ordered)
colnames(av.counts.genes.of.interest.matrix1) <- c("PGCs (early)", "PGCs (late)", "M prospermatogonia", 
                                                   "M_T1 prospermatogonia", "T1 prospermatogonia", "Spermatogonia")
# "Diff. spermatogonia (early)", "Diff. spermatogonia (late)", 
# "Preleptotene spermatocytes")

av.counts.genes.of.interest.matrix.no.na <- na.omit(av.counts.genes.of.interest.matrix1)


l.mat <- rbind(c(0,4), c(0, 3), c(2,1), c(0,0))
l.mat

# Figure 3A ----

#NOTE select matrix1 or matrix2 below to change between unselected and EpCAM


# Assuming av.counts.genes.of.interest.matrix.no.na is your data frame

# Keep rows where at least one value is above 1
av.counts.genes.of.interest.matrix.filtered <- av.counts.genes.of.interest.matrix.no.na[rowSums(av.counts.genes.of.interest.matrix.no.na > 0.1) > 0, ]

dim(av.counts.genes.of.interest.matrix.no.na)
dim(av.counts.genes.of.interest.matrix.filtered)
table(av.counts.genes.of.interest.matrix.filtered)


#Figure 6 F

# tiff(file = "plot.tiff",width = 5000, height = 100000)
dev.off()
heatmap.2(av.counts.genes.of.interest.matrix.no.na, Colv = FALSE,
          Rowv = FALSE,
          trace = "none",
          # col =  colorRampPalette(colors=c("blue","white","red"))(100),
          col = inferno(100),
          density.info="none",
          srtCol = -60, #angle of labels
          adjCol = c(0,0), #offset of labels, should be centered now
          cexRow=1.0, cexCol=1.3,
          lmat= l.mat,
          lwid = c(0.5, 5),
          lhei=c(1, 1, 5, 1),
          scale = "row"
          
)

dev.off()


genes <- c("POU5F1", "NANOG", "TCL1A", "PDPN", "PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1",
           "UTF1", "KIT", "ASB9", "NANOS2", "PIWIL4", "MORC1", "MAGEB2", "TDRD1", "RHOXF1", "MAGEC2", "ZBTB16", "EGR4", "RET", "MAGEA4",
           "SOHLH1", "DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8", "DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9", "ZCWPW1", "RAD51AP2", "SYCP3",
           "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1")

DimPlot(HS_NCG_only_NCG)
NCG_germ_macaque4_subset <- subset(HS_NCG_only_NCG, idents = c("Pachytene spermatocytes", "Diplotene/2ndary spermatocytes", "Spermatids"), invert = T)
DimPlot(NCG_germ_macaque4_subset)
# all.markers <- FindAllMarkers(NCG_germ_macaque4_subset, min.pct = 0, logfc.threshold = 0, features = genes)
PGCs.E.markers <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "PGCs (early)", min.pct = 0, logfc.threshold = 0, features = genes)
PGCs.L.markers <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "PGCs (late)", min.pct = 0, logfc.threshold = 0, features = genes)
M.ProSpermatogonia.markers <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "M prospermatogonia", min.pct = 0, logfc.threshold = 0, features = genes)
M2T1.markers <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "M/T1 prospermatogonia", min.pct = 0, logfc.threshold = 0, features = genes)
T.ProSpermatogonia.markers <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "T1 prospermatogonia", min.pct = 0, logfc.threshold = 0, features = genes)
Spg.markers <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "Spermatogonia", min.pct = 0, logfc.threshold = 0, features = genes)
Diff.Spg.E.markers <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "Diff. spermatogonia (early)", min.pct = 0, logfc.threshold = 0, features = genes)
Diff.Spg.L.markers <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "Diff. spermatogonia (late)", min.pct = 0, logfc.threshold = 0, features = genes)
Pre.Leptotene.markers <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "Prelep/Lep/Zygotene", min.pct = 0, logfc.threshold = 0, features = genes)


PGCs.E.markers.ordered <- PGCs.E.markers %>% 
  arrange(avg_log2FC) %>% 
  slice_head(n = 50)


graph.data <- graph.data %>%
  arrange(Celltype)

order.of.genes <- as.vector(genes)


# graph.data$Pathway2 <- factor(graph.data$Pathway, levels = rev(unique(graph.data$Pathway)))
all.markers$gene <- factor(all.markers$gene, levels = rev(genes))
all.markers$cluster <- factor(all.markers$cluster, levels = c("PGCs (early)", 
                                                              "PGCs (late)", 
                                                              "M prospermatogonia",
                                                              "M/T1 prospermaogonia", 
                                                              "T1 prospermatogonia", 
                                                              "Spermatogonia",
                                                              "Diff. spermatogonia (early)", 
                                                              "Diff. spermatogonia (late)", 
                                                              "Prelep/Lep/Zygotene"))



ggplot(data = all.markers, na.rm = T,
       aes(x = cluster, y = (gene), 
           color = avg_log2FC, size = pct.1))+ 
  geom_point() +
  # scale_color_gradient(high = "red" , low = "blue") +
  scale_color_viridis(option = "inferno")+
  scale_size(range = c(-7, 7)) +
  theme_bw() + 
  ylab("Gene") + 
  xlab("Cell Type") + 
  ggtitle("Genes expressed by cell type") +
  scale_y_discrete(position = "right")+
  theme(axis.text.x = element_text( #face = "bold",
    color = "black", size = 14, angle = -60, hjust = 0
  ),
  axis.text.y = element_text( #face = "bold",
    color = "black",  hjust = 1)
  ) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(legend.position = "left")




#pairwise comparison



all.markers$cluster <- factor(all.markers$cluster, levels = c("PGCs (early)", 
                                                              "PGCs (late)", 
                                                              "M prospermatogonia",
                                                              "M/T1 prospermaogonia", 
                                                              "T1 prospermatogonia", 
                                                              "Spermatogonia",
                                                              "Diff. spermatogonia (early)", 
                                                              "Diff. spermatogonia (late)", 
                                                              "Prelep/Lep/Zygotene"))

###

p.val.cutoff <- 0.05
log2fc.cutoff <- 0.4

###

# NCG_germ_macaque4_subset_original -> NCG_germ_macaque4_subset

# NCG_germ_macaque4_subset <- NCG_germ_macaque4_subset


PGCE.vs.PGCL <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "PGCs (early)", ident.2 = "PGCs (late)", logfc.threshold = log2fc.cutoff, min.pct = 0)
PGCL.vs.M <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "PGCs (late)", ident.2 = "M prospermatogonia", logfc.threshold = log2fc.cutoff, min.pct = 0)
M.vs.MT1 <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "M prospermatogonia", ident.2 = "M/T1 prospermatogonia", logfc.threshold = log2fc.cutoff, min.pct = 0)
MT1.vs.T1 <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "M_T1 prospermatogonia", ident.2 = "T1 prospermatogonia", logfc.threshold = log2fc.cutoff, min.pct = 0)
T1.vs.Spg <- FindMarkers(NCG_germ_macaque4_subset, ident.1 = "T1 prospermatogonia", ident.2 = "Spermatogonia", logfc.threshold = log2fc.cutoff, min.pct = 0)


list.of.comparisons <- list(PGCE.vs.PGCL, PGCL.vs.M, M.vs.MT1, MT1.vs.T1, T1.vs.Spg)
list.of.counts <- list(av.counts.pgcE, av.counts.pgcL, av.counts.m, av.counts.mt1, av.counts.t1, av.counts.undiff)
list.of.cell.types <- c("PGCs (early)", 
                        "PGCs (late)", 
                        "M prospermatogonia",
                        "M/T1 prospermaogonia", 
                        "T1 prospermatogonia", 
                        "Spermatogonia")

for(number in 1:(length(list.of.cell.types)-1)){
  # number <- 1
  current.comparison <- list.of.comparisons[[number]]
  current.counts <- list.of.counts[[number]]
  current.counts2 <- list.of.counts[[number+1]]
  
  current.comparison$gene <- rownames(current.comparison)
  # current.comparison.filtered <- filter(current.comparison, p_val_adj < p.val.cutoff)
  current.comparison.filtered.up <- filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC > 0)
  current.comparison.filtered.down <- filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC < 0)
  current.comparison.dataframe <- data.frame(current.counts, current.counts2)
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  current.comparison.subset <- filter(current.comparison.dataframe, av.counts.diffL > 0 | av.counts.Spermatocytes>0)
  # current.comparison.sig.up <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.up$gene)
  # current.comparison.sig.down <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.down$gene)
  current.comparison.subset <- na.omit(current.comparison.subset)
  current.comparison.subset$DEGs <- "Not significant"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene] <- "Up"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene] <- "Down"
  print(paste0(list.of.cell.types[number], " vs ", list.of.cell.types[number+1], 
               " Up: ", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]), 
               " Down:", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])))
  
  
  assign(paste0("p", number), (ggplot() + geom_point(data = current.comparison.subset, aes(x = current.counts, y = current.counts2, color = DEGs)) + 
                                 scale_color_manual(values=c('#E69F00','#999999', '#56B4E9')) + 
                                 theme_classic() + labs(title=paste(list.of.cell.types[number], " vs. ", list.of.cell.types[number+1], ", q", p.val.cutoff, ", log2fc", log2fc.cutoff),
                                                        x=list.of.cell.types[number], y = list.of.cell.types[number+1]) + theme(plot.title = element_text(size=9)))
  )
}


plot_grid(ncol = 4, p1, p2, p3, p4, 
          p5, p6, p7, p8)





gene.of.interest <- "STRA8"

gene.of.interest.PGCE <- c("SPRR2F", "ABHD12B", "HMGA2", "DDT", "KHDC3L", "NANOG", "SOX17", "YBX1", "MIF")
gene.of.interest.PGCL <- c("GAGE1", "TEX14", "GAGE12H", "GAGE2A", "XIST", "DCAF4L1", "ELAVL2", "TFRAP2C", "UTF1",
                           "CCND1", "ID3", "MAGEB2", "FHL2", "ID1", "MAGED2", "TFRAP2C", "SOX17", "NANOG", "UTF1", "POU5F1")
gene.of.interest.M <- c("GCNA", "HIST1H3B", "FOXM1", "SYCP2", "E2F8", "CDC25A", "XIST", "MKI67",
                        "ACTG2", "DPPA5", "UBE2C", "KHDC3L", "GPHA2", "CENPA", "NANOG", "APOC1", "POU5F1")
gene.of.interest.transitional <- c("NANOS2", "PAGE2", "RHOXF1", "ASB9", "TEX15", "TDRD1", "PAGE2B",
                                   "DGKK", "NANOS3", "MDK", "MAP4K4", "STMN1", "IGFBP2")
gene.of.interest.T1 <- c("FOXF1", "DPPA5", "EGR4", "SIX1", "PIWIL4", "VCX3B", "MAGEC2", "DPPA2", "MORC1", "MAGEB2",
                         "GAGE12H", "NANOS2", "MEG3", "EGR4", "XIST", "HELLPAR", "DCAF4L1")
gene.of.interest.spg <- c("RET", "POU3F1", "HMGB4", "TSPAN33", "FOXD1", "SCT", "KCNQ2", "MAGEA4", "VCX", "TCF3", "UTF1", "PIWIL4", 
                          "UTF1", "MORC1", "SCT", "FOXD1", "POU3F1", "RGS14", "TSPAN33", "PIWIL4", "EGR4", "MAPK3", "HES5", "PHOX2A", "FIGLA")
gene.of.interest.earlydiff <- c("HIST1H1A", "L1TD1", "ASB9", "EIF3K", "EGR1", "HIST1H3B", "HMGB4",
                                "RESF1", "ID3", "DUSP6", "EGR4", "SIX1", "DPPA2", "FGFR3", "ID4", "ID2", "PCSK1N", "ASB9", "MAGEA4")
gene.of.interest.latediff <- c("TEX19", "MEIOC", "SYCP3", "SYCE2", "RAD51AP2", "HORMAD1", "SCML1", "CENPU", "TOP2A", "SMC1B", "HSPA5", "HIST1H1A",
                               "ELAVL2", "CENPF", "TDRD1", "ESX1", "ITGA6", "RHOXF1", "KIT", "UBE2C", "HIST1H2AB", "SOHLH1", "CCND1")
gene.of.interest.prelep <- c("TEX101", "MEIOB", "SPO11", "ZCWPW1", "RAD51AP2", "SYCP1", "SYCP2", "TEX12", "CCNB3", "RNF212B", "FBXO47", "HOTAIR")

list.of.lists.of.genes <- list(gene.of.interest.PGCE,
                               gene.of.interest.PGCL,
                               gene.of.interest.M,
                               gene.of.interest.transitional,
                               gene.of.interest.T1,
                               gene.of.interest.spg,
                               gene.of.interest.earlydiff,
                               gene.of.interest.latediff,
                               gene.of.interest.prelep)




for(number in 1:(length(list.of.cell.types)-1)){
  # number <- 1
  print(number)
  gene.of.interest <- c(list.of.lists.of.genes[[number]], list.of.lists.of.genes[[number+1]])
  current.comparison <- list.of.comparisons[[number]]
  current.counts <- list.of.counts[[number]]
  current.counts2 <- list.of.counts[[number+1]]
  
  current.comparison$gene <- rownames(current.comparison)
  # current.comparison.filtered <- filter(current.comparison, p_val_adj < p.val.cutoff)
  current.comparison.filtered.up <- filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC > 0)
  current.comparison.filtered.down <- filter(current.comparison, p_val_adj < p.val.cutoff & avg_log2FC < 0)
  current.comparison.dataframe <- data.frame(current.counts, current.counts2)
  current.comparison.dataframe$gene.name <- rownames(current.comparison.dataframe)
  current.comparison.subset <- filter(current.comparison.dataframe, av.counts.diffL > 0 | av.counts.Spermatocytes>0)
  # current.comparison.sig.up <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.up$gene)
  # current.comparison.sig.down <- filter(current.comparison.dataframe, gene.name %in% current.comparison.filtered.down$gene)
  current.comparison.subset <- na.omit(current.comparison.subset)
  current.comparison.subset$DEGs <- "Not significant"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene] <- "Up"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene] <- "Down"
  current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% gene.of.interest] <- "GeneOfInterest"
  print(paste0(list.of.cell.types[number], " vs ", list.of.cell.types[number+1], 
               " Up: ", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.up$gene]), 
               " Down:", length(current.comparison.subset$DEGs[current.comparison.subset$gene.name %in% current.comparison.filtered.down$gene])))
  
  colnames(current.comparison.subset) <- c("HS.counts", "NCG.counts", "gene.name", "DEGs")
  
  
  degs_order <- c( "Down", "GeneOfInterest", "Not significant", "Up")
  
  current.comparison.subset_ordered <- current.comparison.subset %>%
    mutate(DEGs = factor(DEGs, levels = degs_order)) 
  
  # Define the order of levels for DEGs
  degs_order <- c("Down", "GeneOfInterest", "Not significant", "Up")
  
  # Create the scatterplot with "GeneOfInterest" points on top
  assign(paste0("r", number), ggplot(data = current.comparison.subset_ordered, aes(x = HS.counts, y = NCG.counts, 
                                                                                   label = gene.name,
                                                                                   color = DEGs)) + 
           geom_point(data = filter(current.comparison.subset_ordered, DEGs == "Not significant"), alpha = 0.5) + # Plot "Not Significant" points first
           geom_point(data = filter(current.comparison.subset_ordered, DEGs != c("Not significant", "GeneOfInterest"))) + # Plot other points
           geom_point(data = filter(current.comparison.subset_ordered, DEGs == "GeneOfInterest")) + 
           geom_text(data = filter(current.comparison.subset_ordered, DEGs == "GeneOfInterest"), hjust = -0.1, vjust = -0.1) +
           scale_color_manual(values = c('#E69F00', "red4", '#999999', '#56B4E9')) + 
           theme_classic() + 
           labs(title = paste(list.of.cell.types[number], " vs. ", list.of.cell.types[number + 1], ", q", p.val.cutoff, ", log2fc", log2fc.cutoff),
                x = list.of.cell.types[number], y = list.of.cell.types[number + 1]) + 
           theme(plot.title = element_text(size = 9)))
}




plot_grid(ncol = 1, r1, r2, r3, r4, 
          r5, r6, r7, r8)




### heatmap ###


genes <- c("POU5F1", "NANOG", "TCL1A", "PDPN", "PRDM1", "NANOS3", "TFAP2C", "SOX17", "L1TD1",
           "UTF1", "KIT", "ASB9", "NANOS2", "PIWIL4", "MORC1", "MAGEB2", "TDRD1", "RHOXF1", "MAGEC2", "ZBTB16", "EGR4", "RET", "MAGEA4",
           "SOHLH1", "DMRT1", "SOHLH2", "ESX1", "DMRTB1", "STRA8", "DAZL", "SCML2", "GAL3ST1", "TEX19", "PRDM9", "ZCWPW1", "RAD51AP2", "SYCP3",
           "DMC1", "PEG10", "SPO11", "MEIOB", "ATR", "SYCP1", "MYBL1", "HORMAD1")

genes <- rownames(topNmarkers)

variable.features <- VariableFeatures(NCG_germ_macaque4_subset)

# "MLH3", "PIWIL1", "MLH1", "SPATA8", "CCDC112", "OVOL1", "OVOL2",
# "CCDC42", "ACR", "CDRT15", "TEKT1", "H1FOO", "TJP3", "FAM24A", "CATSPER3", "SPACA3", "SPACA1", "ACRV1", "SAXO1", "CAPZA3", "HEMGN", "CA2", "PRM3",
# "TNP2", "TSSK6", "SMCP", "SPEM2", "LELP1", "MKI67", "UBE2C", "TOP2A")

HS_NCG_only_NCG2 <- NCG_germ_macaque4_subset

DefaultAssay(HS_NCG_only_NCG2) <- "RNA"


#can I do this with seurat objects? Need to pull out the normalized counts

pt.matrix.all <- HS_NCG_only_NCG2@assays$RNA@data
pseudotime <- as.data.frame(HS_NCG_only_NCG2$pseudotime)

dim(pt.matrix.all)
# gene.names.all <- rownames(pt.matrix.all)

# length(genes)
# length(intersect(genes, gene.names.all))
# setdiff(genes, gene.names.all)

pt.matrix.all[1:4, 1:4]




intersect.genes <- (intersect(genes, rownames(pt.matrix.all)))
genes <- intersect.genes
# genes <- variable.features

pt.matrix.subset <- pt.matrix.all[genes,]
t.pt.matrix.all <- t(pt.matrix.subset)
t.pt.matrix.all[1:4, 1:4]


pseudotime[1:4,1]
t.pt.matrix.all <- as.data.frame(t.pt.matrix.all)
dim(t.pt.matrix.all)
t.pt.matrix.all$pseudotime <- pseudotime[,1]
dim(t.pt.matrix.all)
t.pt.matrix.all[1:5, 1:5]
t.pt.matrix.sorted <- t.pt.matrix.all[order(t.pt.matrix.all$pseudotime),]
t.pt.matrix.sorted[1:5, ]


# pt.matrix <- t(as.matrix(t.pt.matrix.sorted[,1:(length(genes))]))
pt.matrix <- t(as.matrix(t.pt.matrix.sorted))
pt.matrix[, 1:5]

pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=5)$y}))
# pt.matrix <- t(apply(pt.matrix,1,function(x){smooth.spline(x,df=2)$y}))
pt.matrix <- t(apply(pt.matrix,1,function(x){(x-mean(x))/sd(x)}))
dim(pt.matrix)
pt.matrix <- pt.matrix[c(1:length(genes)),]
dim(pt.matrix)

#Ward.D2 Hierarchical Clustering
hthc <- Heatmap(
  pt.matrix,
  name                         = "z-score",
  col                          = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
  show_row_names               = TRUE,
  show_column_names            = FALSE,
  row_names_gp                 = gpar(fontsize = 16),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  row_title_rot                = 0,
  use_raster                   = FALSE,
  cluster_rows                 = TRUE,
  cluster_row_slices           = TRUE,
  cluster_columns              = FALSE)

print(hthc)

VlnPlot(NCG_germ_macaque4, features = c("NANOS2", "ID1", "MEIOB"))



#plot pseudotime values by sample

pseudotime.values <- NCG_germ_macaque4$pseudotime
range((pseudotime.values), na.rm = TRUE)

orig.ident <- NCG_germ_macaque4$orig.ident
stage.values <- NCG_germ_macaque4$stage

pseudotime.df <- data.frame(stage.values, pseudotime.values)

stripplot(pseudotime.values)



table(stage.values)


ggplot(pseudotime.df, aes(x=pseudotime.values, y=0, color = stage.values)) +
  geom_point(size = 4)  +
  # annotate("segment",x=1,xend=20, y=0, yend=0, size=2) +
  # annotate("segment",x=1,xend=1, y=-0.1,yend=0.1, size=2) +
  # annotate("segment",x=20,xend=20, y=-0.1,yend=0.1, size=2) +
  # geom_text(aes(label = x), col="white") +
  scale_x_continuous(limits = c(1,max(pseudotime.df$pseudotime.values))) +
  scale_y_continuous(limits = c(-1,1)) +
  # scale_color_manual(values = unname(colours)) #+ 
  theme(panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank())

q_colors =  12
v_colors =  viridis(q_colors, option = "B")
show_col(v_colors)

stripchart(pseudotime.values ~ stage.values, data = pseudotime.df, col = v_colors[1:9], 
           pch = 16, frame = FALSE, method = "jitter")


hist(pseudotime.df$pseudotime.values, breaks = 100)

#Try with density lines

table(pseudotime.df$stage.values)




DimPlot(NCG_germ_macaque4, split.by = "stage", group.by = "cell.type") &
  # scale_color_viridis(option="plasma") & 
  DarkTheme()


#work out the proportions by sample. 
DimPlot(NCG_germ_macaque4, group.by = "stage")
individual.samples <- SplitObject(NCG_germ_macaque4, split.by = "stage")

for(current.loop in 1:length(individual.samples)){
  # current.loop <- 1
  Idents(individual.samples[[current.loop]]) <- "stage"
  print(table(Idents(individual.samples[[current.loop]])))
  Idents(individual.samples[[current.loop]]) <- "cell.type"
  print(table(Idents(individual.samples[[current.loop]])))
}




# Figure showing proportions of samples.


DimPlot(NCG_germ_macaque4_noPGCLC)
DimPlot(NCG_germ_macaque4_transplant, split.by = "stage")



FeaturePlot(NCG_germ_macaque4_noPGCLC, features = c("TFAP2C", "DDX4", "KIT"), ncol = 3)  & NoLegend() & NoAxes() & DarkTheme() & scale_color_viridis(option = "magma")

DimPlot(NCG_germ_macaque4, label = T, group.by = "cell.type") + NoLegend() + NoAxes()
DimPlot(NCG_germ_macaque4, label = T, group.by = "cell.type")+ NoAxes()




Idents(NCG_germ_macaque4) <- "stage"
DimPlot(NCG_germ_macaque4)
NCG_germ_macaque4_4m <- subset(NCG_germ_macaque4, idents = c("d123", "d127", "d120"))
NCG_germ_macaque4_6m <- subset(NCG_germ_macaque4, idents = c("d180", "d181"))
NCG_germ_macaque4_8m <- subset(NCG_germ_macaque4, idents = c("d251", "d245"))
Idents(NCG_germ_macaque4_4m) <- "cell.type"
Idents(NCG_germ_macaque4_6m) <- "cell.type"
Idents(NCG_germ_macaque4_8m) <- "cell.type"
p1 <- DimPlot(NCG_germ_macaque4_4m) + NoAxes() + NoLegend()
p2 <- DimPlot(NCG_germ_macaque4_6m) + NoAxes() + NoLegend()
p3 <- DimPlot(NCG_germ_macaque4_8m) + NoAxes() + NoLegend()
plot_grid(p1, p2, p3, ncol = 3)




#repeat for cluster heatmap

NCG_germ_macaque4b <- NCG_germ_macaque3_temp
Idents(NCG_germ_macaque4b) <- "cell.type"
DimPlot(NCG_germ_macaque4b)
DimPlot(NCG_germ_macaque4b, group.by = "seurat_clusters", label = T)
NCG_germ_macaque4b_subset <- subset(NCG_germ_macaque4b, idents = c("PGCLCs"), invert = T)
NCG_germ_macaque4b_subset <- NCG_germ_macaque4b
Idents(NCG_germ_macaque4b_subset) <- "seurat_clusters"
DimPlot(NCG_germ_macaque4b_subset, label = T)

order.of.clusters <- c(8, 4, 11, 5, 3, 1, 12, 6, 9, 0, 2, 7, 14, 16, 17)


# all.markers.subset <- filter(all.markers_NCG, avg_log2FC > 1)
cluster0.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 0)
cluster1.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 1)
cluster2.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 2)
cluster3.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 3)
cluster4.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 4)
cluster5.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 5)
cluster6.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 6)
cluster7.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 7)
cluster8.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 8)
cluster9.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 9)
cluster10.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 10) 
cluster11.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 11)
cluster12.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 12)
cluster13.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 13)
cluster14.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 14)
cluster15.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 15)
cluster16.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 16)
cluster17.markers <-FindMarkers(NCG_germ_macaque4b_subset, ident.1 = 17)


# setwd("/Users/ewhelan/Desktop/NCG/NCG_marker_genes")
# write.csv(PGCs.E.markers.all, "PGCs.E.markers.all.csv")



genes.per.cluster <- 100


cluster0.markers.ordered <- cluster0.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster1.markers.ordered <- cluster1.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster2.markers.ordered <- cluster2.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster3.markers.ordered <- cluster3.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster4.markers.ordered <- cluster4.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster5.markers.ordered <- cluster5.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster6.markers.ordered <- cluster6.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster7.markers.ordered <- cluster7.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster8.markers.ordered <- cluster8.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster9.markers.ordered <- cluster9.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster10.markers.ordered <- cluster10.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster11.markers.ordered <- cluster11.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster12.markers.ordered <- cluster12.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster13.markers.ordered <- cluster13.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster14.markers.ordered <- cluster14.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster15.markers.ordered <- cluster15.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster16.markers.ordered <- cluster16.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)
cluster17.markers.ordered <- cluster17.markers %>% arrange(avg_log2FC) %>% slice_tail(n = genes.per.cluster)

cluster0.markers.ordered$cluster <- 0
cluster1.markers.ordered$cluster <- 1
cluster2.markers.ordered$cluster <- 2
cluster3.markers.ordered$cluster <- 3
cluster4.markers.ordered$cluster <- 4
cluster5.markers.ordered$cluster <- 5
cluster6.markers.ordered$cluster <- 6
cluster7.markers.ordered$cluster <- 7
cluster8.markers.ordered$cluster <- 8
cluster9.markers.ordered$cluster <- 9
cluster10.markers.ordered$cluster <- 10
cluster11.markers.ordered$cluster <- 11
cluster12.markers.ordered$cluster <- 12
cluster13.markers.ordered$cluster <- 13
cluster14.markers.ordered$cluster <- 14
cluster15.markers.ordered$cluster <- 15
cluster16.markers.ordered$cluster <- 16
cluster17.markers.ordered$cluster <- 17


topNmarkers <- rbind(cluster10.markers.ordered, cluster15.markers.ordered, cluster13.markers.ordered,
                     cluster8.markers.ordered, cluster4.markers.ordered, cluster11.markers.ordered, cluster5.markers.ordered, 
                     cluster3.markers.ordered, cluster1.markers.ordered,cluster12.markers.ordered, cluster6.markers.ordered, 
                     cluster9.markers.ordered, cluster0.markers.ordered,  cluster2.markers.ordered, cluster7.markers.ordered, 
                     cluster14.markers.ordered, cluster16.markers.ordered, cluster17.markers.ordered)







# Heatmap of select genes by cluster ----

Idents(NCG_germ_macaque4b_subset) <- NCG_germ_macaque4b_subset$seurat_clusters
DimPlot(NCG_germ_macaque4b_subset)

data.0 <- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 0)])
data.1 <- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 1)])
data.2 <- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 2)])
data.3 <- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 3)])
data.4 <- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 4)])
data.5 <- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 5)])
data.6 <- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 6)])
data.7 <- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 7)])
data.8 <- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 8)])
data.9 <- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 9)])
data.10 <- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 10)])
data.11<- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 11)])
data.12 <- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 12)])
data.13 <- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 13)])
data.14 <- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 14)])
data.15 <- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 15)])
data.16 <- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 16)])
data.17 <- as.matrix(GetAssayData(NCG_germ_macaque4b_subset, slot = "data")[, WhichCells(NCG_germ_macaque4b_subset, ident = 17)])

av.counts.0 <- apply(data.0, 1, mean)
av.counts.1 <- apply(data.1, 1, mean)
av.counts.2 <- apply(data.2, 1, mean)
av.counts.3 <- apply(data.3, 1, mean)
av.counts.4 <- apply(data.4, 1, mean)
av.counts.5 <- apply(data.5, 1, mean)
av.counts.6 <- apply(data.6, 1, mean)
av.counts.7 <- apply(data.7, 1, mean)
av.counts.8 <- apply(data.8, 1, mean)
av.counts.9 <- apply(data.9, 1, mean)
av.counts.10 <- apply(data.10, 1, mean)
av.counts.11 <- apply(data.11, 1, mean)
av.counts.12 <- apply(data.12, 1, mean)
av.counts.13 <- apply(data.13, 1, mean)
av.counts.14 <- apply(data.14, 1, mean)
av.counts.15 <- apply(data.15, 1, mean)
av.counts.16 <- apply(data.16, 1, mean)
av.counts.17 <- apply(data.17, 1, mean)

av.counts.HS_NCG_only_NCG <- cbind(av.counts.10, av.counts.15, av.counts.13, 
                                   av.counts.8, av.counts.4,av.counts.11, av.counts.5, av.counts.3,
                                   av.counts.1, av.counts.12, av.counts.6,av.counts.9, 
                                   av.counts.0, av.counts.2, 
                                   av.counts.7, 
                                   av.counts.14, 
                                   av.counts.16, av.counts.17)


av.counts.HS_NCG_only_NCG.df <- as.data.frame(av.counts.HS_NCG_only_NCG)

gene.name.list <- rownames(av.counts.HS_NCG_only_NCG.df)

av.counts.HS_NCG_only_NCG.df$gene.name <- gene.name.list

# genes.of.interest <- Keren.gene.list[c(1:48)] #this is just Keren's list of genes
# genes.of.interest <- (all.markers.subset$gene) #this is a log-fold cutoff
genes.of.interest <- rownames(topNmarkers) #this is just taking the top N markers

av.counts.genes.of.interest <- av.counts.HS_NCG_only_NCG.df %>%
  dplyr::filter(gene.name %in% genes.of.interest)


av.counts.genes.of.interest.ordered <- av.counts.genes.of.interest[match(genes.of.interest, av.counts.genes.of.interest$gene.name),]

integer.col <- ncol(av.counts.genes.of.interest.ordered)

av.counts.genes.of.interest.ordered <- av.counts.genes.of.interest.ordered[,-integer.col]

av.counts.genes.of.interest.matrix1 <- as.matrix(av.counts.genes.of.interest.ordered)
colnames(av.counts.genes.of.interest.matrix1) <- c("cluster10", "cluster15", "cluster13", 
                                                   
                                                   "cluster8", "cluster4", "cluster11", "cluster5", "cluster3", "cluster1", "cluster12", 
                                                   "cluster6", "cluster9", "cluster0", "cluster2", "cluster7", "cluster14", "cluster16", "cluster17")

# c("cluster0", "cluster1", "cluster2", "cluster3", "cluster4", "cluster5", "cluster6", 
#                                                  "cluster7", "cluster8", "cluster9", #"cluster10", 
#                                                  "cluster11", "cluster12", 
#                                                  # "cluster13",
#                                                  "cluster14", 
#                                                  # "cluster15", 
#                                                  "cluster16", "cluster17")





# av.counts.genes.of.interest.matrix1 <- av.counts.genes.of.interest.matrix1[, order.of.clusters2]

av.counts.genes.of.interest.matrix.no.na <- na.omit(av.counts.genes.of.interest.matrix1)

l.mat <- rbind(c(0,4), c(0, 3), c(2,1), c(0,0))
l.mat

# Figure 3A ----

#NOTE select matrix1 or matrix2 below to change between unselected and EpCAM


# tiff(file = "plot.tiff",width = 5000, height = 100000)

dev.off()
heatmap.2(av.counts.genes.of.interest.matrix.no.na, Colv = FALSE,
          Rowv = FALSE,
          trace = "none",
          # col =  colorRampPalette(colors=c("blue","white","red"))(100),
          
          # col = rev(brewer.pal(n = 11, name = "Spectral")),
          # col = rev(brewer.pal(n = 11, name = "RdYlGn")),
          # col = rev(brewer.pal(n = 11, name = "RdYlBu")),
          # col = rev(brewer.pal(n = 11, name = "RdGy")),
          col = rev(brewer.pal(n = 11, name = "RdBu")),
          # col = rev(brewer.pal(n = 11, name = "PuOr")),
          # col = rev(brewer.pal(n = 11, name = "PRGn")),
          # col = rev(brewer.pal(n = 11, name = "PiYG")),
          # col = rev(brewer.pal(n = 11, name = "BrBG")),
          # col = inferno(100),
          # col = viridis(100),
          density.info="none",
          srtCol = -90, #angle of labels
          adjCol = c(0,0), #offset of labels, should be centered now
          key = T,
          cexRow=1.0, cexCol=1.3,
          lmat= l.mat,
          lwid = c(0.5, 5),
          lhei=c(1, 1, 5, 1),
          scale = "row"
          
)



## plot using seurat's default tools
DimPlot(NCG_germ_macaque4 , reduction = "PHATE", group.by = "orig.ident", label = F) + ggtitle(label = "PHATE")+ ggtitle(paste("knn ", knn.value, ", decay ", decay.value, sep = ""))
DimPlot(NCG_germ_macaque4 , reduction = "PHATE",  group.by = "cell.type", label = T) + ggtitle(label = "PHATE")+ ggtitle(paste("knn ", knn.value, ", decay ", decay.value, sep = ""))


FeaturePlot(NCG_germ_macaque4, features = c("SCML2", "SALL4", "REC8", "SPDYA", "PRDM9", "MAGEC2", "ESX1", "MAGEA4", "UCHL1", "DMRT1", "SOHLH2", "SYCP2", "ZBTB16", "ETV5", "MEIOB", "HORMAD1", "TEX101", "LY6K", "FOXO1"), 
            # split.by = "orig.ident", 
            order = T) & DarkTheme() & scale_color_viridis(option = "C")
DimPlot(NCG_germ_macaque4, group.by = "cell.type", split.by = "orig.ident", label = F)
DimPlot(NCG_germ_macaque4, group.by = "seurat_clusters", split.by = "orig.ident", label = F)


FeaturePlot(NCG_germ_macaque4, features = c("POU5F1", "TFAP2C", "NANOS3", "PIWIL4", "REC8", "ASB9"), reduction = "PHATE",
            # split.by = "orig.ident", 
            order = T) & DarkTheme() & scale_color_viridis(option = "C") & xlim (c(-0.04, 0.04)) & ylim (c(-0.02, 0.04))


eoin.palette.light <- c("antiquewhite3",
                        "lightcyan3", "pink3", "brown",  "sienna1", "aquamarine3","dodgerblue1", "yellowgreen", "lightgoldenrod2","darkorchid1", "darkgreen", "darkblue")

eoin.palette.dark <- c("gray28",
                       "cyan4", "pink3", "brown",  "sienna3", "aquamarine3","dodgerblue3", "darkgreen", "mediumpurple3","darkorchid1", "yellowgreen", "lightgoldenrod2", "darkblue")

eoin.palette.dark2 <- c("antiquewhite3","gray28",
                        "cyan4", "pink3", "brown",  "sienna3", "aquamarine3","dodgerblue3", "darkgreen", "mediumpurple3","darkorchid1", "yellowgreen", "lightgoldenrod2", "darkblue")


cell.type.colors.full <- c("#2179b4", "#b4d88b", "#36a047", "#f6999a", "#f15a29", "#fdbf6f", "#f7941d", "#cab3d6", "#6b3f98",
                           "#8c3f30", "#cdc7b1", "#2a2d7c")



DimPlot(NCG_germ_macaque4, group.by = "seurat_clusters", split.by = "orig.ident", ncol = 2)


#Figure 6C
DimPlot(NCG_germ_macaque4, group.by = "cell.type", cols = cell.type.colors.full)






p2 <- FeaturePlot(NCG_germ_macaque4, #split.by = "orig.ident", 
                  features = c("DDX4"), order = T, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "red", "orangered", "yellow", "white"))
p1 <- FeaturePlot(NCG_germ_macaque4, #split.by = "orig.ident",
                  features = c("TFAP2C"), order = T, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "darkgreen", "green", "yellow", "white"))
p3 <- FeaturePlot(NCG_germ_macaque4,#split.by = "orig.ident", 
                  features = c("PIWIL4"), order = T, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "cyan4", "cyan", "white"))

### Figure 6D

plot <- plot_grid(p1, p2, p3, ncol = 1)

ggsave("macaque_TFAP2C_DDX4_PIWIL4.pdf", plot, width = 4, height = 8)
# pdf("macaque_TFAP2C_DDX4_PIWIL4.pdf", plot, width = 4, height = 8)


p4 <- FeaturePlot(HS_NCG, features = c("MAGEB2"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "blue2", "deepskyblue", "skyblue"))
p5 <- FeaturePlot(HS_NCG, features = c("MAGEA4"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "magenta4", "hotpink", "white"))
p6 <- FeaturePlot(HS_NCG, features = c("SYCP1"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "brown", "burlywood"))


p7 <- FeaturePlot(HS_NCG, features = c("MEIOC"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "orange4", "orange", "yellow"))
p8 <- FeaturePlot(HS_NCG, features = c("ACRV1"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "plum4", "plum", "white"))
p9 <- FeaturePlot(HS_NCG, features = c("PRM1"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "olivedrab", "olivedrab1", "yellow"))

p4 <- FeaturePlot(NCG_germ_macaque4,  features = c("POU5F1"), order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "blue2", "deepskyblue", "skyblue", "white"))
p5 <- FeaturePlot(NCG_germ_macaque4, features = c("KIT"), order = T, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "magenta4", "hotpink", "white"))
p6 <- FeaturePlot(NCG_germ_macaque4,features = c("DAZL"), order = T, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "brown", "burlywood"))
p7 <- FeaturePlot(NCG_germ_macaque4,features = c("SALL4"), order = T, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "orange4", "orange", "yellow"))
p8 <- FeaturePlot(NCG_germ_macaque4,features = c("EGR4"), order = T, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "plum4", "plum", "white"))
p9 <- FeaturePlot(NCG_germ_macaque4,features = c("SOHLH1"), order = T, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "olivedrab", "olivedrab1", "yellow"))

### FIGURE 6E ###
plot_grid(p4, p5, p6, p7, p8, p9, ncol = 6)

p8 <- FeaturePlot(NCG_germ4,features = c("DND1"), order = T, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "plum4", "plum", "white"))
p9 <- FeaturePlot(NCG_germ4,features = c("NANOS3"), order = T, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "olivedrab", "olivedrab1", "yellow"))
plot_grid(p8, p9)

p1 <- FeaturePlot(NCG_germ_macaque4,  features = c("STRA8"), split.by = "orig.ident", order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "blue2", "deepskyblue", "skyblue", "white"))
p2 <- FeaturePlot(NCG_germ_macaque4,  features = c("SOHLH2"), split.by = "orig.ident", order = F, min.cutoff = 0) & DarkTheme() & scale_colour_gradientn(colors = c("gray20", "magenta4", "hotpink", "white"))
plot_grid(p1, p2, ncol = 1)







