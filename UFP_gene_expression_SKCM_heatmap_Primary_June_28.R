# 
# install.packages("pheatmap")
# install.packages("d3heatmap")
 # install.packages("survival")
 # install.packages("survminer")
 library(ggpubr)
library(ggplot2)
library(RColorBrewer)
library(sva)
library(matrixStats)
library(data.table)
library(dplyr)
library(reshape)
library(pheatmap)
library(circlize)
library(ComplexHeatmap)
library(factoextra)
library(MASS)
library(mclust)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms

library(dendextend) # for comparing two dendrograms
library(survival)
library(survminer)
 library(tidyr)
 library(matrixStats)
##
## create moveme function to move columns
####  https://stackoverflow.com/questions/3369959/moving-columns-within-a-data-frame-without-retyping

moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}



##### SET WORKING DIRECTORY , use getwd to designate new variable'current_dir' to use in sub directory names

setwd("~/Desktop/R_proteastasis_scripts/TCGA_Data_analysis_21/SKCM")
current_dir<-getwd()

#load log expression data and add new columns to clinical and expression data ####

loadRData <- function(fileName){
  #loads an RData file, and returns it
  load(fileName)
  get(ls()[ls() != "fileName"])
}

#### UFP ####
UFP_exp_info <- loadRData('UFP_Gene_log_exp.RData')

### load table of PN gene expression and create data frame with just sample number and cluster number

PN_table <- read.csv('SKCM_exp_info_prim_April.csv')
primary_cluster_allocation_df <- subset(PN_table, select = c(2,9) )
#Add additional column with new cluster name then remove previous cluster column
primary_cluster_allocation_df$Cluster <-  sapply(primary_cluster_allocation_df$cluster2, function(x)
  ifelse (x == '1' ,"SKCM(P)_A",'SKCM(P)_B'))

primary_clusters_df <- subset(primary_cluster_allocation_df , select = c(1,3) )
primary_clusters_df <- data.frame(primary_clusters_df, row.names = 1)   

# Save cluster table

save(primary_clusters_df, file = 'Primary_SKCM_Samples_Cluster_Allocation.Rdata')

########### Create table of just primary sample UFP expression ####
# # select just primary samples 

UFP_exp_info_prim <-  subset(UFP_exp_info, rownames(UFP_exp_info) %in% rownames(primary_clusters_df))

### Calculate Z scores of expression ###

# UFP_exp_info_prim_Z <- t(UFP_exp_info_prim [1:ncol(UFP_exp_info_prim)])
#(UFP_exp_info_prim_Z-rowMeans(UFP_exp_info_prim_Z))/(rowSds(as.matrix(UFP_exp_info_prim_Z)))[row(UFP_exp_info_prim_Z)]

UFP_exp_info_prim_t <- t(UFP_exp_info_prim)
UFP_exp_info_prim_Z <- apply(UFP_exp_info_prim_t[, 1:103], 1, function(x) (x - mean(x)) / sd(x))

zscore <- function(x) {
  z <- (x - mean(x)) / sd(x)
  return(z)
}

UFP_exp_info_prim_Z_2 <- function(UFP_exp_info_prim_t) {
  z <- (UFP_exp_info_prim_t - mean(UFP_exp_info_prim_t)) / sd(UFP_exp_info_prim_t)
  return(z)
}

UFP_exp_info_prim_Z_3 <- (scale(t(UFP_exp_info_prim_t)))
UFP_exp_info_prim_Z_3<- UFP_exp_info_prim_Z_3[ ,colSums(is.na(UFP_exp_info_prim_Z_3)) == 0]

#remove rows with NA
UFP_exp_info_prim_Z<- UFP_exp_info_prim_Z[ ,colSums(is.na(UFP_exp_info_prim_Z)) == 0]

UFP_exp_info_prim_Z_2<- UFP_exp_info_prim_Z_2[ ,colSums(is.na(UFP_exp_info_prim_Z_2)) == 0]

# # prob not needed Set Sample_Cluster as Factor ####
# 
# SKCM_exp_info_prim$Cluster <- factor(SKCM_exp_info_prim$Cluster,
#                                            levels = c( "SKCM(P)_A",'SKCM(P)_B'))
# save(UFP_exp_info_prim, file ='UFP_exp_info_prim_.RData')
# #

# 
# #### GENE CLUSTERS 
# 
# # Allocate Cluster numbers with hclust agglomerative
# 
# # methods to assess
# m <- c( "average", "single", "complete", "ward")
# names(m) <- c( "average", "single", "complete", "ward")
# 
# # function to compute coefficient
# ac <- function(x) {
#   agnes(heatmap_matrix_all, method = x)$ac
# }
# 
# map_dbl(m, ac)
# 
# 
# ### ward method gives highest clustering coefficient (0.935)
# 
# hc3 <- agnes(heatmap_matrix_all, method = "ward")
# pltree(hc3, cex = 0.1, hang = -1, main = "Dendrogram of Genes")
# 
# # Cut tree into 3 groups
# sub_grp <- cutree(hc3, k = 3)
# 
# # Number of members in each cluster
# table(sub_grp)
# 
# #Add gene cluster number 
# heatmap_matrix_df <-heatmap_matrix_all
# heatmap_matrix_df$Hcluster3  <- sub_grp 
# 
# # convert cluster column to data frame
# gene_cluster_df <- as.data.frame(heatmap_matrix_df$Hcluster3)
# 
# #Rename cluster column
# colnames(gene_cluster_df)[colnames(gene_cluster_df)=="heatmap_matrix_df$Hcluster3"] <- "Ward_Cluster"
# #Convert Row Names to Column
# gene_cluster_df$Gene <- rownames(gene_cluster_df)

### Read in PN Gene list

# 
# PN_Genes <- read.csv('PN_Genes.csv')
# PN_Genes <- as.data.frame(PN_Genes)
# 
# #Merge PN Gene Data Frame and   Cluster Data Frame 
# 
# PN_Genes_Cluster_df <- merge(x = PN_Genes,y = gene_cluster_df, by.x ="Gene.Symbol", by.y="Gene"  )
# 
# 
# 
# # Add new column with Ward cluster renumbered so order is same as heatmap cluster
# 
# PN_Genes_Cluster_df$Ward_Cluster_reordered <- sapply(PN_Genes_Cluster_df$Ward_Cluster, function(x)
#   ifelse (x == '3','Primary_1',
#           ifelse     (x == '1','Primary_2', 'Primary_3')
#   ))
# 
# # Rename PN_Genes_Cluster_df$Ward_Cluster_reordered  PN_Genes_Cluster_df$Gene_Cluster
# colnames(PN_Genes_Cluster_df)[colnames(PN_Genes_Cluster_df)=="Ward_Cluster_reordered"] <- "Primary_Gene_Cluster"
# 
# 
# # specify PN_Genes_Cluster_df$Classification as a factor
# PN_Genes_Cluster_df$Classification <- factor(PN_Genes_Cluster_df$Classification,
#                                                levels = c("Autophagy", "Degradation", "Deubiquitination","Ubiquitination","Folding" ))
# 
# 
# 
# 
# #Save PN_Genes_Cluster_df
# 
# save(PN_Genes_Cluster_df, file = 'SKCM_gene_clusters_May.RData')
# 
# 
# # load Gene clusters df  
# PN_Genes_Cluster_df <- loadRData('SKCM_gene_clusters_May.RData')

# # create new column with new names for primary clusters
# PN_Genes_Cluster_df$Primary_Gene_Cluster_new <- sapply(PN_Genes_Cluster_df$Primary_Gene_Cluster, function(x)
#   ifelse (x == 'Primary_1','PN_1(P)',
#           ifelse     (x == 'Primary_2','PN_2(P)', 'PN_3(P)')
#   ))
# 
# #Save PN_Genes_Cluster_df
# 
# save(PN_Genes_Cluster_df, file = 'SKCM_gene_clusters_June.RData')
# 
# # load Gene clusters df  
# PN_Genes_Cluster_df <- loadRData('SKCM_gene_clusters_June.RData')

######### HEATMAPS ######### 



#create matrix for heatmap (gene expression columns transposed so genes are rows and samples are columns )
heatmap_matrix_all <- as.matrix(t(UFP_exp_info_prim_Z_3))


# ________________________________________
#Create Colour function to distribute colours more evenly ####

heatcolS <- circlize::colorRamp2 (breaks = c(-7,-0.1,0,0.1, 2.5,5, 20),
                                  c("dodgerblue4","#f1fbfe" ,"white", "#ffff55", "#ff4015", "#ff0000", '#5e0406'),
                                  transparency = .3)
heatcolS(seq(-7, 20, 0.01))

# Create top annotation

ann2 <- data.frame (primary_clusters_df$Cluster)
                    

colnames(ann2) <- c('Sample_Cluster_Primary')

                    ann2$Sample_Cluster_Primary <- factor(ann2$Sample_Cluster_Primary,
                                                                levels = c( "SKCM(P)_A","SKCM(P)_B"))
                    
                  

colours2 <- list('Sample_Cluster_Primary' = c('SKCM(P)_A' = 'darkorchid4', 'SKCM(P)_B' = 'gold')
                             
)


colAnn2 <- HeatmapAnnotation(df = ann2,
                             which = 'col',
                             col = colours2,
                             annotation_width = unit(c(2, 10), 'cm'),
                             show_annotation_name = FALSE,
                             annotation_legend_param = list(
                               Sample_Cluster_Primary = list(
                                 title = "Sample Cluster",
                             #  title_position = "lefttop-rot",
                               annotation_legend_side="right",
                               annotation_legend_height = unit(160, "mm"),  
                               annotation_legend_width = unit(60, "mm"), 
                             grid_height = unit(2, "cm"), grid_width = unit(0.5, "cm"),
                             labels_gp = gpar(fontsize = 20),
                             title_gp = gpar(fontsize = 20, fontface = 'bold'))
                               )
                             )
                             

# rowAnn <- data.frame (PN_Genes_Cluster_df$Gene.Symbol ,PN_Genes_Cluster_df$Primary_Gene_Cluster_new)
# rowAnn   <- data.frame(rowAnn, row.names = 1)   
# colnames(rowAnn) <- c('Gene_Cluster_Primary')
# colours3 <- list('Gene_Cluster_Primary' = c('PN_1(P)' = 'slateblue', 'PN_2(P)' = 'pink', 'PN_3(P)' = 'darkcyan'))
# 
# rowAnn2 <- HeatmapAnnotation(df = rowAnn,
#                              which = 'row',
#                              col = colours3,
#                              annotation_width = unit(c(2, 10), 'cm'),
#                              show_annotation_name = FALSE,
#                              annotation_legend_param = list(
#                                Gene_Cluster_Primary = list(
#                               #   title_position = "lefttop-rot",
#                                  title = "Gene Cluster",
#                                annotation_legend_side="right",
#                                annotation_legend_height = unit(160, "mm"),  
#                                annotation_legend_width = unit(60, "mm"),
#                                grid_height = unit(2, "cm"), grid_width = unit(0.5, "cm"),
#                                labels_gp = gpar(fontsize = 20),
#                                title_gp = gpar(fontsize = 20, fontface = 'bold'))
#                              ))
# #
## Heatmap with clustering ward.D2 no names ####

pdf (paste0("UFP_heatmap_cancer_SKCM_3.pdf"),width=17, height=7)
#HM <-
Heatmap (heatmap_matrix_all, width = unit(260, "mm"), height = unit(150, "mm"),
        # name ='',
         rect_gp = gpar(col= FALSE),
         
      #   column_title_gp = gpar(fontsize = 20, fontface = "bold"),
         # column_names_centered = TRUE,
         #column_names_rot = 90,
         clustering_method_rows = "ward.D2",
         # clustering_method_columns = "ward.D2",
         # column_names_gp = gpar(fontsize = 25, fontface = "bold"),
         #row_names_gp = gpar(fontsize = 1, fontface = "bold"),
         row_names_gp = gpar(fontsize = 5, fontface = "bold"),
         column_dend_height = unit(10,'mm'),
         row_dend_width=unit(10,"mm"),
         #column_title = "Expression of Proteostasis Network Genes in Primary Skin Cutaneous Melanoma Samples",
         #column_title = "A",
         show_column_names = FALSE,
         show_row_names = TRUE,
         # row_split = 3, 
         column_split = ann2$Sample_Cluster_Primary,
         row_gap = unit(3, "mm"),
         column_gap = unit(3, "mm"),
         col =  heatcolS,
        top_annotation=colAnn2,
       #  right_annotation = rowAnn2 ,
         cluster_column_slices = FALSE,
      
      
         heatmap_legend_param = list(
           title = "Expression", 
          # title_position = "lefttop-rot",
           color_bar = "continuous",
           legend_height = unit(120, "mm"),  
           legend_width = unit(20, "mm"), 
           # fonts = 20, 
           # fontface = "bold",
          title_gp = gpar(fontsize = 20, fontface = 'bold'),
         labels_gp = gpar(col = "black", fontsize = 20),
           heatmap_legend_side="right")
      
     )

dev.off()

# # Final heatmap with clustered genes with gene names ####
# 
# 
# pdf (paste0("SKCM_heatmap_cancer_exp_Z_Score_all_primary_ward_cluster_with_names.pdf"),width=20, height=10)
# #HM <-
# Heatmap (heatmap_matrix_all, width = unit(260, "mm"), height = unit(150, "mm"),
#          name ='',
#          rect_gp = gpar(col= FALSE),
#          
#          column_title_gp = gpar(fontsize = 20, fontface = "bold"),
#          # column_names_centered = TRUE,
#          #column_names_rot = 90,
#          clustering_method_rows = "ward.D2",
#          clustering_method_columns = "ward.D2",
#          # column_names_gp = gpar(fontsize = 25, fontface = "bold"),
#          row_names_gp = gpar(fontsize = 1, fontface = "bold"),
#          # row_names_gp = gpar(fontsize = 15, fontface = "bold"),
#          column_dend_height = unit(10,'mm'),
#          row_dend_width=unit(10,"mm"),
#         # column_title = "Expression of Proteostasis Network Genes in Primary Skin Cutaneous Melanoma Samples",
#         column_title = "B",
#         show_column_names = FALSE,
#          show_row_names = FALSE,
#          row_split = 3, 
#          column_split = ann2$Sample_Cluster,
#          row_gap = unit(3, "mm"),
#          column_gap = unit(3, "mm"),
#          col =  heatcolS,
#          top_annotation=colAnn2,
#          right_annotation = rowAnn2 ,
#          cluster_column_slices = FALSE,
#          heatmap_legend_param = list(
#            title = "Expression", color_bar = "continuous",
#            legend_height = unit(40, "mm"),  legend_width = unit(40, "mm"), fonts = 10, fontface = "bold",
#            labels_gp = gpar(col = "black", fontsize = 10),
#            heatmap_legend_side="left", annotation_legend_side="right",
#            annotation_legend_height = unit(100, "mm"),  
#            annotation_legend_width = unit(60, "mm"), fonts = 20, fontface = "bold"))
# 
# dev.off()
# 
# # Final heatmap with clustered genes without gene names ####
# 
# 
# pdf (paste0("SKCM_heatmap_cancer_exp_Z_Score_all_primary_ward_cluster_without_names.pdf"),width=15, height=9)
# #HM <-
# Heatmap (heatmap_matrix_all, width = unit(260, "mm"), height = unit(150, "mm"),
#          name ='',
#          rect_gp = gpar(col= FALSE),
#          
#          column_title_gp = gpar(fontsize = 20, fontface = "bold"),
#          # column_names_centered = TRUE,
#          #column_names_rot = 90,
#          clustering_method_rows = "ward.D2",
#          clustering_method_columns = "ward.D2",
#          # column_names_gp = gpar(fontsize = 25, fontface = "bold"),
#          row_names_gp = gpar(fontsize = 1, fontface = "bold"),
#          # row_names_gp = gpar(fontsize = 15, fontface = "bold"),
#          column_dend_height = unit(10,'mm'),
#          row_dend_width=unit(10,"mm"),
#          column_title = "Expression of Proteostasis Network Genes in Primary Skin Cutaneous Melanoma Samples",
#          show_column_names = FALSE,
#          show_row_names = FALSE,
#          row_split = 3, 
#          column_split = ann2$Sample_Cluster,
#          row_gap = unit(3, "mm"),
#          column_gap = unit(3, "mm"),
#          col =  heatcolS,
#          top_annotation=colAnn2,
#          right_annotation = rowAnn2 ,
#          cluster_column_slices = FALSE,
#          heatmap_legend_param = list(
#            title = "Expression", color_bar = "continuous",
#            legend_height = unit(60, "mm"),  legend_width = unit(60, "mm"), fonts = 10, fontface = "bold",
#            labels_gp = gpar(col = "black", fontsize = 10),
#            heatmap_legend_side="left", annotation_legend_side="right",
#            annotation_legend_height = unit(300, "mm"),  
#            annotation_legend_width = unit(120, "mm"), fonts = 30, fontface = "bold"))
# 
# dev.off()
# 
# 
# 
# 


#### MFP ####
MFP_exp_info <- loadRData('MFP_Gene_log_exp.Rdata')


########### Create table of just primary sample MFP_exp_info_prim expression ####
# # select just primary samples 

MFP_exp_info_prim <-  subset(MFP_exp_info, rownames(MFP_exp_info) %in% rownames(primary_clusters_df))

### Calculate Z scores of expression ###

# MFP_exp_info_prim_Z <- t(MFP_exp_info_prim [1:ncol(MFP_exp_info_prim)])
#(MFP_exp_info_prim_Z-rowMeans(MFP_exp_info_prim_Z))/(rowSds(as.matrix(MFP_exp_info_prim_Z)))[row(MFP_exp_info_prim_Z)]

MFP_exp_info_prim_t <- t(MFP_exp_info_prim)


MFP_exp_info_prim_Z<- (scale(t(MFP_exp_info_prim_t)))

#remove rows with NA
MFP_exp_info_prim_Z<- MFP_exp_info_prim_Z[ ,colSums(is.na(MFP_exp_info_prim_Z)) == 0]


# # prob not needed Set Sample_Cluster as Factor ####
# 
# SKCM_exp_info_prim$Cluster <- factor(SKCM_exp_info_prim$Cluster,
#                                            levels = c( "SKCM(P)_A",'SKCM(P)_B'))
# save(MFP_exp_info_prim, file ='MFP_exp_info_prim_.RData')
# #

# 
# #### GENE CLUSTERS 
# 
# # Allocate Cluster numbers with hclust agglomerative
# 
# # methods to assess
# m <- c( "average", "single", "complete", "ward")
# names(m) <- c( "average", "single", "complete", "ward")
# 
# # function to compute coefficient
# ac <- function(x) {
#   agnes(heatmap_matrix_all, method = x)$ac
# }
# 
# map_dbl(m, ac)
# 
# 
# ### ward method gives highest clustering coefficient (0.935)
# 
# hc3 <- agnes(heatmap_matrix_all, method = "ward")
# pltree(hc3, cex = 0.1, hang = -1, main = "Dendrogram of Genes")
# 
# # Cut tree into 3 groups
# sub_grp <- cutree(hc3, k = 3)
# 
# # Number of members in each cluster
# table(sub_grp)
# 
# #Add gene cluster number 
# heatmap_matrix_df <-heatmap_matrix_all
# heatmap_matrix_df$Hcluster3  <- sub_grp 
# 
# # convert cluster column to data frame
# gene_cluster_df <- as.data.frame(heatmap_matrix_df$Hcluster3)
# 
# #Rename cluster column
# colnames(gene_cluster_df)[colnames(gene_cluster_df)=="heatmap_matrix_df$Hcluster3"] <- "Ward_Cluster"
# #Convert Row Names to Column
# gene_cluster_df$Gene <- rownames(gene_cluster_df)

### Read in PN Gene list

# 
# PN_Genes <- read.csv('PN_Genes.csv')
# PN_Genes <- as.data.frame(PN_Genes)
# 
# #Merge PN Gene Data Frame and   Cluster Data Frame 
# 
# PN_Genes_Cluster_df <- merge(x = PN_Genes,y = gene_cluster_df, by.x ="Gene.Symbol", by.y="Gene"  )
# 
# 
# 
# # Add new column with Ward cluster renumbered so order is same as heatmap cluster
# 
# PN_Genes_Cluster_df$Ward_Cluster_reordered <- sapply(PN_Genes_Cluster_df$Ward_Cluster, function(x)
#   ifelse (x == '3','Primary_1',
#           ifelse     (x == '1','Primary_2', 'Primary_3')
#   ))
# 
# # Rename PN_Genes_Cluster_df$Ward_Cluster_reordered  PN_Genes_Cluster_df$Gene_Cluster
# colnames(PN_Genes_Cluster_df)[colnames(PN_Genes_Cluster_df)=="Ward_Cluster_reordered"] <- "Primary_Gene_Cluster"
# 
# 
# # specify PN_Genes_Cluster_df$Classification as a factor
# PN_Genes_Cluster_df$Classification <- factor(PN_Genes_Cluster_df$Classification,
#                                                levels = c("Autophagy", "Degradation", "Deubiquitination","Ubiquitination","Folding" ))
# 
# 
# 
# 
# #Save PN_Genes_Cluster_df
# 
# save(PN_Genes_Cluster_df, file = 'SKCM_gene_clusters_May.RData')
# 
# 
# # load Gene clusters df  
# PN_Genes_Cluster_df <- loadRData('SKCM_gene_clusters_May.RData')

# # create new column with new names for primary clusters
# PN_Genes_Cluster_df$Primary_Gene_Cluster_new <- sapply(PN_Genes_Cluster_df$Primary_Gene_Cluster, function(x)
#   ifelse (x == 'Primary_1','PN_1(P)',
#           ifelse     (x == 'Primary_2','PN_2(P)', 'PN_3(P)')
#   ))
# 
# #Save PN_Genes_Cluster_df
# 
# save(PN_Genes_Cluster_df, file = 'SKCM_gene_clusters_June.RData')
# 
# # load Gene clusters df  
# PN_Genes_Cluster_df <- loadRData('SKCM_gene_clusters_June.RData')

######### HEATMAPS ######### 



#create matrix for heatmap (gene expression columns transposed so genes are rows and samples are columns )
heatmap_matrix_all <- as.matrix(t(MFP_exp_info_prim_Z))


# ________________________________________
#Create Colour function to distribute colours more evenly ####

heatcolS <- circlize::colorRamp2 (breaks = c(-7,-0.1,0,0.1, 2.5,5, 20),
                                  c("dodgerblue4","#f1fbfe" ,"white", "#ffff55", "#ff4015", "#ff0000", '#5e0406'),
                                  transparency = .3)
heatcolS(seq(-7, 20, 0.01))

# Create top annotation

ann2 <- data.frame (primary_clusters_df$Cluster)


colnames(ann2) <- c('Sample_Cluster_Primary')

ann2$Sample_Cluster_Primary <- factor(ann2$Sample_Cluster_Primary,
                                      levels = c( "SKCM(P)_A","SKCM(P)_B"))



colours2 <- list('Sample_Cluster_Primary' = c('SKCM(P)_A' = 'darkorchid4', 'SKCM(P)_B' = 'gold')
                 
)


colAnn2 <- HeatmapAnnotation(df = ann2,
                             which = 'col',
                             col = colours2,
                             annotation_width = unit(c(2, 10), 'cm'),
                             show_annotation_name = FALSE,
                             annotation_legend_param = list(
                               Sample_Cluster_Primary = list(
                                 title = "Sample Cluster",
                                 #  title_position = "lefttop-rot",
                                 annotation_legend_side="right",
                                 annotation_legend_height = unit(160, "mm"),  
                                 annotation_legend_width = unit(60, "mm"), 
                                 grid_height = unit(2, "cm"), grid_width = unit(0.5, "cm"),
                                 labels_gp = gpar(fontsize = 20),
                                 title_gp = gpar(fontsize = 20, fontface = 'bold'))
                             )
)


# rowAnn <- data.frame (PN_Genes_Cluster_df$Gene.Symbol ,PN_Genes_Cluster_df$Primary_Gene_Cluster_new)
# rowAnn   <- data.frame(rowAnn, row.names = 1)   
# colnames(rowAnn) <- c('Gene_Cluster_Primary')
# colours3 <- list('Gene_Cluster_Primary' = c('PN_1(P)' = 'slateblue', 'PN_2(P)' = 'pink', 'PN_3(P)' = 'darkcyan'))
# 
# rowAnn2 <- HeatmapAnnotation(df = rowAnn,
#                              which = 'row',
#                              col = colours3,
#                              annotation_width = unit(c(2, 10), 'cm'),
#                              show_annotation_name = FALSE,
#                              annotation_legend_param = list(
#                                Gene_Cluster_Primary = list(
#                               #   title_position = "lefttop-rot",
#                                  title = "Gene Cluster",
#                                annotation_legend_side="right",
#                                annotation_legend_height = unit(160, "mm"),  
#                                annotation_legend_width = unit(60, "mm"),
#                                grid_height = unit(2, "cm"), grid_width = unit(0.5, "cm"),
#                                labels_gp = gpar(fontsize = 20),
#                                title_gp = gpar(fontsize = 20, fontface = 'bold'))
#                              ))
# #
## Heatmap with clustering ward.D2 no names ####

pdf (paste0("MFP_exp_info_prim_heatmap_cancer_SKCM.pdf"),width=17, height=7)
#HM <-
Heatmap (heatmap_matrix_all, width = unit(260, "mm"), height = unit(150, "mm"),
         # name ='',
         rect_gp = gpar(col= FALSE),
         
         #   column_title_gp = gpar(fontsize = 20, fontface = "bold"),
         # column_names_centered = TRUE,
         #column_names_rot = 90,
         clustering_method_rows = "ward.D2",
         # clustering_method_columns = "ward.D2",
         # column_names_gp = gpar(fontsize = 25, fontface = "bold"),
         #row_names_gp = gpar(fontsize = 1, fontface = "bold"),
         row_names_gp = gpar(fontsize = 5, fontface = "bold"),
         column_dend_height = unit(10,'mm'),
         row_dend_width=unit(10,"mm"),
         #column_title = "Expression of Proteostasis Network Genes in Primary Skin Cutaneous Melanoma Samples",
         #column_title = "A",
         show_column_names = FALSE,
         show_row_names = TRUE,
         # row_split = 3, 
         column_split = ann2$Sample_Cluster_Primary,
         row_gap = unit(3, "mm"),
         column_gap = unit(3, "mm"),
         col =  heatcolS,
         top_annotation=colAnn2,
         #  right_annotation = rowAnn2 ,
         cluster_column_slices = FALSE,
         
         
         heatmap_legend_param = list(
           title = "Expression", 
           # title_position = "lefttop-rot",
           color_bar = "continuous",
           legend_height = unit(120, "mm"),  
           legend_width = unit(20, "mm"), 
           # fonts = 20, 
           # fontface = "bold",
           title_gp = gpar(fontsize = 20, fontface = 'bold'),
           labels_gp = gpar(col = "black", fontsize = 20),
           heatmap_legend_side="right")
         
)

dev.off()

