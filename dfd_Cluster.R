dfd_Cluster <- function(genome_cover_percentage,avg_dfd_length,avg_genes_per_dfd, tot_chrom_cover_percentage, categories_index, plot_path, tsne = 0, dfds_overlap = 0)
{
  data_categories <-read.delim(categories_index, sep="\t", header = F) 
  type_of_tissue = data_categories
  
  if(tsne == 0)
  {
    #RUN XGBOOST NO tsne
    ##XGBoost
    #reference https://xgboost.readthedocs.io/en/latest/R-package/xgboostPresentation.html
    #reference https://stackoverflow.com/questions/46736934/plotting-the-auc-from-an-xgboost-model-in-r
    #Validate
    library("tidyverse")
    library(caret)
    library("xgboost")
    library("caret")
    library(pROC)
    library("ROCit")
    data_all_model <- data.frame(cbind(genome_cover_percentage,avg_dfd_length,avg_genes_per_dfd,type_of_tissue))
    colnames(data_all_model)[length(data_all_model)] = "type_of_tissue"
    #K-FOLD HOMEMADE
    set.seed(3456)
    k1 = 10
    folds <- createFolds(data_all_model$type_of_tissue, k = k1, list = TRUE, returnTrain = FALSE)
    RMSE1 <- vector(mode = "double", length = 10)
    R2_1 <- vector(mode = "double", length = 10)
    MAE_1 <- vector(mode = "double", length = 10)
    AUC1 <- vector(mode = "double", length = 10)
    for(i in seq(1,k1,1))
    {
      kanTrain <- data_all_model[ folds[[i]],]
      kanTest  <- data_all_model[-folds[[i]],]
      train <- kanTrain
      test <- kanTest
      kanTrainmat <- data.matrix(kanTrain, rownames.force = NA)
      kanTestmat <- data.matrix(kanTest, rownames.force = NA)
      bstTree <- xgboost(data = kanTrainmat, label = kanTrainmat[,length(kanTrainmat[1,])], eta = 0.4, min_split_loss = 0.1, max_depth = 6, min_child_weight = 2, subsample = 0.5, sampling_method = "uniform", nrounds = 10, tree_method = "hist")
      pred <- predict(bstTree, kanTestmat)
      prediction <- round(pred)
      RMSE1[i] = RMSE(prediction, kanTestmat[,length(kanTestmat[1,])])/mean(kanTestmat[,length(kanTestmat[1,])])
      R2_1[i] = R2(prediction, kanTestmat[,length(kanTestmat[1,])])
      MAE_1[i] = MAE(prediction, kanTestmat[,length(kanTestmat[1,])])
      g <- multiclass.roc(type_of_tissue ~ prediction, data = kanTest, levels = as.factor(kanTest$type_of_tissue))
      AUC1[i] = g$auc
    }
    RMSEt = mean(RMSE1)
    R2t = mean(R2_1)
    MAEt = mean(MAE_1)
    AUCt = mean(AUC1)
  
    ###END OF XGBOOST
    #output <- list(RMSE = RMSEt, R2 = R2t, MAE = MAEt, AUC = AUCt)
    #return(output)
  }
  
  if(tsne == 1)
  {
    #RUN XGBOOST WITH tsne
    ##XGBoost
    #reference https://xgboost.readthedocs.io/en/latest/R-package/xgboostPresentation.html
    #reference https://stackoverflow.com/questions/46736934/plotting-the-auc-from-an-xgboost-model-in-r
    #Validate
    library("tidyverse")
    library(caret)
    library("xgboost")
    library("caret")
    library(pROC)
    library("ROCit")
    library(tsne)
    ####Collapse matrix with tsne
    data_overlap <-read.delim(dfds_overlap, sep=" ")
    tsne_overlap = tsne(data_overlap, initial_config = NULL, k = 3,initial_dims = dim(data_overlap[1]), perplexity = 35)
    data_all_model <- data.frame(cbind(genome_cover_percentage,avg_dfd_length,avg_genes_per_dfd, tsne_overlap,type_of_tissue))
    colnames(data_all_model)[length(data_all_model)] = "type_of_tissue"
    #K-FOLD HOMEMADE
    set.seed(3456)
    k1 = 10
    folds <- createFolds(data_all_model$type_of_tissue, k = k1, list = TRUE, returnTrain = FALSE)
    RMSE1 <- vector(mode = "double", length = 10)
    R2_1 <- vector(mode = "double", length = 10)
    MAE_1 <- vector(mode = "double", length = 10)
    AUC1 <- vector(mode = "double", length = 10)
    for(i in seq(1,k1,1))
    {
      kanTrain <- data_all_model[ folds[[i]],]
      kanTest  <- data_all_model[-folds[[i]],]
      train <- kanTrain
      test <- kanTest
      kanTrainmat <- data.matrix(kanTrain, rownames.force = NA)
      kanTestmat <- data.matrix(kanTest, rownames.force = NA)
      bstTree <- xgboost(data = kanTrainmat, label = kanTrainmat[,length(kanTrainmat[1,])], eta = 0.4, min_split_loss = 0.1, max_depth = 6, min_child_weight = 2, subsample = 0.5, sampling_method = "uniform", nrounds = 10, tree_method = "hist")
      pred <- predict(bstTree, kanTestmat)
      prediction <- round(pred)
      RMSE1[i] = RMSE(prediction, kanTestmat[,length(kanTestmat[1,])])/mean(kanTestmat[,length(kanTestmat[1,])])
      R2_1[i] = R2(prediction, kanTestmat[,length(kanTestmat[1,])])
      MAE_1[i] = MAE(prediction, kanTestmat[,length(kanTestmat[1,])])
      g <- multiclass.roc(type_of_tissue ~ prediction, data = kanTest, levels = as.factor(kanTest$type_of_tissue))
      AUC1[i] = g$auc
    }
    RMSEt = mean(RMSE1)
    R2t = mean(R2_1)
    MAEt = mean(MAE_1)
    AUCt = mean(AUC1)
    
    ###END OF XGBOOST
    #output <- list(RMSE = RMSEt, R2 = R2t,MAE = MAEt, AUC = AUCt)
    #return(output)
  }
  library(viridis)
  library("dendextend")
  library("stats")
  library(gplots)
  
  dchrom_cover_percentage <- dist(genome_cover_percentage)
  chrom_cover_percentage_clust <- hclust(dchrom_cover_percentage,method="ward.D2")
  cl <- cutree(chrom_cover_percentage_clust,k=length(levels(as.factor(type_of_tissue$V1))))
  categories_cl <- data.frame("Categories" = type_of_tissue$V1, "Cluster" = factor(cl))
  cluster_genome_cover = table(categories_cl)
  filename_clusterbar = paste0(plot_path,"Genome_Cover_percentage_clustering_barplot.png")
  png(filename_clusterbar, width = 1500, height = 1000)
  par(cex = 1.5)
  plot(table(categories_cl),color = viridis(length(levels(as.factor(type_of_tissue$V1)))), main = "Clusters by Genome Cover Percentage")
  dev.off()

  dchrom_cover_percentage <- dist(tot_chrom_cover_percentage)
  chrom_cover_percentage_clust <- hclust(dchrom_cover_percentage,method="ward.D2")
  cl <- cutree(chrom_cover_percentage_clust,k=length(levels(as.factor(type_of_tissue$V1))))
  categories_cl <- data.frame("Categories" = type_of_tissue$V1, "Cluster" = factor(cl))
  cluster_chromosome_cover = table(categories_cl)
  filename_clusterbar = paste0(plot_path,"Chromosome_Cover_percentage_clustering_barplot.png")
  png(filename_clusterbar, width = 1500, height = 1000)
  par(cex = 1.5)
  plot(table(categories_cl),color = viridis(length(levels(as.factor(type_of_tissue$V1)))), main = "Clusters by Chromosome Cover Percentage")
  dev.off()
  clcols2 <- c(viridis(length(levels(as.factor(categories_cl$Cluster)))))[categories_cl$Cluster]
  hcols <- c(magma(length(levels(as.factor(categories_cl$Categories)))))[categories_cl$Cluster]
  filename_clusterheat = paste0(plot_path,"Chromosome_Cover_percentage_clustering_heatmap.png")
  png(filename_clusterheat, width = 1500, height = 1000)
  par(cex = 1.5)
  heatmap.2(tot_chrom_cover_percentage, hclustfun =
              function(x) hclust(x,method = 'ward.D2'), 
            scale="none", mar=c(5, 10), trace="none",
            RowSideColors = clcols2, labRow = type_of_tissue$V1, main = "Clusters by Chromosome Cover Percentage")
  dev.off()
  
  output <- list(model = bstTree, RMSE = RMSEt, R2 = R2t,MAE = MAEt, AUC = AUCt, cluster_genome_cover = cluster_genome_cover, cluster_chromosome_cover = cluster_chromosome_cover)
  return(output)
}
