bipartite_cluster <- function(count_matrix, bipartite_path, categories_index, plot_path, xgboost = 0)
{
  library(ggplot2)
  library(igraph)
  library(stats)
  library(viridis)
  library("xgboost")
  library("caret")
  library(pROC)
  library("ROCit")
  library("tidyverse")
  library(caret)
  data2<-read.delim(count_matrix, sep="\t")
  ######Cycle through each graph to get modularity and out degree (undirected graph in = out = total/2) and connected components etc
  modularity <- vector(mode = "double", length = length(data2)-5)
  degree_vector <- vector(mode = "double", length = length(data2)-5)
  components_vector <- vector(mode = "double", length = length(data2)-5)
  comp_size <- vector(mode = "double", length = length(data2)-5)
  comp_size_tot <- vector(mode = "double", length = length(data2)-5)
  for(i in seq(1,length(data2)-5,1))
  {
    #read bipartite graphs of each sequencing experiment
    filename = paste(bipartite_path,"bipartite_",i,".txt",sep="")
    g4 = read_graph(filename, format = "gml")
    wtc <- cluster_walktrap(g4)
    #modularity(wtc)
    modularity[i] = modularity(g4, membership(wtc))
    #print(modularity(g4, membership(wtc)))
    degree_vector[i] = mean(degree(g4,V(g4),mode = "all"))
    components_vector[i] = count_components(g4, mode = "weak")
    b = components(g4, "weak")
    comp_size[i] = mean(b$csize)
    comp_size_tot[i] = sum(b$csize)
  }
  # modularity = modularity[-1]  
  # degree_vector = degree_vector[-1]
  # components_vector = components_vector[-1]
  # comp_size = comp_size[-1]
  # comp_size_tot = comp_size_tot[-1]
  
  
  ###Compare modularity between cell types
  data_categories <-read.delim(categories_index, sep="\t", header = F) 
  type_of_tissue = data_categories
  type_of_tissue = type_of_tissue[type_of_tissue!=""]
  
  df_bipartites_modularity = as.data.frame(cbind(modularity,type_of_tissue))
  colnames(df_bipartites_modularity) <- c("values","idx")
  df_bipartites_modularity$values = as.double(df_bipartites_modularity$values)
  mkrk = kruskal.test(values~idx,data=df_bipartites_modularity)
  mkrk$p.value
  
  df_bipartites_degree = as.data.frame(cbind(degree_vector,type_of_tissue))
  colnames(df_bipartites_degree) <- c("values","idx")
  df_bipartites_degree$values = as.double(df_bipartites_degree$values)
  dkrk = kruskal.test(values~idx,data=df_bipartites_degree)
  dkrk$p.value
  
  df_bipartites_components = as.data.frame(cbind(components_vector,type_of_tissue))
  colnames(df_bipartites_components) <- c("values","idx")
  df_bipartites_components$values = as.double(df_bipartites_components$values)
  krk = kruskal.test(values~idx,data=df_bipartites_components)
  krk$p.value
  
  df_bipartites_components_size = as.data.frame(cbind(comp_size,type_of_tissue))
  colnames(df_bipartites_components_size) <- c("values","idx")
  df_bipartites_components_size$values = as.double(df_bipartites_components_size$values)
  c1krk = kruskal.test(values~idx,data=df_bipartites_components_size)
  c1krk$p.value
  
  df_bipartites_components_size_tot = as.data.frame(cbind(comp_size_tot,type_of_tissue))
  colnames(df_bipartites_components_size_tot) <- c("values","idx")
  df_bipartites_components_size_tot$values = as.double(df_bipartites_components_size_tot$values)
  c2krk = kruskal.test(values~idx,data=df_bipartites_components_size_tot)
  c2krk$p.value
  
  #########Do boxplot of degree per major cell type####
  # Create the ggplot object
  ggplot(df_bipartites_degree, aes(x = idx, y = values)) +

    # Add geometry for boxplots
    geom_boxplot(aes(fill = idx)) +

    # Set labels and title
    labs(title = "Degree by Category",
         x = "Category",
         y = "Degree") +

    # Rotate x-axis labels for better readability (optional)
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename = "Degree Distributions.png", path = plot_path)

  #########Do boxplot of modularity per major cell type####
  # Create the ggplot object
  ggplot(df_bipartites_modularity, aes(x = idx, y = values)) +

    # Add geometry for boxplots
    geom_boxplot(aes(fill = idx)) +

    # Set labels and title
    labs(title = "Modularity by Category",
         x = "Category",
         y = "Modularity") +

    # Rotate x-axis labels for better readability (optional)
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename = "Modularity Distributions.png", path = plot_path)

  #########Do boxplot of components size per major cell type####
  # Create the ggplot object
  ggplot(df_bipartites_components_size, aes(x = idx, y = values)) +

    # Add geometry for boxplots
    geom_boxplot(aes(fill = idx)) +

    # Set labels and title
    labs(title = "Connected Component Mean Size by Category",
         x = "Category",
         y = "Mean Size of Connected Components") +

    # Rotate x-axis labels for better readability (optional)
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(filename = "Mean Connected Component Size Distributions.png", path = plot_path)

  df_metrics <- as.data.frame(cbind(modularity, degree_vector, components_vector, comp_size, comp_size_tot, type_of_tissue))
  colnames(df_metrics) <- c("modularity", "mean_degree", "number_of_components", "component_mean_size", "total_component_size", "category")
  
  #######try connected components mean size first
  #Calculate pairwise distances using Ward's method
  distance_matrix <- dist(df_metrics$component_mean_size)
  #Perform hierarchical clustering using hclust
  ward.clusters <- hclust(distance_matrix, method = "ward.D")
  #plot(ward.clusters, main = "Ward Dendrogram")
  #Assign cluster labels to each data point (cut the dendrogram)
  # Choose number of clusters (k) 
  k <- length(levels(as.factor(type_of_tissue))) 
  cluster_labels <- cutree(ward.clusters, k = k)
  #Add the cluster labels as a new column to data frame
  df_metrics$comp_mean_size_cluster <- cluster_labels
  # Check cluster frequency
  table(df_metrics$category,df_metrics$comp_mean_size_cluster)
  # Create a frequency table
  freq_table <- table(df_metrics$category,df_metrics$comp_mean_size_cluster)
  # Generate colors based on the number of unique categories
  colors <- viridis(n = nlevels(factor(df_metrics$category)))
  # Create the bar plot with ggplot2
  ggplot(as.data.frame(freq_table)) +
    aes(x = Var1, y = Freq, fill = Var2) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors, breaks = names(levels(factor(df_metrics$category))), labels = names(levels(factor(df_metrics$category)))) +
    labs(title = "Category by Cluster based on Connected Components Mean Size", x = "Category", y = "Frequency", fill = "category") +
    theme(legend.title = element_text(size = 12, face = "bold"),  # Adjust legend title options
          legend.position = "topright",  # Position legend at top right
          legend.text = element_text(size = 10))  # Adjust legend text size
  ggsave(filename = "Mean Connected Component Size Clusters.png", path = plot_path)
  
  ####try all 3 features (modularity, degree, connected components mean size)
  data_matrix = as.matrix(cbind(df_metrics$modularity,df_metrics$mean_degree,df_metrics$component_mean_size))
  # Set the number of clusters (k)
  num_clusters <- length(levels(as.factor(df_metrics$category)))
  # Perform k-means clustering
  kmeans_model <- kmeans(data_matrix, centers = num_clusters, nstart = 20)
  # Get cluster labels for each data point
  cluster_labels <- kmeans_model$cluster
  df_metrics$all_features_cluster <- cluster_labels
  table(df_metrics$category,df_metrics$all_features_cluster)
  freq_table <- table(df_metrics$category,df_metrics$all_features_cluster)
  colors <- viridis(n = nlevels(factor(df_metrics$category)))
  # Create the bar plot for clusters all features
  ggplot(as.data.frame(freq_table)) +
    aes(x = Var1, y = Freq, fill = Var2) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colors, breaks = names(levels(factor(df_metrics$category))), labels = names(levels(factor(df_metrics$category)))) +
    labs(title = "All Features Frequency of Category by Cluster", x = "Category", y = "Frequency", fill = "category") +
    theme(legend.title = element_text(size = 12, face = "bold"),  # Adjust legend title options
          legend.position = "topright",  # Position legend at top right
          legend.text = element_text(size = 10))  # Adjust legend text size
  ggsave(filename = "All Features Clusters.png", path = plot_path)
  
  
  ##########XGBoost for Bipartites######
  if(xgboost == 1)
  {
    Cell_Type = type_of_tissue
    data_all_model <- data.frame(cbind(df_metrics$component_mean_size,Cell_Type))
    #K-FOLD HOMEMADE
    set.seed(3456)
    k1 = 10
    folds <- createFolds(data_all_model$Cell_Type, k = k1, list = TRUE, returnTrain = FALSE)
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
      bstTree <- xgboost(data = kanTrainmat, label = as.integer(as.factor(as.character(kanTrain$Cell_Type))), eta = 0.4, min_split_loss = 0.1, max_depth = 6, min_child_weight = 2, subsample = 0.5, sampling_method = "uniform", nrounds = 10, tree_method = "hist")
      pred <- predict(bstTree, kanTestmat)
      prediction <- round(pred)
      #df1 = data.frame( R2 = R2(prediction, kanTest$Cell_Type),
      #RMSE = RMSE(prediction, kanTest$Cell_Type),
      #MAE = MAE(prediction, kanTest$Cell_Type))
      RMSE1[i] = RMSE(prediction, as.integer(as.factor(kanTest$Cell_Type)))/mean(as.integer(as.factor(kanTest$Cell_Type)))
      R2_1[i] = R2(prediction, as.integer(as.factor(kanTest$Cell_Type)))
      MAE_1[i] = MAE(prediction, as.integer(as.factor(kanTest$Cell_Type)))
      g <- multiclass.roc(Cell_Type ~ prediction, data = kanTest, levels = as.factor(kanTest$Cell_Type))
      AUC1[i] = g$auc
    }
    RMSE2 = mean(RMSE1)
    R22 = mean(R2_1)
    MAE2 = mean(MAE_1)
    AUC2 = mean(AUC1)
    boost_metrics = as.data.frame(cbind(RMSE2, R22, MAE2, AUC2))
    colnames(boost_metrics) <- c("RMSE", "R2", "MAE", "AUC")
  }
  if(xgboost == 2)
  {
    Cell_Type = type_of_tissue
    data_all_model <- data.frame(cbind(df_metrics$modularity,df_metrics$mean_degree,df_metrics$component_mean_size,Cell_Type))
    #K-FOLD HOMEMADE
    set.seed(3456)
    k1 = 10
    folds <- createFolds(data_all_model$Cell_Type, k = k1, list = TRUE, returnTrain = FALSE)
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
      bstTree <- xgboost(data = kanTrainmat, label = as.integer(as.factor(as.character(kanTrain$Cell_Type))), eta = 0.4, min_split_loss = 0.1, max_depth = 6, min_child_weight = 2, subsample = 0.5, sampling_method = "uniform", nrounds = 10, tree_method = "hist")
      pred <- predict(bstTree, kanTestmat)
      prediction <- round(pred)
      #df1 = data.frame( R2 = R2(prediction, kanTest$Cell_Type),
      #RMSE = RMSE(prediction, kanTest$Cell_Type),
      #MAE = MAE(prediction, kanTest$Cell_Type))
      RMSE1[i] = RMSE(prediction, as.integer(as.factor(kanTest$Cell_Type)))/mean(as.integer(as.factor(kanTest$Cell_Type)))
      R2_1[i] = R2(prediction, as.integer(as.factor(kanTest$Cell_Type)))
      MAE_1[i] = MAE(prediction, as.integer(as.factor(kanTest$Cell_Type)))
      g <- multiclass.roc(Cell_Type ~ prediction, data = kanTest, levels = as.factor(kanTest$Cell_Type))
      AUC1[i] = g$auc
    }
    RMSE2 = mean(RMSE1)
    R22 = mean(R2_1)
    MAE2 = mean(MAE_1)
    AUC2 = mean(AUC1)
    boost_metrics = as.data.frame(cbind(RMSE2, R22, MAE2, AUC2))
    colnames(boost_metrics) <- c("RMSE", "R2", "MAE", "AUC")
  } 
  
  ###END OF XGBOOST for bipartites
  if(xgboost != 0)
  {
    output <- list(metrics = df_metrics, modularity_pval = mkrk$p.value, degree_pval = dkrk$p.value, mean_components_size_pval = c1krk$p.value, total_components_size_pval = c2krk$p.value, boost_metrics = boost_metrics, model = bstTree)
  }
  if(xgboost == 0)
  {
    output <- list(metrics = df_metrics, modularity_pval = mkrk$p.value, degree_pval = dkrk$p.value, mean_components_size_pval = c1krk$p.value, total_components_size_pval = c2krk$p.value)
  }
  return(output)
  
}

bmetrics <- bipartite_cluster("/home/loukos/Desktop/phd/scripttogit/zeta_log_table_scale_demo.txt", "/home/loukos/Desktop/phd/scripttogit/bipartites/", "/home/loukos/Desktop/phd/scripttogit/major_cell_type_demo.txt", "/home/loukos/Desktop", 1)
bmetrics <- bipartite_cluster("/home/loukos/Desktop/phd/scripttogit/zeta_log_table_scale_demo.txt", "/home/loukos/Desktop/phd/scripttogit/bipartites/", "/home/loukos/Desktop/phd/scripttogit/major_cell_type_demo.txt", "/home/loukos/Desktop", 2)
bmetrics <- bipartite_cluster("/home/loukos/Desktop/phd/scripttogit/zeta_log_table_scale_demo.txt", "/home/loukos/Desktop/phd/scripttogit/bipartites/", "/home/loukos/Desktop/phd/scripttogit/major_cell_type_demo.txt", "/home/loukos/Desktop")

table(bmetrics$metrics$comp_mean_size_cluster,bmetrics$metrics$category)









