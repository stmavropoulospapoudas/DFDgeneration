
bipartite_metrics <- function(count_matrix, bipartite_path, categories_index, plot_path)
{
  library(ggplot2)
  library(igraph)
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
  #modularity = modularity[-1]  
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
  colnames(df_metrics) <- c("modularity", "mean degree", "number_of_components", "component_mean_size", "total_component_size", "category")
  output <- list(metrics = df_metrics, modularity_pval = mkrk$p.value, degree_pval = dkrk$p.value, mean_components_size_pval = c1krk$p.value, total_components_size_pval = c2krk$p.value)
  return(output)

}

