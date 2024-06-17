plot_bipartite <- function(count_matrix, bipartite_path, degree_trim = 5, plot_path, plot_width = 8000, plot_height = 8000)
{
  #########Trim and Plot all bipartite networks#############
  library(igraph)
  data<-read.delim(count_matrix, sep="\t")
  data2 = data[is.na(data$gene_name) == 0,] 
  for(i in seq(1,length(data2)-5,1))
  {
    progress = paste("This is experiment: ",i,sep="")
    print(progress)
    #read bipartite graphs of each sequencing experiment
    filename = paste(bipartite_path,"bipartite_",i,".txt",sep="")
    g3 = read_graph(filename, format = "gml")
    # Minimum degree threshold
    min_degree <- degree_trim
    # Filter vertices with degree above threshold
    filtered_vertices <- which(degree(g3) >= min_degree)
    # Trim the graph based on filtered vertices
    trimmed_graph <- induced_subgraph(g3, filtered_vertices)
    ##due to lack of enrichments, some DFD vertexes have a weight of NaN. Replace with 1 to allow for community detection algorithms###
    if(length(which(E(trimmed_graph)$weight == "NaN"))>0)
    {
      E(trimmed_graph)$weight[which(E(trimmed_graph)$weight == "NaN")] <- 1
    }
    #layout <- layout.kamada.kawai(trimmed_graph)
    layout <- layout.fruchterman.reingold
    #layout <- layout_with_fr(trimmed_graph)
    filename2 = paste(plot_path,"bipartite_",i,".png",sep = "")
    main_title = paste0("Plot of Bipartite Network for Experiment ",i)
    png(filename2, width = plot_width, height = plot_height)  # Replace with your desired filename and dimensions
    #plot(trimmed_graph, layout = layout, vertex.label.cex=4)
    k_louvain<-cluster_louvain(trimmed_graph)
    plot(k_louvain, trimmed_graph, vertex.size=V(trimmed_graph)$degree, vertex.label.color="black", vertex.label.family="sans", vertex.label.cex=4, vertex.label.dist=0.5, edge.color="black", main = main_title)
    title(main_title, cex.main = 8)
    dev.off()
  }

}
  
plot_bipartite("/home/loukos/Desktop/phd/scripttogit/zeta_log_table_scale.txt", "/home/loukos/Desktop/encode/DFDs2/DFDs/bipartite/", 8, "/home/loukos/Desktop/",8000, 8000)
plot_bipartite("/home/loukos/Desktop/phd/scripttogit/zeta_log_table_scale_demo.txt", "/home/loukos/Desktop/phd/scripttogit/bipartites/", 8, "/home/loukos/Desktop/",8000, 8000)