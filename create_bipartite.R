count_matrix = "/home/loukos/Desktop/encode/zeta_log_table_scale.txt"
DFD_path = "/home/loukos/Desktop/encode/DFDs2/DFDs/total/"
Organism = "hsapiens"
bipartite_path = "/home/loukos/Desktop/"

create_bipartite <- function(count_matrix, DFD_path, Organism = "hsapiens", bipartite_path)
{
  ############Create bipartite functional-positional  networks##########
  library(igraph)
  library(gprofiler2)
  
  data2<-read.delim(count_matrix, sep="\t")
  for(i in seq(1,length(data2)-5,1))
  #for(i in seq(3,3,1))
  {
    progress = paste("This is experiment: ",i,sep="")
    print(progress)
    #read DFDs of each sequencing experiment
    filename = paste(DFD_path,i,"_tot.txt",sep="")
    data = read.delim(filename, sep=" ")
    g3 <- make_empty_graph(n = 0, directed = F)
    for(j in seq(1,length(data$start),1))
    {
      genes_index = seq(data$start[j],data$end[j],1)
      genes = data2$gene_name[genes_index[1]:genes_index[length(genes_index)]]
      GO_genes_of_DFD <- gost(query = genes, 
                              organism = Organism, ordered_query = FALSE, 
                              multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                              measure_underrepresentation = FALSE, evcodes = FALSE, 
                              user_threshold = 0.05, correction_method = "g_SCS", 
                              domain_scope = "annotated", custom_bg = NULL, 
                              numeric_ns = "", sources = NULL, as_short_link = FALSE)
      if(length(GO_genes_of_DFD$result$term_name) > 0)
      {
        edges = c()
        for(k in seq(1,length(GO_genes_of_DFD$result$term_name),1))
        {
          edges = c(edges,paste(i,"_",j,sep=""),GO_genes_of_DFD$result$term_name[k])
        }
        gtmp <- graph(edges, directed=F)
        E(gtmp)$weight <- GO_genes_of_DFD$result$p_value
        g3 <- union(g3,gtmp, byname = "auto")
      }
    }
    filename3 = paste(bipartite_path,"bipartite_",i,".txt",sep = "")
    write_graph(g3,filename3, format = "gml")
  }
}

create_bipartite("/home/loukos/Desktop/phd/scripttogit/zeta_log_table_scale_demo.txt", "/home/loukos/Desktop/phd/scripttogit/DFDs/", "hsapiens", "/home/loukos/Desktop/phd/scripttogit/bipartites/")
  