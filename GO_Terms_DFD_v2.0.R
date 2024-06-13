go_terms_DFD <- function(count_matrix, DFD_path, GO_path, Organism)
{
  ###########Get all GO term Enrichments for each DFD in each experiment#################
  print("Now getting GO enrichments...")
  library(gprofiler2)
  data2<-read.delim(count_matrix, sep="\t")
  for(i in seq(2,length(data2)-4,1))
  {
    progress = paste0("Now doing experiment: ",i-1)
    print(progress)
    #read DFDs of each sequencing experiment
    filename = paste(DFD_path,i,"_tot.txt",sep="")
    data = read.delim(filename, sep=" ")
    #create empty vector to gather all terms
    GO_total = c()
    for(j in seq(1,length(data$chromosome),1))
    {
      print("This is iteration")
      print(j)
      #get gene list for genes inside each DFD
      genes_index = seq(data$start[j],data$end[j],1)
      genes_in_DFDs = data2$gene_name[genes_index[1]:genes_index[length(genes_index)]]
      #get GO enrichments for DFD
      if(length(genes_in_DFDs) > 0)
      {
        GO_of_DFD <- gost(query = genes_in_DFDs, 
                          organism = Organism, ordered_query = FALSE, 
                          multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                          measure_underrepresentation = FALSE, evcodes = FALSE, 
                          user_threshold = 0.05, correction_method = "g_SCS", 
                          domain_scope = "annotated", custom_bg = NULL, 
                          numeric_ns = "", sources = NULL, as_short_link = FALSE)
        #get enrichments of all DFDs in experiment
        GO_total = c(GO_total, GO_of_DFD$result$term_name)
      }
      
    }
    filename2 = paste(GO_path,"GO_terms_",i-1,".txt",sep = "")
    write.table(GO_total,filename2, quote = F, col.names = F, row.names = F,sep = "\t")
  }
}

go_terms_DFD("/home/loukos/Desktop/phd/scripttogit/zeta_log_table_scale.txt", "/home/loukos/Desktop/encode/DFDs2/DFDs/total/", "/home/loukos/Desktop/", "hsapiens")
