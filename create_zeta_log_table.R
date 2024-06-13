
create_zeta_log_table <- function(raw_count_file, gene_names_file, count_table_path)
{
  #get gene names
  gene_names = read.delim(gene_names_file,sep = "\t") 
  dataset_1 = read.delim(raw_count_file,sep = "\t") 
  num_of_experiments = dim(dataset_1)[2]
  
  #initialise vectors
  gene_name <- vector(mode = "character", length = length(dataset_1$gene_id))
  start <- vector(mode = "character", length = length(dataset_1$gene_id))
  end <- vector(mode = "character", length = length(dataset_1$gene_id))
  chromosome <- vector(mode = "character", length = length(dataset_1$gene_id))
  
  #get gene name and coordinates for each gene
  for(i in seq(1,length(dataset_1$gene_id),1))
  {
    print(paste0(((i/length(dataset_1$gene_id))) * 100 , " %"))
    tst1 = strsplit(toString(dataset_1$gene_id[i]), "")[[1]]
    a = ""
    if(tst1[1] == "E")
    {  for(j in seq(1,length(tst1),1))
    {
      if(tst1[j] != ".")
      {
        a = c(a,tst1[j])
      }
      if(tst1[j] == ".")
      {
        break
      }
    }
      a = a[a != ""]
      a = paste(a,collapse = "")
    }
    gene_names_tmp = gene_names[gene_names$Gene.stable.ID == a,]
    gene_name[i] = toString(gene_names_tmp$Gene.name[1])
    start[i] = toString(gene_names_tmp$Gene.start..bp.[1])
    end[i] = toString(gene_names_tmp$Gene.end..bp.[1])
    chromosome[i] = toString(gene_names_tmp$Chromosome.scaffold.name[1])
  }
  
  #add new info to original dataframe
  dataset_1 = cbind(dataset_1, gene_name, start, end, chromosome)
  dataset_1[dataset_1 == 0] <- NA
  
  #convert values to log
  dataset_1[,2:num_of_experiments] = log10(dataset_1[,2:num_of_experiments])
  
  dataset_2 = dataset_1
  dataset_2[,2:num_of_experiments] = scale(dataset_1[,2:num_of_experiments])
  # dataset_3 = dataset_1
  # for(i in seq(2,num_of_experiments,1))
  # {
  #   dataset_3[,i] = scale(dataset_1[,i])
  # }
  filename = paste0(count_table_path,"zeta_log_table_scale.txt")
  write.table(dataset_2, filename, quote = FALSE, row.names = F, sep = "\t")
}

create_zeta_log_table("/home/loukos/Desktop/encode/Supplemental_table_S2.csv", "/home/loukos/Desktop/encode/encode_genes.txt", "/home/loukos/Desktop/")
data<-read.delim("/home/loukos/Desktop/zeta_log_table_scale.txt", sep="\t", header = TRUE)



