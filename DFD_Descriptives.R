
DFD_descriptives <- function(num_of_experiments, chrom_length_path, chromosome_levels, DFD_path, plot_path, count_matrix_file, categories_index, deg_threshold)
{  
  #get average length, total number of sig breakpoints, average number per chromosome
  library(ggplot2)
  chrom_lengths = read.delim(chrom_length_path, sep="\t")
  total_genome_length = sum(chrom_lengths$length)
  avg_dfd_length = vector(mode = "double", length = 107)
  dfd_num = vector(mode = "double", length = 107)
  genome_cover_percentage = vector(mode = "double", length = 107)
  chrom_cover_percentage = vector(mode = "double", length = 23)
  tot_chrom_cover_percentage = matrix(0,nrow = 107, ncol = 23)
  #dfd_per_chrom = vector(mode = "double", length = 106)
  data_chrom_init = table(as.factor(chromosome_levels))
  sm1 = data_chrom_init - 1
  gene_number = vector(mode = "double", length = 107)
  avg_genes_per_dfd = vector(mode = "double", length = 107)
  for(i in seq(2,num_of_experiments+1,1))
  {
    filename = paste(DFD_path,i,"_tot.txt",sep="")
    data = read.delim(filename, sep=" ")
    #change columns first coords then gene start end
    data_tmp = data.frame(data[,c(1,4,5)],stringsAsFactors = FALSE)
    data_tmp$chromosome = paste("chr",data_tmp$chromosome,sep = "")
    colnames(data_tmp) <- c("chrom","chromStart","chromEnd")
    dfd_num[i] = length(data$chromosome)
    avg_dfd_length[i] = mean(data$end_coord - data$start_coord)
    tmp1 = 0
    for(a in seq(1,length(data$start_coord),1))
    {
      tmp1 = tmp1 + ((data$end_coord[a] - data$start_coord[a]))
    }
    genome_cover_percentage[i] = tmp1 / total_genome_length
    sm1 = sm1 + summary(as.factor(data$chromosome))
    g = 1
    data$chromosome = as.factor(data$chromosome)
    for(j in levels(data$chromosome))
    {
      data2 = data[data$chromosome == j,]
      chrom_lengths2 = chrom_lengths[chrom_lengths$chromosome == paste(j," ", sep = ""),]
      tmp = 0
      for(b in seq(1,length(data2$start_coord),1))
      {
        tmp = tmp + ((data2$end_coord[b] - data2$start_coord[b]))
      }
      tot_chrom_cover_percentage[i,g] = tmp / chrom_lengths2$length[1]
      g = g + 1
    }
    gene_number[i] = 0
    for(k in seq(1,length(data$end),1))
    {
      gene_number[i] = gene_number[i] + ((data$end[k] - data$start[k]))
    }
  }
  dfd_per_chrom = sm1/106
  dfd_num = dfd_num[dfd_num!=0]
  gene_number = gene_number[gene_number!=0]
  avg_genes_per_dfd = gene_number / dfd_num
  tot_chrom_cover_percentage = tot_chrom_cover_percentage[-1,]
  
  summary(dfd_num)
  filename_hist1 = paste0(plot_path,"DFD_per_experiment_histogram.png")
  png(filename_hist1, width = 2000, height = 1500)
  hist(dfd_num,breaks = 20, main = "DFDs per experiment distribution", xlab = "Number of DFDs")
  dev.off()
  filename_box1 = paste0(plot_path,"DFD_per_experiment_boxplot.png")
  png(filename_box1, width = 1000, height = 1500)
  boxplot(dfd_num, main = "DFD numbers per experiment", ylab = "Number of DFDs")
  dev.off()
  filename_box2 = paste0(plot_path,"DFD_per_chromosome_boxplot.png")
  png(filename_box2, width = 1500, height = 1000)
  boxplot(sm1/106~levels(data$chromosome),cex.axis=1.2, cex.main = 1.2,ylab="Average DFD number", xlab = "Chromosome", main = "Average DFD Numbers per Chromosome")
  dev.off()
  filename_box3 = paste0(plot_path,"DFD_length_boxplot.png")
  png(filename_box3, width = 1500, height = 1000)
  avg_dfd_length = avg_dfd_length[avg_dfd_length!=0]
  boxplot(avg_dfd_length, main = "Average DFD length", ylab = "DFD length in bp")
  dev.off()
  genome_cover_percentage = genome_cover_percentage[genome_cover_percentage!=0]
  filename_hist2 = paste0(plot_path,"Genome_Cover _Percentage_histogram.png")
  png(filename_hist2, width = 2000, height = 1500)
  hist(genome_cover_percentage, breaks = 20, main = "Genome Coverage Percentage", xlab = "DFD length Vs overall genome length")
  dev.off()
  filename_box4 = paste0(plot_path,"Genome_Cover _Percentage_boxplot.png")
  png(filename_box4, width = 1500, height = 1000)
  boxplot(genome_cover_percentage, ylab = "Proportion of Genome in DFDs", main = "Genome Coverage Percentage")
  dev.off()
  
  
  data2<-read.delim(count_matrix_file, sep="\t")
  gene_number <- vector(mode="integer", length = 107)
  deg_number <- vector(mode="integer", length = 107)
  dfd_number <- vector(mode="integer", length = 107)
  for(i in seq(2,107,1))
  {
    progress = paste("This is experiment: ",i-1," / ",num_of_experiments,sep="")
    print(progress)
    #read DFDs of each sequencing experiment
    filename = paste(DFD_path,i,"_tot.txt",sep="")
    data = read.delim(filename, sep=" ")
    GO_total = c()
    gene_list <- list()
    DEGS_list <- list()
    for(j in seq(1,length(data$start),1))
    {
      #isolate DEGS
      #data2[,i] = data2[,i][abs(data2[,i] ) >= 1.5]
      dfd_number[i] = length(data$start)
      DEGs_index = intersect(seq(data$start[j],data$end[j],1),which(abs(data2[,i]) >= deg_threshold))
      genes_index = seq(data$start[j],data$end[j],1)
      if(length(DEGs_index) > 0)
      {
        DEGs = data2$gene_name[DEGs_index[1]:DEGs_index[length(DEGs_index)]] 
      }
      genes = data2$gene_name[genes_index[1]:genes_index[length(genes_index)]]
      gene_list = append(gene_list,genes)
      DEGS_list = append(DEGS_list,DEGs)
    }
    gene_list = unlist(gene_list)
    DEGS_list = unlist(DEGS_list)
    gene_number[i] = length(gene_list)
    deg_number[i] = length(DEGS_list)
  }
  dfd_number = dfd_number[-1]
  deg_number = deg_number[-1]
  gene_number = gene_number[-1]
  
  summary(dfd_number)
  summary(gene_number)
  summary(deg_number)
  ###VALUES TO BE RETURNED
  dfd_num = dfd_number
  deg_num = deg_number
  gene_num = gene_number
  
  ########DO Odds Ratio#########
  odds_ratio = (deg_number/gene_number) / (genome_cover_percentage)
  data_categories <-read.delim(categories_index, sep="\t", header = F) 
  type_of_tissue = data_categories
  or_cell_type = as.data.frame(cbind(type_of_tissue, odds_ratio))
  colnames(or_cell_type) = c("type_of_tissue","odds_ratio")
  krk = kruskal.test(odds_ratio ~ type_of_tissue, or_cell_type)
  krk$p.value
  ###VALUES TO BE RETURNED
  odds_rat = odds_ratio
  odds_rat_pval = krk$p.value
  
  
  #########DO boxplot odds Ratio#########
  filename_box5 = paste0(plot_path,"Odds_Ratio_boxplot.png")
  png(filename_box5, width = 1500, height = 1000)
  boxplot(odds_ratio ~ type_of_tissue, data = or_cell_type, 
          notch = FALSE, # Add notch for median comparison
          varwidth = TRUE, # Adjust box width based on sample size
          xlab = "Category", 
          ylab = "Odds Ratio", 
          main = "Odds Ratio Distribution by Category") 
  mtext(paste("Kruskal-Wallis p-value:", odds_rat_pval), xpos = 0.2, ypos = max(or_cell_type$odds_ratio) + 0.5, adj = 1)
  dev.off()

  #######Boxplot genome cover percentage per category####
  genome_cover_percentage_df = as.data.frame(cbind(type_of_tissue, genome_cover_percentage))
  colnames(genome_cover_percentage_df) = c("type_of_tissue", "genome_cover_percentage")
  filename_box6 = paste0(plot_path,"Genome_Cover_percentage_per_category_boxplot.png")
  png(filename_box6, width = 1500, height = 1000)
  boxplot(genome_cover_percentage ~ type_of_tissue, data = genome_cover_percentage_df, 
          notch = FALSE, # Add notch for median comparison
          varwidth = TRUE, # Adjust box width based on sample size
          xlab = "Category", 
          ylab = "Genome Cover Percentage", 
          main = "Genome Cover Percentage Distribution by Category") 
  dev.off()
  #######Boxplot chromosome cover percentage per category####
  tot_chrom_cover_percentage_df = as.data.frame(cbind(type_of_tissue, tot_chrom_cover_percentage))
  colnames(tot_chrom_cover_percentage_df) = c("type_of_tissue", "tot_chrom_cover_percentage")
  pvalue <- vector(mode="double", length = dim(tot_chrom_cover_percentage)[2])
  for(h in seq(2,length(tot_chrom_cover_percentage_df),1))
  {
    krk = kruskal.test(tot_chrom_cover_percentage_df[,h] ~ tot_chrom_cover_percentage_df[,1])
    pvalue[h-1] = krk$p.value 
  }
  filename_box7 = paste0(plot_path,"Chromosome_Cover_percentage_boxplot.png")
  png(filename_box7, width = 1500, height = 1000)
  boxplot(tot_chrom_cover_percentage,
          notch = FALSE, # Add notch for median comparison
          names = levels(as.factor(chromosome_levels)),
          varwidth = TRUE, # Adjust box width based on sample size
          xlab = "Chromosome", 
          ylab = "Chromosome Cover Percentage", 
          main = "Chromosome Cover Percentage Distribution ")
  dev.off()
  m = 1
  for(h in levels(as.factor(type_of_tissue$V1)))
  {
    tot_chrom_cover_percentage_df_tmp = tot_chrom_cover_percentage_df[tot_chrom_cover_percentage_df[,1] == h,]
    tot_chrom_cover_percentage_tmp = tot_chrom_cover_percentage_df_tmp[,-1]
    filename_boxtmp = paste0(plot_path,"Chromosome_Cover_percentage_boxplot_c", m, ".png")
    png(filename_boxtmp, width = 1500, height = 1000)
    boxplot(tot_chrom_cover_percentage_tmp,
            notch = FALSE, # Add notch for median comparison
            names = levels(as.factor(chromosome_levels)),
            varwidth = TRUE, # Adjust box width based on sample size
            xlab = "Chromosome", 
            ylab = "Chromosome Cover Percentage", 
            main = paste0("Chromosome Cover Percentage Distribution ",h))
    dev.off()
    m = m + 1
  }
  
   
  output <- list(dfd_num = dfd_num, gene_num = gene_num, deg_num = deg_num, cover_percentage = genome_cover_percentage, tot_chrom_cover_percentage = tot_chrom_cover_percentage, tot_chrom_cover_percentage_pval = pvalue, odds_rat = odds_rat, odds_rat_pval = odds_rat_pval, tot_chrom_cover_percentage = tot_chrom_cover_percentage)
  return(output)
}


descriptives <- DFD_descriptives(106, "/home/loukos/Desktop/encode/ncbi_g38_human_chromosome_lengths.txt", c("1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X"), "/home/loukos/Desktop/encode/DFDs2/DFDs/total/", "/home/loukos/Desktop/", "/home/loukos/Desktop/encode/zeta_log_table_scale.txt", "/home/loukos/Desktop/phd/major_cell_type.txt", 1.5)
