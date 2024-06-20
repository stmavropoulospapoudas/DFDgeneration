create_DFDs <- function(count_matrix_file, breakpoint_path, DFD_path, chromosome_levels, significance_threshold = 0.1)
{
  ###Get dataset in format #geneid, #experiment 1, #experiment 2, ... , #experiment n, #genename, #start, #end, #chromosome
  ##### the #experiment X columns should have differential expression levels from the genes
  #If differential expression data is not available, define DEGs as more than 1.5 SD far from mean (z-score >= 1.5, line 23
  data<-read.delim(count_matrix_file, sep="\t")
  #create dataset duplicates
  data5 = data
  data3 = data
  data2 = data[is.na(data$gene_name) == 0,] 
  library("TTR")
  library("strucchange")
  n = 1
  
  ##start from columns 2 (exclude #geneid column) to total columns - 4 (exclude #genename, #start, #end, #chromosome columns)
  for(i in seq(2,length(data2)-4,1))
  {
   progress = paste0("Now doing experiment: ",n)
   print(progress)
   #isolate all relevant values (chromosome, coordinates, gene names, expression values) for each experiment
   data3 = data2[,c("chromosome","start","end","gene_name")]   
   data3[,1] = data2$chromosome
   data3[,2] = data2$start
   data3[,3] = data2$end
   data3[,4] = data2$gene_name
   data3[,5] = data2[,i]
   colnames(data3)<-c("chrom", "start", "end", "gene","exp_val")
   #isolate "DEGs" based on expression threshold
   #define DEGs as more than 1.5 SD far from mean
   degs = data3[abs(data3$exp_val) >= 1.5,]
   #set genes with undefined expression values to 0
   data3[is.na(data3)] <- 0 
   #create DFDs
   #set sliding window for smooth moving average
   win=50
   #circle through al chromosomes 
   #seq1 = c(1:22,"X")
   #seq1 = levels(as.factor(data3$chromosome))
   for(chrom in chromosome_levels)
   { 
     #chrom = 1
     print(chrom)
     #select values per chromosome
     data_chrom = data3[data3$chrom == chrom,]
     #sort genes by start coordinates
     data_chrom = data_chrom[order(data_chrom$start),]
     #convert expression values to timeseries
     exp_1 = ts(data_chrom$exp_val)
     #get smooth moving average
     exp_1_SMA <- SMA(exp_1,n=win)
     exp_1_SMA[is.na(exp_1_SMA)] <- 0 ###discuss convert NAs to 0s
     #get breakpoints per chromosome 
     ##ERROR HANDLER TO AVOID COLLAPSE
     tryCatch({
       breakpoints(exp_1_SMA ~ 1, h=0.05,  dynamic = FALSE, rescale = FALSE, type="fluctuation")->bp;confint(bp)->ci;
     }, error=function(e){print("error");breakpoints(exp_1_SMA ~ 1, h=0.1,  dynamic = FALSE, rescale = FALSE, type="fluctuation")->bp;confint(bp)->ci;
  })
     breakpoints <- data.frame(chromosome=character(length(bp$breakpoints)),
                                           start=integer(length(bp$breakpoints)),
                                           end=integer(length(bp$breakpoints)),
                                           start_coords=integer(length(bp$breakpoints)),
                                           end_coords=integer(length(bp$breakpoints)),
                                           score=double(length(bp$breakpoints)),
                                           stringsAsFactors=FALSE)
     breakpoints$chromosome[1] = chrom
     breakpoints$start[1] = 1
     breakpoints$end[1] = bp$breakpoints[1]
     breakpoints$start_coords[1] = data_chrom$start[1]
     breakpoints$end_coords[1] = data_chrom$end[bp$breakpoints[1]]
     breakpoints$score[1] = mean(exp_1_SMA[1:bp$breakpoints[1]])
     g = 2
     for(p in seq(1,length(bp$breakpoints)-1,1))
     {
       breakpoints$chromosome[g] = chrom
       breakpoints$start[g] = bp$breakpoints[p]
       breakpoints$end[g] = bp$breakpoints[p+1]
       breakpoints$start_coords[g] = data_chrom$start[bp$breakpoints[p]]
       breakpoints$end_coords[g] = data_chrom$end[bp$breakpoints[p+1]]
       breakpoints$score[g] = mean(exp_1_SMA[bp$breakpoints[p]:bp$breakpoints[p+1]])
       g = g + 1
     }
     breakpoints$chromosome[length(bp$breakpoints)] = chrom
     breakpoints$start[p+1] = bp$breakpoints[p+1]
     breakpoints$end[p+1] = length(exp_1_SMA)
     breakpoints$start_coords[p+1] = data_chrom$start[bp$breakpoints[p+1]]
     breakpoints$end_coords[p+1] = data_chrom$end[length(exp_1_SMA)]
     breakpoints$score[p+1] = mean(exp_1_SMA[bp$breakpoints[p+1]:length(exp_1_SMA)])
     ####produce breakpoint bed files
     filename = paste(breakpoint_path,i,"chromosome",chrom,".bed",sep="")
     write.table(breakpoints, file=filename, quote = FALSE, row.names = F)
     
     #breakpoints(exp_1_SMA ~ 1, h=0.05,  dynamic = FALSE, rescale = FALSE, type="fluctuation")->bp;confint(bp)->ci;plot(exp_1_SMA);lines(ci)
     #create dataframe to fill with DFDs
     significant_breakpoints <- data.frame(chromosome=character(length(bp$breakpoints)),
                                           start=integer(length(bp$breakpoints)),
                                           end=integer(length(bp$breakpoints)),
                                           start_coord=integer(length(bp$breakpoints)),
                                           end_coord=integer(length(bp$breakpoints)),
                                           value=double(length(bp$breakpoints)),
                                           stringsAsFactors=FALSE)
     #set significance threshold
     threshold = significance_threshold
     #see if first breakpoint is significant 
     if(mean(abs(exp_1_SMA[1:bp$breakpoints[1]])) >= threshold)
     {
       significant_breakpoints$chromosome[1] = chrom
       significant_breakpoints$start[1] = 1
       significant_breakpoints$end[1] = bp$breakpoints[1]
       significant_breakpoints$start_coord[1] = data_chrom$start[1]
       significant_breakpoints$end_coord[1] = data_chrom$end[bp$breakpoints[1]]
       significant_breakpoints$value[1] = mean(exp_1_SMA[1:bp$breakpoints[1]])
     }
     #see if the rest of breakpoints are significant
     g = 2
     for(k in seq(1,length(bp$breakpoints)-1,1))
     {
       if(mean(abs(exp_1_SMA[bp$breakpoints[k]:bp$breakpoints[k+1]])) >= threshold)
       {
         significant_breakpoints$chromosome[g] = chrom
         significant_breakpoints$start[g] = bp$breakpoints[k]
         significant_breakpoints$end[g] = bp$breakpoints[k+1]
         significant_breakpoints$start_coord[g] = data_chrom$start[bp$breakpoints[k]]
         significant_breakpoints$end_coord[g] = data_chrom$end[bp$breakpoints[k+1]]
         significant_breakpoints$value[g] = mean(exp_1_SMA[bp$breakpoints[k]:bp$breakpoints[k+1]])
         g = g + 1
       }
     }
     #see if last breakpoint is significant
     if(mean(abs(exp_1_SMA[bp$breakpoints[k+1]:length(exp_1_SMA)])) >= threshold)
     {
       significant_breakpoints$chromosome[k+1] = chrom
       significant_breakpoints$start[k+1] = bp$breakpoints[k+1]
       significant_breakpoints$end[k+1] = length(exp_1_SMA)
       significant_breakpoints$start_coord[k+1] = data_chrom$start[bp$breakpoints[k+1]]
       significant_breakpoints$end_coord[k+1] = data_chrom$end[length(exp_1_SMA)]
       significant_breakpoints$value[k+1] = mean(exp_1_SMA[bp$breakpoints[k+1]:length(exp_1_SMA)])
     }
     filename = paste(DFD_path,i,"chromosome",chrom,".bed",sep="")
     significant_breakpoints = significant_breakpoints[significant_breakpoints$start != 0,]
     write.table(significant_breakpoints, file=filename, quote = FALSE, row.names = F)
   }
   DFDs_total = paste(DFD_path,i,"_tot.txt",sep="")
   command <- paste0("cd ", DFD_path, " && awk 'FNR>1' *.bed > ",DFDs_total)
   system(command, intern = FALSE)
   command <- paste0("rm ", DFD_path, "*.bed")
   system(command, intern = FALSE)
   breakpoints_total = paste(breakpoint_path,"breakpoints",i,".txt",sep="")
   command <- paste0("cd ", breakpoint_path, " && awk 'FNR>1' *.bed > ",breakpoints_total)
   system(command, intern = FALSE)
   command <- paste0("rm ", breakpoint_path, "*.bed")
   system(command, intern = FALSE)
   n = n + 1
  }
}
