convert_to_bed <-function(DFD_path, bed_path, number_of_files)
{
  for(i in seq(1,number_of_files,1))
  {
    filename = paste(DFD_path,i,"_tot.txt",sep="")
    data = read.delim(filename, sep=" ")
    #change columns first coords then gene start end
    data_tmp = data.frame(data[,c(1,4,5)],stringsAsFactors = FALSE)
    data_tmp$chromosome = paste("chr",data_tmp$chromosome,sep = "")
    colnames(data_tmp) <- c("chrom","chromStart","chromEnd")
    filenamew = paste0(bed_path,"DFDs_",i,".bed")
    write.table(data_tmp,filenamew,col.names = F, quote = F)
  }
}

