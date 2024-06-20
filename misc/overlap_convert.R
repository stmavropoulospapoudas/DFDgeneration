overlap_convert <- function(dfds_overlap_input, dfds_overlap_output, create_overlap_matrix = 0)
{
  command = paste0("sed '/^File/!d' ", dfds_overlap_input, " > ", dfds_overlap_output)
  system(command, intern = FALSE)
  if(create_overlap_matrix == 1)
  {
    data_overlap = read.delim(dfds_overlap_output, sep=" ", header = F)
    data_overlap_tmp <- split(data_overlap$V10, cut(seq_along(data_overlap$V10), sqrt(length(data_overlap$V10)), labels=FALSE))
    overlap_matrix = c()
    for(i in seq(1,length(sqrt(data_overlap$V10)),1))
    {
      #print(as.double(unlist(data_overlap_tmp[i])))
      overlap_matrix = rbind(overlap_matrix, as.double(unlist(data_overlap_tmp[i])))
    }
    return(overlap_matrix) 
  }
}
