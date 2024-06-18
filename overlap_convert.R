dfds_overlap = "/home/loukos/Desktop/results_DFDS.txt"
dfds_overlap2 = "/home/loukos/Desktop/results_DFDS2.txt"
command = paste0("sed '/^File/!d' ", dfds_overlap, " > ", dfds_overlap2)
system(command, intern = FALSE)
#filename = paste(dfds_overlap,i,"_tot.txt",sep="")
data_overlap = read.delim(dfds_overlap2, sep=" ", header = F)
data_overlap$V10
#data_overlap$V12
for(i in seq(1,length(data_overlap$V10),1))
{
  
}
overlap_convert <- function(dfds_overlap_input, dfds_overlap_output)
{
  command = paste0("sed '/^File/!d' ", dfds_overlap_input, " > ", dfds_overlap_output)
  system(command, intern = FALSE)
}
  