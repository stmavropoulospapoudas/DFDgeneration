
overlap_convert <- function(dfds_overlap_input, dfds_overlap_output)
{
  command = paste0("sed '/^File/!d' ", dfds_overlap_input, " > ", dfds_overlap_output)
  system(command, intern = FALSE)
}
  
