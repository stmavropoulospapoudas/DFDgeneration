

get_jaccard_matrix <- function(count_matrix, GO_terms_path)
{
  ####function for jaccard similarity
  jaccard_similarity_index <- function(a, b) {
    intersection = length(intersect(a, b))
    union = length(a) + length(b) - intersection
    return (intersection/union)
  }
  
  data2<-read.delim(count_matrix, sep="\t")
  ############get grid all Vs all similarity
  jaccard_similarity_matrix <- matrix(0,nrow = length(data2)-5, ncol = length(data2)-5)
  for(i in seq(1,length(data2)-5,1))
  {
    progress = paste("This is experiment: ",i,"/",length(data2)-5,sep="")
    print(progress)
    #read DFDs of each sequencing experiment
    filename = paste(GO_terms_path,"GO_terms_",i,".txt",sep="")
    data_go = read.delim(filename)
    for(j in seq(1,length(data2)-5,1))
    {
      if(i != j)
      {
        filename = paste(GO_terms_path,"GO_terms_",j,".txt",sep="")
        data_go_tmp = read.delim(filename)
        jaccard_similarity_matrix[i,j] = jaccard_similarity_index(data_go[,1],data_go_tmp[,1])
      }
    }
  }  
  return(jaccard_similarity_matrix)
}  
  

