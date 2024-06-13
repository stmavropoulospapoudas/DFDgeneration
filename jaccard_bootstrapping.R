jaccard_bootstrapping <- function(jaccard_similarity_matrix, categories_index, num_permutations = 10000)
{
  data_categories <-read.delim(categories_index, sep="\t", header = F) 
  type_of_tissue = data_categories
  type_of_tissue = type_of_tissue[type_of_tissue!=""]
  #bind cell type to dataframe/matrix
  df_jaccard = as.data.frame(jaccard_similarity_matrix)
  df_jaccard = cbind(df_jaccard,type_of_tissue)
  df_jaccard = rbind(df_jaccard, type_of_tissue)
  #colnames(df_jaccard) = df_jaccard$type_of_tissue
  categories = levels(as.factor(type_of_tissue))
  pvalue_vector <- vector(mode = "double", length = length(categories))
  expected_similarity <- vector(mode = "double", length = length(categories))
  observed_similarity <- vector(mode = "double", length = length(categories))
  m = 1
  for(k in categories)
  {
    df_jaccard_tmp = df_jaccard[(df_jaccard$type_of_tissue == k),]
    df_jaccard_tmp = df_jaccard_tmp[, !colSums(df_jaccard_tmp == k)]
    #remove last line wich will be text character category values
    df_jaccard_tmp = df_jaccard_tmp[-dim(df_jaccard_tmp),]
    df_jaccard_tmp = df_jaccard_tmp[,-dim(df_jaccard_tmp)]
    num_jaccard = dim(df_jaccard_tmp)[1] * dim(df_jaccard_tmp)[2]
    median_jaccard = median(as.double(as.matrix(df_jaccard_tmp)))
    simulated_jaccard_matrix <- vector(mode = "double",length = num_jaccard)
    simulated_jaccard_median <- vector(mode = "double",length = 10000)
    for(i in seq(1,num_permutations,1))
    {
      print("This is Simulation: ")
      print(i)
      simulated_jaccard_matrix <- vector(mode = "double",length = num_jaccard)
      for(j in seq(1,num_jaccard,1))
      {
        indexi = as.integer(runif(1)*length(df_jaccard$V1))
        indexj = as.integer(runif(1)*length(df_jaccard$V1))
        if(indexi != 0 && indexj != 0)
        {
          simulated_jaccard_matrix[j] = df_jaccard[indexi,indexj]
        }
        if(indexi == 0 || indexj == 0)
        {
          j = j - 1
        }
      }
      simulated_jaccard_median[i] = median(as.double(simulated_jaccard_matrix))
    }  
    pvalue=length(which(simulated_jaccard_median>median_jaccard))/10000
    if(pvalue > 0.95)
    {
      pvalue=-length(which(simulated_jaccard_median<median_jaccard))/10000
    }
    pvalue_vector[m] = pvalue
    expected_similarity[m] = mean(simulated_jaccard_median)
    observed_similarity[m] = median_jaccard
    m = m + 1
  }
  pvalue_vector = as.data.frame(rbind(categories, expected_similarity, observed_similarity, pvalue_vector))
  return(pvalue_vector)
}

pvalues = jaccard_bootstrapping(jaccard_similarity_matrix, "/home/loukos/Desktop/phd/major_cell_type.txt", 1000)
pvalues = jaccard_bootstrapping(jaccard_similarity_matrix, "/home/loukos/Desktop/phd/major_cell_type.txt")