jaccard_bootstrapping <- function(jaccard_similarity_matrix, categories_index, num_permutations = 10000)
{
  data_categories <-read.delim(categories_index, sep="\t", header = F) 
  type_of_tissue = data_categories
  type_of_tissue = type_of_tissue[type_of_tissue!=""]
  #bind cell type to dataframe/matrix
  df_jaccard = as.data.frame(jaccard_similarity_matrix)
  df_jaccard = cbind(df_jaccard,type_of_tissue)
  #colnames(df_jaccard) = df_jaccard$type_of_tissue
  categories = levels(as.factor(type_of_tissue))
  pvalue_vector <- vector(mode = "double", length = length(categories))
  expected_similarity <- vector(mode = "double", length = length(categories))
  observed_similarity <- vector(mode = "double", length = length(categories))
  m = 1
  for(k in categories)
  {
    print(k)
    df_jaccard_tmp = df_jaccard[(df_jaccard$type_of_tissue == k),]
    df_jaccard_tmp = rbind(df_jaccard_tmp, type_of_tissue)
    #remove last row wich will be text character category values
    df_jaccard_tmp = df_jaccard_tmp[, !(names(df_jaccard_tmp) %in% "type_of_tissue")]
    last_row = df_jaccard_tmp[nrow(df_jaccard_tmp),]
    select_columns_index = which(last_row == k)
    colnamesjac = colnames(df_jaccard_tmp)
    df_jaccard_tmp = df_jaccard_tmp[,select_columns_index]
    colnames(df_jaccard_tmp) = colnamesjac[select_columns_index]
    df_jaccard_tmp = df_jaccard_tmp[-nrow(df_jaccard_tmp),]
    num_jaccard = dim(df_jaccard_tmp)[1] * dim(df_jaccard_tmp)[2]
    median_jaccard = median(as.double(as.matrix(df_jaccard_tmp)))
    simulated_jaccard_matrix <- vector(mode = "double",length = num_jaccard)
    simulated_jaccard_median <- vector(mode = "double",length = num_permutations)
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
    pvalue=length(which(simulated_jaccard_median>median_jaccard))/num_permutations
    if(pvalue > 0.95)
    {
      pvalue=-length(which(simulated_jaccard_median<median_jaccard))/num_permutations
    }
    pvalue_vector[m] = pvalue
    expected_similarity[m] = mean(simulated_jaccard_median)
    observed_similarity[m] = median_jaccard
    m = m + 1
  }
  pvalue_vector = as.data.frame(rbind(categories, expected_similarity, observed_similarity, pvalue_vector))
  return(pvalue_vector)
}
