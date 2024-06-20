deg_enrichments_DFDs <- function(DEG_DFD_Enrichment_path, categories_index)
{
  ########Do DEG Enrichment in DFDs##########
  DEG_DFD_Enrichment <-read.delim(DEG_DFD_Enrichment_path, header = F, sep=" ")
  data_categories <-read.delim(categories_index, sep="\t", header = F) 
  type_of_tissue = data_categories
  ##try with all overlap
  DEG_DFD_Enrichment_df = as.data.frame(cbind(type_of_tissue$V1, DEG_DFD_Enrichment))
  colnames(DEG_DFD_Enrichment_df)[1] <- "type_of_tissue"
  krk = kruskal.test(V10 ~ type_of_tissue, DEG_DFD_Enrichment_df)
  a = krk$p.value
  ##try only with significant overlap
  DEG_DFD_Enrichment_df = DEG_DFD_Enrichment_df[DEG_DFD_Enrichment_df$V12 < 0.05,]
  krk = kruskal.test(V10 ~ type_of_tissue, DEG_DFD_Enrichment_df)
  b = krk$p.value
  output <- list(all_enrichments_pval = a, significant_enrichments_pval = b)
  return(output) 
}
