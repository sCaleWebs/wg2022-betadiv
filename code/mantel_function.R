library(vegan)


mantel_loo_function <- function(abundance_mat, environment_mat) {
  
  environment_dist <- dist(environment_mat) 
  
  beta_div <- vegdist(abundance_mat, binary=FALSE, method = "bray")
  
  mantel_res <- mantel(beta_div, environment_dist)
  
  mantel_loo_all <- c()
  
  for(i in 1:ncol(abundance_mat)){
    
    new_abundance <- abundance_mat[, -i]
    
    beta_loo <- vegdist(new_abundance, binary=FALSE, method = "bray")
    
    mantel_loo <- mantel(beta_loo, environment_dist)
    
    mantel_loo_all <- c(mantel_loo_all, mantel_loo$statistic)
    
  }
  
  mantel_all_df <- data.frame(bwg_name =   colnames(abundance_mat), mantel_r_loo = mantel_loo_all, mantel_original_r = mantel_res$statistic)
  
  return(mantel_all_df)
}
