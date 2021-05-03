
##### create covariance matrix for random error terms #####

create_sigma2 <- function(datlist, R2){
  
  rowlist <- sapply(datlist, nrow)
  
  # no observed value for all measures
  if (all(rowlist == 0)){
    
    Sigma2 <- NULL
    
  }else{
    
    # if there is no observed value for any one of the measures, remove empty dataframe
    if (any(rowlist == 0)){
      
      datlist <- datlist[rowlist != 0]
      
    }
    
    Sigma2 <- NULL
    
    for (i in 1:length(datlist)){
      
      Cmat <- NULL
      
      for (j in 1:length(datlist)){
        
        time1 <- datlist[[i]][, 2]; time2 <- datlist[[j]][, 2]
        
        if (i == j){
          
          # diagonal matrices
          cmat <- R2[i, j] * diag(1, length(time1))
          
        }else{
          
          # off diagonal matrices
          covmat <- outer(1:length(time1), 1:length(time2), FUN = "paste", sep = ",")
          
          comb <- cbind(1:length(time1), match(time1, time2, nomatch = NA))
          
          
          if(all(is.na(comb[, 2]))){
            
            cmat <- matrix(0, nrow = length(time1), ncol = length(time2))
            
          }else{
            
            covmat[covmat %in% apply(comb, 1, paste, collapse = "," )] <- 1
            
            covmat[covmat != 1] <- 0
            
            cmat <- (covmat %>% as.data.frame %>%
                       sapply(function(x){as.numeric(levels(as.factor(x))[as.factor(x)])})) * R2[i, j] 
          }
        }
        
        if(is.null(dim(cmat))){
          
          Cmat <- matrix(c(Cmat, cmat), nrow = 1)
          
        }else{
          
          Cmat <- cbind(Cmat, cmat)
        }
        
      }
      
      Sigma2 <- rbind(Sigma2, Cmat)
      
    }
    
    
  }
  
  
  return(Sigma2)
  
}
