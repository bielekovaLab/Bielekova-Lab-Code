get.vals <-  function(data, chunk) {
  
  correlation <- cor.test(data$`Trial 1`,data$`Trial 2`, method = "spearman")
  
  cor.val <- round(correlation[["estimate"]],4)
  p.val <- correlation[["p.value"]]
  
  time <- ((chunk)*90000)/1000
  
  tbl <- cbind(time, cor.val,p.val)
  
  tbl <- as.data.frame(tbl)
  
  return(tbl)
  
}
