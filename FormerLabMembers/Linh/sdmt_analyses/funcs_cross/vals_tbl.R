vals.tbl <- function(data,chunk,cohort) {
  
  chunked.data <- find_vals(data,chunk,cohort)
  
  tbl <- get.vals(chunked.data,chunk)
  
  return(tbl)
}