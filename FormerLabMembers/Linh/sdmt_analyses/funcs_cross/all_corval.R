all_corval <- function(split.data, cohort) { 
  
max.time <- 90000

chunk.interval <- 5000

chunk.list <- seq(5000, max.time, by=chunk.interval)

chunk.list <- chunk.list/90000

cor.val <- bind_rows(lapply(1:length(chunk.list), 
                               function(x) vals.tbl(split.data,chunk.list[x], cohort)))

return(cor.val)

}
