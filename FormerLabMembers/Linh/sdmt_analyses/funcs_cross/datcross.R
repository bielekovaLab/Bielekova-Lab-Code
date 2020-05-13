# function to create a working dataset for cross-sectional data
# with outliers removed, if necessary
datcross <- function(input) {
 
  
  output <- input %>% group_by(PID) %>% 
            filter(Date == min(Date)) %>% ungroup()
 
  return(output)

  
}
