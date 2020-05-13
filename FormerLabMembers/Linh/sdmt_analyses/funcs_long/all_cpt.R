all_cpt <- function(patientcode, vals_table) {
  
  c <- subset(vals_table, PID == patientcode)$c
  b0 <- subset(vals_table, PID == patientcode)$b0
  b1 <- subset(vals_table, PID == patientcode)$b1
  
  find_cpt(transformed, patientcode, c, b0, b1)
  
}
