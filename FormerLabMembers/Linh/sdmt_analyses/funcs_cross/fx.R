# Functional form for NLS, NLME
fx <- function(x,c,b0,b1){
  out <- b0 + as.numeric(x < c)*b1*x + as.numeric(x>=c)*b1*c
  return(out)
}