quant_bin <- function(x,quants=c(0,.5,1),...){
  quantiles<-quantile(x,quants)
  quantiles[1]<-quantiles[1]-0.001
  cut(x,breaks=quantiles,...)
}
