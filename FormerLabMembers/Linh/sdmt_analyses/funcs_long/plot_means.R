plot_means <- function(){
  
  p <- ggplot(means_diff, aes(x=grouping, y=diff, group=grouping)) + 
    geom_violin(trim=FALSE, fill="gray87") + geom_boxplot(width=0.1) + theme_minimal() +
         theme(axis.ticks.length = unit(-1.4,"mm"),
                           axis.text.x = element_text(size = 8, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")),
                           axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm"))) +
        
        # adding 1.5*IQR +/- quantile lines to the plot. grouping below: 1
         annotate('segment',x = 0, y = quantile(subset(means_diff, grouping == 1)$diff)[2] - 1.5*IQR(subset(means_diff, grouping == 1)$diff),
                                                 xend = 1, yend = quantile(subset(means_diff, grouping == 1)$diff)[2] - 1.5*IQR(subset(means_diff, grouping == 1)$diff), 
                                        colour = "darkgreen", alpha = 0.3) +
         annotate('segment',x = 0, y = quantile(subset(means_diff, grouping == 1)$diff)[4] + 1.5*IQR(subset(means_diff, grouping == 1)$diff), 
                                                xend = 1, yend = quantile(subset(means_diff, grouping == 1)$diff)[4] + 1.5*IQR(subset(means_diff, grouping == 1)$diff),
                                        colour = "darkgreen", alpha = 0.3) + 
        
         # adding 1.5*IQR +/- quantile lines to the plot. grouping below: 2
       annotate('segment', x = 0, y = quantile(subset(means_diff, grouping == 2)$diff)[2] - 1.5*IQR(subset(means_diff, grouping == 2)$diff),
                                                xend = 2, yend = quantile(subset(means_diff, grouping == 2)$diff)[2] - 1.5*IQR(subset(means_diff, grouping == 2)$diff), 
                                       colour = "blue", alpha = 0.3) +
       annotate('segment', x = 0, y = quantile(subset(means_diff, grouping == 2)$diff)[4] + 1.5*IQR(subset(means_diff, grouping == 2)$diff), 
                                                xend = 2, yend = quantile(subset(means_diff, grouping == 2)$diff)[4] + 1.5*IQR(subset(means_diff, grouping == 2)$diff), 
                                       colour = "blue", alpha = 0.3) + 
         
        # adding 1.5*IQR +/- quantile lines to the plot. grouping below: 3
         annotate('segment',x = 0, y = quantile(subset(means_diff, grouping == 3)$diff)[2] - 1.5*IQR(subset(means_diff, grouping == 3)$diff),
                                                 xend = 3, yend = quantile(subset(means_diff, grouping == 3)$diff)[2] - 1.5*IQR(subset(means_diff, grouping == 3)$diff), 
                                         colour = "purple", alpha = 0.3) +
         annotate('segment',x = 0, y = quantile(subset(means_diff, grouping == 3)$diff)[4] + 1.5*IQR(subset(means_diff, grouping == 3)$diff), 
                                               xend = 3, yend = quantile(subset(means_diff, grouping == 3)$diff)[4] + 1.5*IQR(subset(means_diff, grouping == 3)$diff), 
                                        colour = "purple", alpha = 0.3) +
         
        # adding 1.5*IQR +/- quantile lines to the plot. grouping below: 4
        annotate('segment',x = 0, y = quantile(subset(means_diff, grouping == 4)$diff)[2] - 1.5*IQR(subset(means_diff, grouping == 4)$diff),
                                               xend = 4, yend = quantile(subset(means_diff, grouping == 4)$diff)[2] - 1.5*IQR(subset(means_diff, grouping == 4)$diff),
                                      colour = "black", alpha = 0.3) +
        annotate('segment', x = 0, y = quantile(subset(means_diff, grouping == 4)$diff)[4] + 1.5*IQR(subset(means_diff, grouping == 4)$diff), 
                                                xend = 4, yend = quantile(subset(means_diff, grouping == 4)$diff)[4] + 1.5*IQR(subset(means_diff, grouping == 4)$diff),
                                       colour = "black", alpha = 0.3) + scale_x_continuous(expand = c(0,0))
  
  print(p)
  
  print(paste("difference range for grouping of 1 sitting:", 
              quantile(subset(means_diff, grouping == 1)$diff)[2] - 1.5*IQR(subset(means_diff, grouping == 1)$diff), "-",
              quantile(subset(means_diff, grouping == 1)$diff)[4] + 1.5*IQR(subset(means_diff, grouping == 1)$diff)))
  
  print(paste("difference range for grouping of 2 sittings:", 
              quantile(subset(means_diff, grouping == 2)$diff)[2] - 1.5*IQR(subset(means_diff, grouping == 2)$diff), "-",
              quantile(subset(means_diff, grouping == 2)$diff)[4] + 1.5*IQR(subset(means_diff, grouping == 2)$diff)))
  
  print(paste("difference range for grouping of 3 sittings:", 
              quantile(subset(means_diff, grouping == 3)$diff)[2] - 1.5*IQR(subset(means_diff, grouping == 3)$diff), "-",
              quantile(subset(means_diff, grouping == 3)$diff)[4] + 1.5*IQR(subset(means_diff, grouping == 3)$diff)))
  
  print(paste("difference range for grouping of 4 sittings:", 
              quantile(subset(means_diff, grouping == 4)$diff)[2] - 1.5*IQR(subset(means_diff, grouping == 4)$diff), "-",
              quantile(subset(means_diff, grouping == 4)$diff)[4] + 1.5*IQR(subset(means_diff, grouping == 4)$diff)))
}

