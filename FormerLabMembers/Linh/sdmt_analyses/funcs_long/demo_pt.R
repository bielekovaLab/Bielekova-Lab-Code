demo_pt <- function() { 
  
pt <- all_means %>% filter(pid == "ACW365") 

group1 <- pt %>% filter(grouping == 1) 

group2 <- pt %>% filter(grouping == 2) %>% dplyr::select(grouping,pid,y)
x <- seq(1,63,2)
group2 <- data.frame(group2,x) %>% dplyr::select(grouping,pid,x,y)

group3 <- pt %>% filter(grouping == 3) %>% dplyr::select(grouping,pid,y)
x <- seq(1,61,3)
group3 <- data.frame(group3,x) %>% dplyr::select(grouping,pid,x,y)

group4 <- pt %>% filter(grouping == 4) %>% dplyr::select(grouping,pid,y)
x <- seq(1,61,4)
group4 <- data.frame(group4,x) %>% dplyr::select(grouping,pid,x,y)

fin <- bind_rows(group1,group2,group3,group4)

print(fin)

ggplot(fin,aes(x = x, y = y)) + geom_line(alpha = 0.4) + geom_point(shape = 1) + facet_wrap(~grouping, ncol = 4) + theme_bw() +
  theme(axis.ticks.length = unit(-1.4,"mm"),
        axis.text.x = element_text(size = 8, margin = unit(c(t = 2.5, r = 0, b = 0, l = 0), "mm")),
        axis.text.y = element_text(size = 8, margin = unit(c(t = 0, r = 2.5, b = 0, l = 0), "mm")),
        strip.text = element_blank()) 

}
