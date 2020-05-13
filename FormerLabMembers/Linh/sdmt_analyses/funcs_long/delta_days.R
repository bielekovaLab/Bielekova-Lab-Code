delta.days <- function(input) {
  
dat <- input %>% group_by(PID) %>% arrange(Date) %>% mutate(diff = Date - min(Date))
  
p <- ggplot(data = dat, aes(x = diff, y = y)) + geom_line() + facet_wrap(~PID) + theme_bw() +
  theme(axis.ticks.length=unit(-0.20, "cm"), axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
        axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), strip.text.x = element_blank())

print(p)

}
