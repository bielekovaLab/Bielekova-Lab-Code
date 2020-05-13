plot_cpt <- function(print_out) {
  regress_vals <- data.frame(c("AAO02", "ABX699", "ABX711",
                               "ACW072", "ACW204", "ACW254",
                               "ACW269", "ACW282", "ACW365",
                               "ACW406", "ACW407", "ACW484", "ACW570","ABX710"),
                             c(9,3,3,3,3,4,3,3,3,4,3,3,10,14),
                             rep(60,14), c(2,2,5,5,5,5,5,5,5,4,4,4,3,9),
                             c("TRUE", "FALSE", "FALSE", "FALSE",
                               "TRUE", "FALSE", "FALSE", "FALSE",
                               "TRUE", "FALSE", "FALSE", "FALSE",
                               "TRUE","TRUE"),
                             c("FALSE", "FALSE", "FALSE", "FALSE",
                               "FALSE", "FALSE", "FALSE", "FALSE",
                               "FALSE", "FALSE", "FALSE", "FALSE",
                               "TRUE", "TRUE"))
  
  colnames(regress_vals ) <- c("PID", "c", "b0", "b1", "axis1", "axis2")
  
  cpt_data <- bind_rows(lapply(1:nrow(regress_vals), function(x) all_cpt(regress_vals$PID[x], regress_vals)))
  
  addtl <- transformed %>% filter(PID == "ACW245" | PID == "ACW559") %>%
           dplyr::select(PID, x,y) %>% group_by(PID) %>% nest() %>% 
           mutate(model = data %>% map(~lm(y ~ x, data = .))) %>% 
           mutate(pred_y = map2(model, data, predict)) %>% 
           unnest(pred_y, data) %>% mutate(cpt = rep(NA, n())) %>%
           dplyr::select(PID, x, y, pred_y, cpt)
  
  out <- cpt_data %>% group_by(PID) %>% distinct(cpt, .keep_all = TRUE)
  
  print(out)
  
  cpt_data <- bind_rows(cpt_data, addtl)
  
  cpt_pre <- cpt_data %>% filter(
      PID == "AAO02" & x <= 6 |
      PID == "ABX699" & x <= 3 |
      PID == "ABX711" & x <= 10 |
      PID == "ACW072" & x <= 4  |
      PID == "ACW204" & x <= 12 |
      PID == "ACW254" & x <= 6 |
      PID == "ACW269" & x <= 3 |
      PID == "ACW282" & x <= 5  |
      PID == "ACW365" & x <= 5  |
      PID == "ACW406" & x <= 13|
      PID == "ACW407" & x <= 8 |
      PID == "ACW484" & x <= 4 |
      PID == "ACW570" & x <= 12|
      PID == "ABX710" & x <= 15 |
      PID == "ACW245" | PID == "ACW559")
  
  cpt_post <- cpt_data %>% filter(
    PID == "AAO02" & x > 6 |
      PID == "ABX699" & x > 3 |
      PID == "ABX711" & x > 10 |
      PID == "ACW072" & x > 4  |
      PID == "ACW204" & x > 12 |
      PID == "ACW254" & x > 6 |
      PID == "ACW269" & x > 3 |
      PID == "ACW282" & x > 5  |
      PID == "ACW365" & x > 5  |
      PID == "ACW406" & x > 13|
      PID == "ACW407" & x > 8 |
      PID == "ACW484" & x > 4 |
      PID == "ACW570" & x > 12|
      PID == "ABX710" & x > 15)
  
  # data for putting in regression line to the data after change point
  after_cpt <- all.long.rm %>% group_by(PID, Date, Time) %>%
               ungroup() %>% group_by(PID, Date) %>%
               summarise(score = mean(Correct)*90) %>%
               ungroup() %>% group_by(PID) %>% mutate(sitting = rank(Date)) %>%
               filter(PID == "AAO02" & sitting > 6 |
             PID == "ABX699" & sitting > 3 |
             PID == "ABX711" & sitting > 10 |
             PID == "ACW072" & sitting > 4  |
             PID == "ACW204" & sitting > 12 |
             PID == "ACW254" & sitting > 6 |
             PID == "ACW269" & sitting > 3 |
             PID == "ACW282" & sitting > 5  |
             PID == "ACW365" & sitting > 5  |
             PID == "ACW406" & sitting > 13 |
             PID == "ACW407" & sitting > 8 |
             PID == "ACW484" & sitting > 4 |
             PID == "ACW570" & sitting > 12|
             PID == "ABX710" & sitting > 15)
  
  if(print_out == TRUE) { 
  p <- ggplot(data = cpt_data, aes(x = x, y = y)) + geom_line() +
    geom_line(data = cpt_pre, aes(x = x, y = pred_y), color = "chocolate3", size = 0.5) + 
    geom_line(data = cpt_post, aes(x = x, y = pred_y), color = "chocolate3", size = 0.5,
              linetype = "dashed", size = 2)+
    stat_smooth(data = after_cpt, aes(x = sitting, y = score), method = "lm", 
                color = "blue", size = 0.5, geom='line', alpha = 0.5, se = FALSE) + ylab("Smartphone SDMT Score")
  
  q <- p + 
    facet_wrap(~PID) + theme_bw() + 
    geom_vline(aes(xintercept = cpt), linetype = "dotted") + 
    theme(axis.ticks.length=unit(-0.20, "cm"), axis.text.x = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), 
          axis.text.y = element_text(margin=unit(c(0.5,0.5,0.5,0.5), "cm")), strip.text.x = element_blank())
  
  print(q)
  }
  
  else if(print_out == FALSE) {
    
    corvals <- after_cpt %>% group_by(PID) %>% mutate(rho = (cor.test(sitting,score, method = "spearman")$estimate),
                                                      p.vals = cor.test(sitting,score, method = "spearman")$p.value) %>%
      filter(Date == min(Date))
    
  return(corvals)
    
  }
  
  else(NULL)
}

