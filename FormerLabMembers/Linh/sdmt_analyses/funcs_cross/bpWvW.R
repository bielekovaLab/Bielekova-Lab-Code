bpWvW <- function(data) {
  
  matched <- fig4B
  
  print(head(matched))
  print(wilcox.test(Paper_Correct~Diagnosis, data = matched))
  
  print(boxplot(Paper_Correct~Diagnosis, data = matched,
                ylab = "Written SDMT Score"))
  print(median(subset(matched, Diagnosis == "HV")$Paper_Correct))
  print(median(subset(matched, Diagnosis == "MS")$Paper_Correct))
  
}
  