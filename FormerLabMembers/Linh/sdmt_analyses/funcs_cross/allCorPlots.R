allCorPlots <- function(split.data) { 

cor.val.HV <- all_corval(split.data, "HV")
subset.HV <- subset(cor.val.HV, time >= 75)

cor.val.MS <- all_corval(split.data, "MS")
subset.MS <- subset(cor.val.MS, time >= 75)

y <- cor.val.HV$cor.val
x <- cor.val.HV$time

mod.hv <- nls(y~fx(x,c,b0,b1), start=c("c"=60,"b0"= 60,"b1"=min(cor.val.MS$time)),
           control=nls.control(warnOnly = TRUE,minFactor = 1e-20,maxiter=1000))
stable.hv <- lm(cor.val~time, data = subset.HV)

out.hv <- summary(mod.hv)

y <- cor.val.MS$cor.val
x <- cor.val.MS$time

mod.ms <- nls(y~fx(x,c,b0,b1), start=c("c"=60,"b0"= 60,"b1"=min(cor.val.MS$time)),
              control=nls.control(warnOnly = TRUE,minFactor = 1e-20,maxiter=1000))
stable.ms <- lm(cor.val~time, data = subset.MS)

out.ms <- summary(mod.ms)

yrange <- range(cor.val.MS$cor.val, cor.val.HV$cor.val)

par(mfrow = c(1,2))

plot(cor.val~time, data = cor.val.HV, pch = 16, col = "black",
     xlab = "Test Progression Time (seconds)", ylab = "Spearman Correlation Coefficient",
     main = "HV Trial 1 v. 2 Correlation Over Time", ylim = yrange, tck = 0.02)

abline(v = round(out.hv$coefficients[1]), col = "black", lty = 3, lwd = 2)
lines(x, predict(mod.hv), col = "chocolate3", lwd = 2)
lines(subset.HV$time, predict(stable.hv), col = "deepskyblue4", lwd = 2)
abline(v = 75, lty = 2, lwd = 2)
text(out.hv$coefficients[1]+6, min(cor.val.HV$cor.val), paste0(round(out.hv$coefficients[1]),"s"), font = 2)

plot(cor.val~time, data = cor.val.MS, pch = 16, col = "black",
     xlab = "Test Progression Time (seconds)", ylab = "Spearman correlation Coefficient",
     main = "MS Trial 1 v. 2 Correlation Over Time", ylim = yrange, tck = 0.02)

abline(v = round(out.ms$coefficients[1]), col = "black", lty = 3, lwd = 2)
lines(x, predict(mod.ms), col = "chocolate3", lwd = 2)
lines(subset.MS$time, predict(stable.ms), col = "deepskyblue4",  lwd = 2)
abline(v = 75, lty = 2, lwd = 2)
text(out.ms$coefficients[1]+5, min(cor.val.HV$cor.val), paste0(round(out.ms$coefficients[1]),"s"), font = 2)

colnames(cor.val.HV) <- c("Time", "Cor.HV", "P.HV")
colnames(cor.val.MS) <- c("Time", "Cor.MS", "P.MS")

return(cbind(cor.val.HV, cor.val.MS))

}
