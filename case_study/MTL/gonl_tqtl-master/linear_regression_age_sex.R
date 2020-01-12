########################################################

# This script will produce statistics and figures for 
# linear regression of MTL on age and sex
# of all individuals across the filtered 248 families 

########################################################





#########################################
###         Load data                 ###
#########################################

load("total.df.RData") # phenotype data of all individuals



###########################################
###   Linear regression on age and sex  ###
###########################################

fit = lm(total.df$MTL/1000 ~ total.df$Age + total.df$Sex)
sf = summary(fit)
sex.color = c("M"="deepskyblue4", "F" = "lightpink3")

par(mar = c(5.1, 5.1, 4.1, 2))
plot(total.df$Age, total.df$MTL/1000, cex.axis = 2,
     pch = 19, col = sex.color[total.df$Sex], 
     xlab = "Age, years", ylab = "MTL, kbp", cex.lab = 2, 
     main = "MTL is negatively correlated with age", cex.main = 2)
abline(fit$coefficients[1:2], col = "gray17", lwd = 2) 

legend(x = 65, y = 12, cex = 1.5,
       legend = c(paste0("Adjusted R2 = ", round(sf$adj.r.squared, 2)), 
                  paste0("p value < 2.2e-16")),
       bty = "n")

legend(x = 20, y = 4, col = sex.color, legend = c("Male", "Female"), pch = 16, cex = 1)

###########################################
###           MTL sex difference        ###
###########################################

sex.split = split(total.df$MTL/1000, total.df$Sex)

boxplot(sex.split, 
        names = c("Female", "Male"), cex.axes = 2, 
        ylab = "MTL, kbp", cex.lab = 2,
        main = "MTL difference between sexes", cex.main =2)

stripchart(sex.split, method = "jitter", add = T, vertical = T, 
           pch = 19, col = c("indianred4", "dodgerblue4"))

legend(x = 1.1, y = 8.4, bty = "n", cex = 1.2,
       legend = c(
         paste0(round(median(sex.split$F), 1), " +/- ", 
                round(sd(sex.split$F), 1), " kbp")))

legend(x = 2.1, y = 8.2, bty = "n", cex = 1.2,
       legend = c(
         paste0(round(median(sex.split$M), 1), " +/- ", 
                round(sd(sex.split$M), 1), " kbp")))

wt = wilcox.test(sex.split$F, sex.split$M)
wt$p.value

# adjustment for age

age.split = split(total.df, total.df$Age)

diffs = unlist(lapply(age.split, function(x) {
  ss = split(x, x$Sex)
  if(length(ss$F) > 0 & length(ss$M > 0))
    diff = mean(ss$F$MTL) - mean(ss$M$MTL)
  else 
    diff = NA
  return(diff)
}  ))

na.ind = which(is.na(diffs))
diffs = diffs[-na.ind]
t.sexdiff.age = t.test(diffs)
age.split = age.split[-na.ind]

# age-wise plot of sex differences in MTL 

plot(as.numeric(names(age.split)), diffs, pch = 16, 
     xlab = "age", ylab = "female MTL - male MTL", 
     main = "Sex difference in MTL per age group")

abline(h = 1)

legend(x = "topright", legend = paste0("p value = ", round(t.sexdiff.age$p.value, 2)))
