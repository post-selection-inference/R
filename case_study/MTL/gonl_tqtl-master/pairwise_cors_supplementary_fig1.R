########################################################

# This script will produce pairwise correlation plots
# for the metadata available for the filtered 248 families 
# this script was used to produce Supplementary figure 1 

########################################################



#########################################
###         Load data                 ###
#########################################

load("child.df.RData") # phenotype data for children and parents 248 families


###########################################
###       Pairwise correlations         ###
###########################################
# select numeric columns
child.df = child.df[,c("MTL", "mMTL", "fMTL", "Age", "MAC", "PAC", "mAge", "fAge")]

# scale
child.df = as.data.frame(scale(child.df))

# plot pairwise correlations and residuals



par(mfrow = c(ncol(child.df), ncol(child.df)))
par(mar = c(2, 1, 2, 1))


for(j.ind in 1:ncol(child.df)){
  j = colnames(child.df)[j.ind]
  for(k.ind in 1:ncol(child.df)){
    k = colnames(child.df)[k.ind]
    if( k == j){
      plot(1,1, type = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")
      text(x = 1,y = 1, labels = j, col = "black", cex = 2)
    }
    else if (k.ind > j.ind){
      fit = lm(child.df[,j] ~ child.df[,k])
      sf = summary(fit)
      plot(child.df[,k], child.df[,j], pch = 19, cex.main = 1.3,
           main = paste0("R.sqr ", round(sf$r.squared,2),
                         "  b1 ", round(sf$coefficients[2,1],2)))
      if(sf$coefficients[2,4] > 0.05)
        color = "gray"
      else if (sf$coefficients[2,4] < 0.001)
        color = "red"
      else if (sf$coefficients[2,4] < 0.01)
        color = "orange"
      else 
        color = "green"
      abline(a = sf$coefficients[1,1], b= sf$coefficients[2,1], col = color, lwd = 2)
    } else {
      fit = lm(child.df[,k] ~ child.df[,j])
      sf = summary(fit)
      plot(child.df[,j],sf$residuals, pch = 19, col = "gray")
      if(sf$coefficients[2,4] > 0.05)
        color = "gray"
      else if (sf$coefficients[2,4] < 0.001)
        color = "red"
      else if (sf$coefficients[2,4] < 0.01)
        color = "orange"
      else 
        color = "green"
      abline(h = 0, col = color)
      
    }
  }
}

# plot legend
par(mfrow = c(1,1))
plot.new()
legend(x = 0, y = 1, 
       legend = c("p > 0.05", "p < 0.05", "p < 0.01", "p < 0.001"), 
       col = c("gray", "green", "orange", "red"), 
       lwd = 3, bty = "n", cex = 2)

