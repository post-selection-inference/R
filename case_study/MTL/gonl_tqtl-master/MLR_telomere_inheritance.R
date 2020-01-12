########################################################

# This script will produce pairwise correlation plots
# for the metadata available for the filtered 248 families 
# this script was used to produce Supplementary figure 1 

########################################################





#########################################
###         Load data                 ###
#########################################

load("child.df.RData") # phenotype data for children and parents 248 families


###################################################
###     MLR with differing inheritance models   ###
###################################################


models = c( 
  "MTL ~ Age ", 
  "MTL ~ Age + Sex ", 
  "MTL ~ Age + mMTL ", 
  "MTL ~ Age + fMTL ", 
  "MTL ~ Age + mMTL + fMTL ",
  "MTL ~ Age + mMTL * fMTL ",
  "MTL ~ mMTL ", 
  "MTL ~ fMTL ", 
  "MTL ~ mMTL * fMTL ", 
  "MTL ~ Age + MAC ", 
  "MTL ~ Age + PAC ", 
  "MTL ~ Age + MAC * PAC ", 
  "MTL ~ Age + mMTL + fMTL + MAC + PAC ", 
  "MTL ~ Age + mMTL + fMTL + MAC + PAC + Sex", 
  "MTL ~ Age + mMTL + fMTL + MAC ", 
  "MTL ~ Age + mMTL + fMTL + PAC ", 
  "MTL ~ Age + mMTL * fMTL + MAC + PAC ", 
  "MTL ~ Age + mMTL * fMTL + MAC ", 
  "MTL ~ Age + mMTL * fMTL + PAC ", 
  "MTL ~ Age + mMTL * fMTL + MAC * PAC ", 
  "MTL ~ Age + mMTL + fMTL + MAC * PAC ", 
  "MTL ~ MAC ", 
  "MTL ~ PAC ", 
  "MTL ~ MAC + PAC ", 
  "MTL ~ MAC * PAC ", 
  "MTL ~ Age * mMTL + Age * fMTL ",
  "MTL ~ Age + mMTL + fMTL + MAC + PAC + Age : mMTL + Age : fMTL ", 
  "MTL ~ Age + mMTL + fMTL + MAC + PAC + Age : mMTL + Age : fMTL + MAC * PAC ")

sf.list = vector(mode = "list", length = length(models))
fit.list = vector(mode = "list", length = length(models))

names(sf.list) = models
names(fit.list) = models

for(i in 1:length(models)){
  model = models[i]
  fit = lm(model, child.df)
  sf = summary(fit)
  sf.list[[i]] = sf
  fit.list[[i]] = fit
}


# summarize data into a rsqr matrix

fit.mat.cols = c("R.adj", "b1.Age", "b1.mMTL", "b1.fMTL", 
                 "b1.MAC", "b1.PAC")
rsqr.mat = matrix(nrow = length(fit.list), ncol = length(fit.mat.cols))
rownames(rsqr.mat) = names(sf.list)
colnames(rsqr.mat) = fit.mat.cols
pval.symbol = function(pval) {
  sig = ""
  if(pval < 0.001)
    sig = "***"
  else if(pval < 0.01)
    sig = "**"
  else if(pval < 0.05)
    sig = "*"
  return(sig)
}

for(i in 1:length(sf.list)){
  rsqr.mat[i,"R.adj"] = round(sf.list[[i]]$adj.r.squared, 3)
  if(grepl(" Age ", names(sf.list)[i])){
    var = "Age"
    rsqr.mat[i, "b1.Age"] = paste(pval.symbol(sf.list[[i]]$coefficients[var,4]), 
                                  round(fit.list[[i]]$coefficients[var],2))
  }
  if(grepl(" mMTL ", names(sf.list)[i])){
    var = "mMTL"
    rsqr.mat[i, "b1.mMTL"] = paste(pval.symbol(sf.list[[i]]$coefficients[var,4]), 
                                   round(fit.list[[i]]$coefficients[var],2))
  }
  if(grepl(" fMTL ", names(sf.list)[i])){
    var = "fMTL"
    rsqr.mat[i, "b1.fMTL"] = paste(pval.symbol(sf.list[[i]]$coefficients[var,4]), 
                                   round(fit.list[[i]]$coefficients[var],2))
  }
  if(grepl(" MAC ", names(sf.list)[i])){
    var = "MAC"
    rsqr.mat[i, "b1.MAC"] = paste(pval.symbol(sf.list[[i]]$coefficients[var,4]), 
                                    round(fit.list[[i]]$coefficients[var],2))
  }
  if(grepl(" PAC ", names(sf.list)[i])){
    var = "PAC"
    rsqr.mat[i, "b1.PAC"] = paste(pval.symbol(sf.list[[i]]$coefficients[var,4]), 
                                    round(fit.list[[i]]$coefficients[var],2))
  }
}


View(rsqr.mat)


# summarize best model 

best.fit  = fit.list[["MTL ~ Age + mMTL + fMTL + MAC + PAC "]]
sf = summary(best.fit)


# plot residuals predictions on the best model

plot(predict.lm(best.fit), sf$residuals)
abline(h = 0)
fit = lm(child.df$MTL ~ predict.lm(best.fit))
sf = summary(fit)

plot(child.df$MTL, predict.lm(best.fit), 
     xlim = c(3500, 12000), ylim  = c(3500, 12000),
     xlab = "Child MTL", ylab = "Predicted MTL", 
     pch = 19, col = "slategray4", cex= 1.5,
     main = "Best fit model: MTL ~ Age + mMTL + fMTL + MAC + PAC")

abline(fit$coefficients[1:2], col = "plum4", lwd = 3) 

legend(x = "topleft", cex = 1.5,
       legend = c(paste0("Adjusted R2 = ", round(sf$adj.r.squared, 2)), 
                  paste0("p value < 2.2e-16"), 
                  "mMTL*** fMTL*** MAC** Age*"),
       bty = "n")



###################################################
###       Partial correlations matrix           ###
###################################################

library(ppcor)

pcor.result = pcor(child.df[, c("MTL", "Age", "mMTL", "fMTL", "MAC", "PAC")])

pcor.mat = pcor.result$estimate
for(i in 1:nrow(pcor.mat)){
  for( j in 1:ncol(pcor.mat)){
    sig = pval.symbol(pcor.result$p.value[i,j])
    pcor.mat[i,j] = paste(sig, round(pcor.result$estimate[i,j], 2))
  }
}


 
View(pcor.mat)
