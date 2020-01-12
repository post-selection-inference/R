########################################################

# This script performs regression analysis for each SNV
# The response variable in the regression models is child's MTL
# the predictors are derived from the MLR inheritance analysis
# and combined with 0, 1, 2 values for each SNV
# (depending on the number of minor alleles in a genotype): 
#
# MTL ~ Age + mMTL + fMTL + MAC + PAC + SNV

# The script performs the calculations in chunks of SNVs for parallelization
# The outputs from different chunks should be merged manually at the end

########################################################





#########################################
###         Load data                 ###
#########################################


args = commandArgs(trailingOnly = T)

cat("Input arguments: ", args, "\n")

child.df.f = args[1] # provided as: "child.df.RData"
load(child.df.f)

ped.f = args[2] # ped file with SNV genotypes - available upon request
cat("Input ped file ", ped.f, "\n")

map.f = args[3] # map file with SNV descriptions - available upon request
cat("Input map file ", map.f, "\n")

from = args[4] # SNV range startpoint, used to process data in chunks
to = args[5] # SNV range endpoint, used to process data in chunks
cat("SNPs start and end specified as ", from, "-", to, "\n")

out.file = args[6]
cat("Output file specified as ", out.file, "\n")

cat("Reading map file ", map.f, "\n")
map = read.table(map.f, sep = "\t")




#########################################
###           FUNCTIONS                ##
#########################################

call.cmd <- function(command){
  OS.name = tolower(Sys.info()["sysname"])
  if (OS.name == "windows")
    shell(command, wait = T)
  else 
    system(paste("bash -c", "'", command, "'"), wait = T)
}

snp.lm = function(child.df, snp.ped, snp.ind){
  # performs linear regression for the SNV specified at the snp.ind index
  
  
  # subset ped file for children only
  snp.ped.ch = snp.ped[child.df$GoNL_Sample_ID,]
  
  # determine children genotypes 
  genotype = c()
  na.ind = c()

  for(i in 1:nrow(snp.ped.ch)){
  
    if(sum(snp.ped.ch[i, c(snp.ind, snp.ind + 1)] == 0) > 0)
      na.ind = c(na.ind, i)
    else 
      genotype = c(genotype, sum(snp.ped.ch[i, c(snp.ind, snp.ind + 1)] == 1))
  
  }
  
  # attach the genotypes to medatada 
  if(length(na.ind) > 0)
    child.snp.df = cbind(child.df[-na.ind,], genotype)
  else
    child.snp.df = cbind(child.df, genotype)
  
  # perform regression
  fit.snp = lm("MTL ~ Age + mMTL + fMTL + mAgeC + fAgeC + genotype", child.snp.df) 
  
  return(fit.snp)
}


get.lm.snp.mat = function(map, ped.f, child.df, 
                          out.file, chunk.size = 100, 
                          base.rsq, base.mMTL.est, base.fMTL.est, 
                          from=0, to=0){
  # performs regression for each SNV in the given range 
  # Returns: 
  #   A matrix of regression coefficients for each SNV (per row)
  #   Also includes difference of coefficient values with or without 
  #   SNV inclusion
  
  out.con = file(out.file, open = "w")
  
  # names of coefficients 
  lm.snp.names = c("Rsq", "delta.Rsq", 
                   "mMTL.est", "mMTL.delta.est", "mMTL.pval", 
                   "fMTL.est", "fMTL.delta.est", "fMTL.pval", 
                   "SNP.est", "SNP.pval")
  
  # determine chunk boundaries
  num.snp = nrow(map)
  if(from != 0)
    chunk.start = from 
  else
    chunk.start = 1
  
  if(to !=0)
    end = to
  else 
    end = num.snp
  
  
  while(chunk.start < end){
    
    cat("Processed ", chunk.start - 1, " snps out of ", end, "\n")
    
    if(chunk.start + chunk.size - 1 <= end)
      this.chunk.size = chunk.size 
    else 
      this.chunk.size = end - chunk.start + 1
    
    chunk.end = chunk.start + this.chunk.size - 1
    
    # create empty matrix
    lm.snp.mat = matrix(nrow = this.chunk.size, ncol = length(lm.snp.names))
    rownames(lm.snp.mat) = map$V2[chunk.start:chunk.end]
    colnames(lm.snp.mat) = lm.snp.names
    
    # subset the ped file according to chunk boundaries 
    temp.file = paste0("temp.", from, "-",to, ".ped")
    
    ped.command = paste0("cat ", ped.f, " | cut -d \" \" -f ", 
                         "2,", 5+chunk.start*2, "-", 6+chunk.end*2, " > ", temp.file) 
    cat("Redirecting ped chunk to file with command\n", ped.command, "\n")
    
    call.cmd(ped.command)
    
    ped = read.table(temp.file, sep = " ", row.names = 1)
    
    
    # perform regression for each SNV (call to snp.lm function)
    # populate the matrix of coefficients
    for(i in 1:this.chunk.size){
      
      sf.snp = summary(snp.lm(child.df, ped, -1 + i*2))
      
      lm.snp.mat[i, "Rsq"] = sf.snp$adj.r.squared
      lm.snp.mat[i, "delta.Rsq"] = sf.snp$adj.r.squared - base.rsq
      
      lm.snp.mat[i, "mMTL.est"] = sf.snp$coefficients[3,1]
      lm.snp.mat[i, "mMTL.delta.est"] = sf.snp$coefficients[3,1] - base.mMTL.est
      lm.snp.mat[i, "mMTL.pval"] = sf.snp$coefficients[3,4]
      
      lm.snp.mat[i, "fMTL.est"] = sf.snp$coefficients[4,1]
      lm.snp.mat[i, "fMTL.delta.est"] = sf.snp$coefficients[4,1] - base.fMTL.est
      lm.snp.mat[i, "fMTL.pval"] = sf.snp$coefficients[4,4]
      
      if (nrow(sf.snp$coefficients) > 6){ 
        lm.snp.mat[i, "SNP.est"] = sf.snp$coefficients[7,1]
        lm.snp.mat[i, "SNP.pval"] = sf.snp$coefficients[7,4]
      } else {
        lm.snp.mat[i, "SNP.est"] = NA
        lm.snp.mat[i, "SNP.pval"] = NA
      }
    } 
    
    # write the header if the this is the first chunk 
    if(chunk.start == 1)
      write.table(lm.snp.mat, file = out.con, sep = "\t", col.names = NA, quote = F, append = F)
    else 
      write.table(lm.snp.mat, file = out.con, sep = "\t", col.names = F, quote = F, append = T)
    chunk.start = chunk.start + this.chunk.size 
  }
  
  cat("Finished! Results written to ", out.file, "\n")
  close(out.con)
}



#########################################
###            PIPELINE               ###
#########################################

#  Base regression without the SNV 
# (to test the effect introduced by each SNV)

fit = lm("MTL ~ Age + mMTL + fMTL + mAgeC + fAgeC ", child.df)
sf = summary(fit)
base.rsq = sf$adj.r.squared
base.mMTL.est = sf$coefficients["mMTL",1]
base.fMTL.est = sf$coefficients["fMTL",1]


# Perform the regression
cat("Starting regression...\n")
get.lm.snp.mat(map, ped.f, child.df, out.file, chunk.size = 400, 
               base.rsq, base.mMTL.est, base.fMTL.est, 
               as.numeric(from), as.numeric(to))


