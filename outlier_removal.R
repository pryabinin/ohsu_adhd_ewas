# load data from probe filtering script
load("filtered_norm_data.Rdata")
beta.values <- methylated(final.lumiMethy)/(methylated(final.lumiMethy)+unmethylated(final.lumiMethy)+100)

# Find outlier samples and probes
# Outlier samples:
num.outlier.samples <- sapply(1:ncol(beta.values),function(x) {
  diff.samples <- sample((1:ncol(beta.values))[-x],100)
  median.diff <- rep(NA,ncol(beta.values))
  for (i in diff.samples) {
    temp <- abs(beta.values[,x] - beta.values[,i])
    temp <- sum(temp>0.4,na.rm = T)
    median.diff[i] <- temp
  }
  return(median(median.diff,na.rm=T))
})

hist(num.outlier.samples,main="Distribution of Median Counts of Differences > 0.4 between 100 Unrelated Pairings for each Samples",
     xlab="Median Counts of Difference > 0.4", "Number of Samples")

# Outlier Probes:
num.outlier.probes <- sapply(1:nrow(beta.values),function(x) {
  diff.samples.1 <- sample((1:ncol(beta.values)),100)
  diff.samples.2 <- sample((1:ncol(beta.values)),100)
  while(sum(diff.samples.1==diff.samples.2)>0)
    diff.samples.2 <- sample((1:ncol(beta.values)),100)
  median.diff <- rep(NA,nrow(beta.values))
  for (i in diff.samples) {
    temp <- abs(beta.values[x,diff.samples.1[i]] - beta.values[x,diff.samples.2[i]])
    temp <- sum(temp>0.4,na.rm = T)
    median.diff[i] <- temp
  }
  return(median(median.diff,na.rm=T))
})

hist(num.outlier.probes,main="Distribution of Number of Unrelated Sample Pairs (out of 100 pairs) that Differed by > 0.4 for a Probe",
     xlab="Number of Unrelated Sample Pairs Differing by > 0.4", "Number of Probe")

# Remove any outlier samples and outlier probes

save(beta.values,sample.annot,file = "filtered_norm_beta_values.Rdata")





