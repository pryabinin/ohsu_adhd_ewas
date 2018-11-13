library(RefFreeEWAS)
library(ggplot2)
library(reshape2)
library(limma)

load("filtered_norm_beta_values.Rdata")

# Create Cell Mixture Models
data.sd <- apply(beta.values,1,sd)
names(data.sd) <- 1:length(data.sd)
most.variable.idx <- as.integer(names(sort(data.sd,decreasing = T)))[1:20000]
beta.value.var <- beta.values[most.variable.idx,]

rf.nnmf  <- RefFreeCellMixArray(as.matrix(beta.value.var),Klist=1:10,verbose=T,iters=25)

# Run Bootstraps
rf.deviance <- RefFreeCellMixArrayDevianceBoots(rf.nnmf,Y = as.matrix(beta.value.var),R = 900,bootstrapIterations = 10)

# Graph the bootstraps:
dat <- rbind(rf.deviance,apply(rf.deviance,2,median))
dat <- as.data.frame(dat)
dat$name <- c("full data set deviance",1:nrow(all.bootstraps),"median boostrap deviance")
dat$color <- c(alpha("black",1),rep(alpha("green",.1),nrow(all.bootstraps)),alpha("red",1))
dat <- melt(dat,id.vars=c("name","color"),value.name = "deviances",variable.name="K")
p <- ggplot(dat, aes(K, deviances)) + labs(title="Deviances of cell-type factorization models",x="K (number of cell types)",y="deviance")+ geom_line(aes(colour=color,group = name)) +
  scale_color_manual(name="Legend",values=c(alpha("black",1), alpha("green",.05),alpha("red",1)),labels=c("Full data set","Bootstraps","Median Deviance\nof Bootstraps")) +
  guides(colour = guide_legend(override.aes = list(alpha = 1)))
p

# Adjust the beta values according to the value that minimizes the median bootstrapped deviance
khat.dat <- rf.nnmf[[10]]
cell.profiles <- khat.dat$Mu
cell.prop <- khat.dat$Omega
sample.means <- apply(cell.prop, 2, mean)
sample.means

centered.cell.prop <- sweep(cell.prop, 2, sample.means, "-")
apply(centered.cell.prop, 2, mean) 

X.cellprop <- model.matrix(~centered.cell.prop) 
lm.betas <- lmFit(beta.values, X.cellprop)

adj.betas <- beta.values - lm.betas$coef[,-1] %*% t(X.cellprop[,-1])
sum(adj.betas>=1 | adj.betas<=0)
trunc.snps <- row(adj.betas)[which(adj.betas>=1 | adj.betas<=0)]
trunc.samps <- col(adj.betas)[which(adj.betas>=1 | adj.betas<=0)]
# plot the number of beta values that lie outside of the 0-1 range
qplot(as.vector(table(trunc.snps)), geom="bar",main="Number of times a CpG needed to be truncated (more than 0)",xlab="Number of times a CpG's value needed to be truncated",ylab="Number of CpGs")
qplot(as.vector(table(trunc.samps)), geom="bar",main="Number of times a sample needed to be truncated (more than 0)",xlab="Number of times a sample's value needed to be truncated",ylab="Number of samples")

# truncate the beta values that lie outside of the range 0-1
adj.betas[adj.betas>=1] <- 0.999999
adj.betas[adj.betas<=0] <- 0.000001

# plot the final beta value distribution
qplot(value, data=melt(t(adj.betas)), geom="histogram",main = "Beta values adjusted for cell-type proportions",xlab = "Beta")

save(adj.betas,sample.annot,file = "adjusted_betas.Rdata")






