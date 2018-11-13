library(limma)
# Load cell type adjusted beta values
load("adjusted_betas.Rdata")

# Perform analysis of full model:
# Methylation (beta values adjusted for cell type) ~ 1 + age + sex + status + 3 Genomic PC's
sample.annot <- sample.annot[match(colnames(adj.betas),sample.annot$patientID),]

temp.status <- factor(sample.annot$status, levels=c("control","case"))
temp.sex <- factor(sample.annot$sex, levels=c("Female","Male"))
temp.design <- model.matrix(~temp.status+temp.sex+sample.annot$Age+sample.annot$PC1+sample.annot$PC2+sample.annot$PC3)
results.fit <- lmFit(adj.betas, temp.design)
results.fit <- eBayes(results.fit)

# Save the results object
save(results.fit,file="full_model_fit_object.Rdata")

# Save the list of top 15 methylation probes, for use in the mQTL analysis
top15 <- row.names(topTable(results.fit,coef="temp.statuscase",number = 15))
write.table(top15,"top15_methyl_probes.txt",sep="\t",quote=F,row.names=F)
