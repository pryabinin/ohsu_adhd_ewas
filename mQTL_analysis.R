# Setup and Perform mQTL analysis using the GEM package
library(data.table)
# Load cell type adjusted beta values
load("adjusted_betas.Rdata")
sample.annot <- sample.annot[match(colnames(adj.betas),sample.annot$patientID),]

# Create text files used by GEM:
# convert sex variable from string to integer
temp.sex <- ifelse(sample.annot$sex)

# create covariate file for GEM
cov.df <- as.data.frame(t(data.frame(ID=sample.annot$patientID,
                                     sex=sample.annot$sex,
                                     age=sample.annot$age,
                                     pc1=sample.annot$PC1,
                                     pc2=sample.annot$PC2,
                                     pc3=sample.annot$PC3)))

write.table(cov.df,"cov_file.txt",quote = F,row.names = T,col.names = F,sep="\t")

# create methylation file for GEM
meth.data <- adj.betas[row.names(adj.betas) %in% top15[,1],]
meth.data <- rbind(sample.annot$patientID,meth.data)
meth.data <- cbind(c("ID",row.names(meth.data)[-1]),meth.data)
meth.data$`c("ID", row.names(meth.data)[-1])` <- as.character(meth.data$`c("ID", row.names(meth.data)[-1])`)
write.table(meth.data,file="methyl_file.txt",quote = F,row.names = F,col.names = F,sep="\t")

# create genotype file for GEM
library(GWASTools)
gds <- GdsGenotypeReader("ADHD_GWAS_genotypes_new_manifest_cleaned.gds")
sample.annotation <- read.delim("ADHD_GWAS_Sample_Annotation_Cleaned.txt",stringsAsFactors = F)
snp.annotation <-read.delim("ADHD_GWAS_SNP_Annotation_Cleaned.txt",stringsAsFactors = F)

scanAnnot <- ScanAnnotationDataFrame(sample.annotation)
snpAnnot <- SnpAnnotationDataFrame(snp.annotation)

genoData <- GenotypeData(gds,snpAnnot,scanAnnot)

all.geno <- getGenotypeSelection(genoData,scanID=as.integer(as.character(meth.data[1,-1])))
row.names(all.geno) <- snpAnnot$snpName
all.geno <- all.geno+1
all.geno <- as.data.table(all.geno,keep.rownames = T)
names(all.geno)[1] <- "ID"
setcolorder(all.geno, c(1,match(as.character(meth.data[1,-1]),colnames(all.geno))))
all.equal(as.character(meth.data[1,]),names(all.geno))
fwrite(all.geno,file="geno_file.txt",quote=F,sep="\t",col.names=T,row.names=F)
close(genoData)

# Perform the mQTL analysis using GEM:
library(GEM)

GEM_Gmodel(snp_file_name = "geno_file.txt",
           covariate_file_name = "cov_file.txt",
           methylation_file_name = "methyl_file.txt",
           Gmodel_pv = 1e-5,
           output_file_name = "G_model_results.txt")

