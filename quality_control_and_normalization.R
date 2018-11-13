# Perform QC and normalization on data
# First, analyze the data with the Illumina BeadArray Controls Reporter and remove any samples according to that tool
# Secondly, visually inspect the beta value distributions (we used GenomeStudio) and remove any samples that appear to have an abnormal beta value distribution

# load the data into lumi and minfi objects in order to perform various normalization and quality control procedures
library(lumi)
library(minfi)
library(ggplot2)
library(ggrepel)

# Load the data into a MethyLumiM object
master.lumiMethy <- lumiMethyR("FinalReport.txt")
# Load your sample information for this experiment (the sample annotation should be in the same order as the methylation data itself)
sample.annot <- read.delim("sampleAnnotation.txt",stringsAsFactors = F)

# Load the data into minfi object
targets <- read.metharray.sheet("./IDATfolder",
                                pattern="SampleSheet.csv")
rgSet <- read.metharray.exp(targets=targets)

# Plot global intensity values
MSet <- preprocessRaw(rgSet)
U.medians <- matrixStats::colMedians(getUnmeth(MSet))
M.medians <- matrixStats::colMedians(getMeth(MSet))
samplelabels <- sample.annot$patientID
samplelabels[U.medians>2000 | M.medians>2000] <- ""

ggplot(mapping=aes(M.medians,U.medians,label=samplelabels,color=sample.annot$Sex))+geom_point(pch=19)+labs(x="Median Methylated Intensity",y="Median Unmethylated Intensity",title="Median Intensities for Each Sample")+
  geom_text_repel(cex=2.75,max.iter = 20000,min.segment.length = 0)+
  geom_segment(aes(x = 2000, y = 0, xend = 2000, yend = 2000,color=NULL))+geom_segment(aes(x = 0, y = 2000, xend = 2000, yend = 2000,color=NULL))+
  theme(plot.title = element_text(hjust = 0.5))+scale_color_discrete("Sex")

# Cluster the samples with MDS clustering
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
pal = gg_color_hue(2)
plotMDS(getM(MSet), top=10000, gene.selection="common",col=pal[as.factor(sample.annot$Sex)],labels=sample.annot$patientID,main="MDS of Raw Beta Values\nUsing 10,000 most variable probes")

# Cluster the samples using hierarchical clustering, using the same probe set as the MDS clustering
top.probes <- rowMeans((getM(MSet) - rowMeans(getM(MSet)))^2)
top.probes <- order(top.probes, decreasing = TRUE)
top.probes <- top.probes[1:10000]

hc <- plotSampleRelation(master.lumiMethy, method="cluster",subset=top.probes)
dend <- as.dendrogram(hc)
labels(dend) <- sample.annot$patientID[order.dendrogram(dend)]
labels_colors(dend) <- pal[factor(sample.annot$Sex[order.dendrogram(dend)])]
plot(dend,main="Clustering of Raw Beta Values on 10,000 Most Variable Probes")

# Check the predicted sex of the samples:
GMSet <- mapToGenome(MSet)
predicted.sex <- getSex(GMSet)

plot(predicted.sex$xMed, predicted.sex$yMed, type = "n", xlab = "X chr, median total intensity (log2)", 
     ylab = "Y chr, median total intensity (log2)",main="Sex Prediction: X and Y Median Intensities")
text(predicted.sex$xMed, predicted.sex$yMed, sample.annot$patientID, col = ifelse(sample.annot$Sex == 
                                                                                   "Male", "deepskyblue", "deeppink3"))

# Plot the beta value distributions before normalization:
raw.betas <- getBeta(MSet, type = "Illumina")
raw.betas <- melt(raw.betas,varnames = "row.names")
names(raw.betas)[2]<-"sample"
raw.betas$sample <- as.character(raw.betas$sample)
raw.betas$Sex <- sample.annot$Sex[match(raw.betas$sample,sample.annot$patientID)]
rm(raw.betas)

ggplot(raw.betas,aes(value,color=Sex))+geom_density(aes(group=sample))+labs(title="Raw Beta Values\n(with 100 offset for denominator)")+theme(plot.title = element_text(hjust = 0.5))

# Plot the CpG Intensity distribution before normalization:
density(estimateIntensity(master.lumiMethy),xlab="log2(CpG-site Intensity)", main="Density of CpG Intensity",addLegend=F,col=pal[as.factor(sample.annot$Sex)])

# Plot the sample detection rate
detP <- detectionP(rgSet)
sample.detection <- apply(detP,2,function(x) sum(x>0.01))
num.detected.sites <- apply(detP,2,function(x) sum(x<=0.01))
sample.detection.rate <- num.detected.sites/nrow(detP)*100
tempdodge <- position_jitter(seed = 1,width = 0.1)

ggplot(mapping = aes(x=sample.annot$Sex,y=sample.detection.rate))+geom_boxplot(outlier.shape = NA)+labs(x="Sex",title="Sample Detection Rate\n(0.01 Detection p-value threshold)",y="Sample Detection Rate") +
  theme(plot.title = element_text(hjust = 0.5))+geom_point(pch=21,position = tempdodge)

# Remove all probes that are missing data in at least one sample
keep.sites <- apply(detP,1,function(x) sum(x>=0.01)==0)
keep.probe.names <- row.names(detP)[keep.sites]
final.lumiMethy <- master.lumiMethy[featureNames(master.lumiMethy) %in% keep.probe.names |
                                        grepl("rs",featureNames(master.lumiMethy)),]

# normalize the data using the lumi package
final.lumiMethy <- lumiMethyC(final.lumiMethy)
final.lumiMethy <- lumiMethyB(final.lumiMethy, method="bgAdjust2C")
final.lumiMethy <- lumiMethyN(final.lumiMethy, method='quantile')

# Plot the beta value distributions after normalization:
norm.betas <- methylated(final.lumiMethy)/(methylated(final.lumiMethy)+unmethylated(final.lumiMethy)+100)
norm.betas <- melt(norm.betas,varnames = "row.names")
names(norm.betas)[2]<-"sample"
norm.betas$sample <- as.character(norm.betas$sample)
norm.betas$Sex <- sample.annot$Sex[match(norm.betas$sample,sample.annot$patientID)]

ggplot(norm.betas,aes(value,color=Sex))+geom_density(aes(group=sample))+labs(title="Normalized Beta Values\n(with 100 offset for denominator)")+theme(plot.title = element_text(hjust = 0.5))

# Plot the CpG Intensity distribution after normalization:
density(estimateIntensity(final.lumiMethy),xlab="log2(CpG-site Intensity)", main="Density of CpG Intensity after Normalization",addLegend=F,col=pal[as.factor(sample.annot$Sex)])

# DNA methylation age prediction by Horvath
# Must decide if to put into script and how

# Compare SNP probes to genotypes? Might complicate the script, but certainly possible
norm.betas <- methylated(final.lumiMethy)/(methylated(final.lumiMethy)+unmethylated(final.lumiMethy)+100)
snp.methy.probes <- norm.betas[grepl("rs",row.names(norm.betas)),]
library(IlluminaHumanMethylationEPICmanifest)
data(IlluminaHumanMethylationEPICmanifest)
snpi.probes <- getProbeInfo(IlluminaHumanMethylationEPICmanifest,type=c("SnpI"))
snpi.probes$Infinium_Design_Type <- "I"
snpii.probes <- getProbeInfo(IlluminaHumanMethylationEPICmanifest,type=c("SnpII"))
snpii.probes$Infinium_Design_Type <- "II"
snpi.probes <- snpi.probes[,c("Name","Infinium_Design_Type")]
snpii.probes <- snpii.probes[,c("Name","Infinium_Design_Type")]
snp.probs <- rbind(snpi.probes,snpii.probes)

# load genotype data
library(GWASTools)
gds <- GdsGenotypeReader("ADHD_GWAS_genotypes_new_manifest_cleaned.gds")
sample.annotation <- read.delim("ADHD_GWAS_Sample_Annotation_Cleaned.txt",stringsAsFactors = F)
snp.annotation <-read.delim("ADHD_GWAS_SNP_Annotation_Cleaned.txt",stringsAsFactors = F)

scanAnnot <- ScanAnnotationDataFrame(sample.annotation)
snpAnnot <- SnpAnnotationDataFrame(snp.annotation)

genoData <- GenotypeData(gds,snpAnnot,scanAnnot)
sum(snpAnnot$snpName %in% row.names(snp.methy.probes))
geno <- getGenotypeSelection(genoData,snp=c(snpAnnot$snpName %in% row.names(snp.methy.probes)),scan=c(scanAnnot$subjectID %in% as.character(sample.annot$patientID)),sort=F)
row.names(geno) <- snpAnnot$snpName[match(row.names(geno),as.character(snpAnnot$snpID))]
colnames(geno) <- scanAnnot$subjectID[match(colnames(geno),as.character(scanAnnot$scanID))]

geno.ordered <- geno[,match(as.character(sample.annot$patientID),colnames(geno))]

snp.methy.probes.ordered <- snp.methy.probes[match(row.names(geno.ordered),row.names(snp.methy.probes)),!is.na(colnames(geno.ordered))]
sample.annot.ordered <- sample.annot[!is.na(colnames(geno.ordered)),]
snp.probs.ordered <- snp.probs[match(row.names(geno.ordered),snp.probs$IlmnID),]
geno.ordered <- geno.ordered[,!is.na(colnames(geno.ordered))]

# ensure that all of the genotypes, methylation data, and annotation are in the same order
dim(snp.methy.probes.ordered)==dim(geno.ordered)
all.equal(row.names(snp.methy.probes.ordered),row.names(geno.ordered))
all.equal(as.character(sample.annot.ordered$patientID),colnames(geno.ordered))
all.equal(sample.annot.ordered$label,colnames(snp.methy.probes.ordered))

# plot the methylation at each probe and the genotyped data:
par(mai = c(1,1,1.1,1), pin = c(5, 5)) 
for (i in 1:ncol(snp.methy.probes.ordered)) {
  plot.title <- paste(paste0("Patient ",sample.annot.ordered$patientID[i]),"red = Infinium I","black = Infinium II",sep="\n")
  temp <- snp.methy.probes.ordered[,i]
  temp[snp.probs.ordered$Infinium_Design_Type=="I"] <- 1 - temp[snp.probs.ordered$Infinium_Design_Type=="I"]
  plot(2-geno.ordered[,i],temp,pch=20,col=color.prob,xlim=c(0,2),ylim=c(0,1),
       main=plot.title,xlab="Genotype (number of G/C alleles)",ylab="Methylation\n(high methylation corresponds to G/C alleles)")
}

# calculate correlations between methylation of SNP probes and genotypes:
corr.dat <- sapply(1:ncol(snp.methy.probes.ordered), function(i) { 
  temp <- snp.methy.probes.ordered[,i]
  temp[snp.probs.ordered$Infinium_Design_Type=="I"] <- 1 - temp[snp.probs.ordered$Infinium_Design_Type=="I"]
  cor(temp,2-geno.ordered[,i],use = "complete.obs")
})

hist(corr.dat,main="Histogram of correlations between Methylation QC Probes\nand true genotypes",xlab="Correlation",ylab="Number of Samples",breaks=20)
corr.dat[corr.dat<0.8]
sum(is.na(corr.dat))
which(is.na(corr.dat))
which(corr.dat<0.8)
# remove any samples for which the SNP and methylation SNP probes to do not agree

# Save your data
save(final.lumiMethy,sample.annot,file = "normalized_data.Rdata")








