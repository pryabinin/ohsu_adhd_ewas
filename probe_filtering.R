# load data from quality control and normaliation
load("normalized_data.Rdata")

# remove multimapping and cross reactive probes
# is it OK to include someone else's results? If not, how to include these results in our github?
evans.filt <- read.delim("mmc2.txt",header=F,stringsAsFactors = F)
evans.filt.2 <- read.delim("mmc3.txt",header=F,stringsAsFactors = F)

keep <- !(featureNames(final.lumiMethy) %in% c(evans.filt$V1,evans.filt.2$V1))
table(keep)
final.lumiMethy <- final.lumiMethy[keep,]

# remove probes which have a SNP with MAF >= 1% according to the Evans paper:
evans.snps <- read.delim("Evans_probe_SNPs.txt",stringsAsFactors = F)
keep <- !(featureNames(final.lumiMethy.all.filt) %in% evans.snps$IlmnID[evans.snps$Max_MAF >= 0.01])
table(keep)
final.lumiMethy <- final.lumiMethy[keep,]

# remove probes which have a MAF >= 1% according to the Illumina manifest:
manifest.file <- "MethylationEPIC_v-1-0_B2.csv"
manifest <- read.csv(manifest.file,skip=7,header = T,stringsAsFactors = F)
manifest <- manifest[1:866895,]

manifest$Max_SNP_MinorAlleleFrequency <- sapply(sapply(strsplit(as.character(manifest$SNP_MinorAlleleFrequency),split = ";"),as.numeric),max)
keep <- !(featureNames(final.lumiMethy) %in% manifest$IlmnID[manifest$Max_SNP_MinorAlleleFrequency >= 0.01])
table(keep)
final.lumiMethy <- final.lumiMethy[keep,]

# Remove SNP probes
keep <- !(featureNames(final.lumiMethy.all.filt) %in% manifest$IlmnID[grepl("rs",manifest$IlmnID)])
table(keep)
final.lumiMethy <- final.lumiMethy[keep,]

# Remove probes that are not in the newest manifest
manifest.b4 <- read.csv("MethylationEPIC_v-1-0_B4.csv",stringsAsFactors=F,skip = 7)
manifest.b4 <- manifest.b4[manifest.b4$Infinium_Design_Type!="",]
keep <- featureNames(final.lumiMethy) %in% manifest.b4$IlmnID
table(keep)
final.lumiMethy <- final.lumiMethy[keep,]

# Remove non-autosomal probes
keep <- !(featureNames(final.lumiMethy) %in% manifest$IlmnID[!(manifest$CHR %in% c("X","Y"))])
table(keep)
final.lumiMethy <- final.lumiMethy[keep,]

save(final.lumiMethy,sample.annot,file = "filtered_norm_data.Rdata")