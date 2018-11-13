This folder describes the quality control, normalization, and analyses in this experiment. The scripts describe the workflow, note that the scripts will not run successfully since they reference data files that are not provided.

The order of the workflow is:
quality_control_and_normalization.R
probe_filtering.R
outlier_removal.R
cell_correction.R
full_model_analysis.R
mQTL_analysis.R

Script explanations:
quality_control_and_normalization.R
	This script goes over quality control checks for filtering out poorly performing samples and probes, as well as the normalization method we used on this experiment.[1,2

probe_filtering.R
	This script removes probes for the following reasons:
		1. Probe was multimapping or cross reactive
		2. Probe contained a SNP underlying it with a MAF >= 1%
		3. Probe was a SNP interrogating probe
		4. Probe does not exist in the newest version (v4) of the Illumina manifest
		5. Probe is non-autosomal

outlier_removal.R
	Detect outlier samples by calculating the median number of probes which differed by more than 0.4 (40% beta difference) between each sample and 100 other samples in the data set
	Detect outlier probes by calculating the median number of samples which differed by more than 0.4 (40% beta difference) between 100 random sample pairings in the data set

cell_correction.R
	Utilize reference-free cell type prediction implemented in the R package RefFreeEWAS. First, cell mixture models for 1 to 10 different cell types are generated. Then, bootstrap models and subsequent deviances are calculated for each model. The model that the minimizes the median of the bootstrapped deviances is selected as the optimal model to proceed with.
	Adjust the beta values by first calculating the linear model: Beta Values ~ 1 + B * centered cell proportions
		and then subsequently calculate the adjusted beta values: Adjusted Beta Values = Beta Value - B * centered cell proportions

full_model_analysis.R
	Perform probe level differential methlyation analysis on the cleaned and cell type adjusted beta values, using the limma framework.
	The model calculated in this analysis is:
		Adjusted Beta Value ~ 1 + ADHD status + Sex + Age + Genomic PCs1-3

mQTL_analysis.R
	Perform mQTL analysis using the GEM framework for the top most differentially methylation. First, text files for use in GEM are generated (covariate file, genotype file, and methylation file). Then, run the GEM G model using the generated files.

References:
	
[1] Du, P., Kibbe, W.A., Lin, S.M. (2008). “lumi: a pipeline for processing Illumina microarray.” Bioinformatics.
[2] Aryee MJ, Jaffe AE, Corrada-Bravo H, Ladd-Acosta C, Feinberg AP, Hansen KD, Irizarry RA (2014). “Minfi: A flexible and comprehensive Bioconductor package for the analysis of Infinium DNA Methylation microarrays.” Bioinformatics, 30(10), 1363–1369.
[3] Gogarten SM, Bhangale T, Conomos MP, Laurie CA, McHugh CP, Painter I, Zheng X, Crosslin DR, Levine D, Lumley T, Nelson SC, Rice K, Shen J, Swarnkar R, Weir BS, Laurie CC (2012). “GWASTools: an R/Bioconductor package for quality control and analysis of genome-wide association studies.” Bioinformatics, 28(24), 3329-3331