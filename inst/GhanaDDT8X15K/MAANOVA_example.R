###########################################
# Differential Expression Testing Tutorial
# ============ ========== ======= ========
#
#           Ricardo A. Verdugo
#             Gary Churchill
#
#          JAX, Bar Harbor 2009
#
#
# Data description
# ---- -----------
# 
# This is an expression profile experiment done in the
# Illumina Mouse-Ref8 platform.
#
# Objective: To assess the effect of genetic variation
# in chromosome Y on mice on the size of cardiomyocytes.
# 
# Experimental design: Eight adult male mice from 
# two strains were profiled, C57BL/6J and 
# C57BL/6J-chrY<A/J/NaJ>, referred to as B and BY herein.
# From each strain (genotype), four animals were 
# castrated and four were sham operated.
#
# Aims: 
#       1) To determine differential expression between genotypes
#       2) To determina differential expression between treatments
#       3) To assess differences in the response to treatment between
#          the two genotypes
#
# References:
#
# For more information about the experimental samples:
#   Llamas, Bastien, Ricardo Verdugo, Gary Churchill, and Christian Deschepper. 2009.
#   Chromosome Y variants from different inbred mouse strains are linked to differences in the
#   morphologic and molecular responses of cardiac cells to postpubertal testosterone. BMC
#   Genomics 10, no. 1 (April 7): 150. doi:10.1186/1471-2164-10-150.
# 
# For informations about the analysis of this data: Verdugo, Ricardo A., Christian F.
#   Deschepper, Gloria Munoz, Daniel Pomp, and Gary A. Churchill. 2009. Importance of
#   randomization in microarray experimental designs with Illumina platforms. Nucl. Acids Res.
#   37, no. 17 (September): 5610-8. doi:10.1093/nar/gkp573.
# 
# To speed up the tutorial, only the first 3000 of the probes in the microarray
# are included. The full dataset is available at the GEO database, record GSE15354.
# http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15354


# Preliminaries
# =============

## Define some variables (it is a good practice to declare them at the beginning of the scripts)
outdir <- "output"
annot_file <- file.path(outdir, "MouseRef-8_annot.txt") # file.path() is a robust way to define file paths
# Proportion of false discoveries that is acceptable
fdr_th <- 0.2 

# Define some useful functions (it is a good idea to always document what the intend to do)
logdiff2FC=function(x, base=2) {
  # Transforms x from a log difference (logdiff) to a fold change (FC) scale
  # logdiff is the expression difference in log scale of base = 'base'
  # FC is the ratio between two expression values in the original scale
  # The absolute value of the returned FC is always greater than or equal to 1
  # The sign is kept to indicate direction of change
  # When the the log difference is 0, the FC is always 1.
  signs = sign(x)
  out   = base^abs(x)
  if(any(x==0)) {
    out[x==0] = 1
  }
  return(out*signs)
} # logdiff2FC

# Create an output directory
if(!file.exists(outdir)) {
  dir.create(outdir, mode = "0755", recursive=T)
}

# Readin raw data
# ====== === ====
Data.Raw=read.delim("Illum_data_sample.txt")
signal=grep("AVG_Signal", colnames(Data.Raw))
sdevs=grep("BEAD_STDERR", colnames(Data.Raw))
detection=grep("Detection.Pval", colnames(Data.Raw))
annot=c(1,2,99:125)

# Read hybridization design
design=read.csv("YChrom_design.csv")
# Lets take a look
print(design)

# Save annotations in case they are needed for other applications
# In Illumina experiment, the data comes with several columns of probe annotation
write.table(Data.Raw[,annot], file=annot_file, sep="\t", quote=F, col.names=T, row.names=F)


# Quality Control
# ======= =======

# Boxplots
palette(rainbow(4))
# Color by genotype
png(file.path(outdir,"boxplot_raw_genotype.png"), width=4, height=4, unit="in", res=150)
  par(xpd=NA, mar= c(6.1, 4.1, 4.1, 2.1), cex=.7)
  boxplot(as.data.frame(log2(Data.Raw[,signal])), horiz=T, main="Raw log2 values Boxplot", las=1, col=design$Genotype, names=design$Sentrix_Position, cex.axis=.9)
  legend(8, 2.5, legend=levels(design$Genotype), fill=1:2, ncol=2, xjust=.5)
dev.off()
# Color by treatment
png(file.path(outdir,"boxplot_raw_treatment.png"), width=4, height=4, unit="in", res=150)
  par(xpd=NA, mar= c(6.1, 4.1, 4.1, 2.1), cex=.7)
  boxplot(as.data.frame(log2(Data.Raw[,signal])), horiz=T, main="Raw log2 values Boxplot", las=1, col=design$Treatment, names=design$Sentrix_Position, cex.axis=.9)
  legend(8, 2.5, legend=levels(design$Treatment), fill=1:2, ncol=2, xjust=.5)
dev.off()

# Create matrix of raw data
rawdata=as.matrix(Data.Raw[,signal])
rownames(rawdata)=Data.Raw$PROBE_ID
colnames(rawdata)=design$Sample_Name

## Scatter plots of raw data
# Untransformed
png(file.path(outdir,"Pairs_scatter_raw.png"), width=8, height=8, unit="in", res=150)
  par(cex=.2, mar=c(2.1,2.1,2.1,1.1))
  pairs(rawdata, main="Raw Intensity Values", pch=".",  gap=.5, cex.labels=.5)
dev.off()
# log2 transformed
png(file.path(outdir,"Pairs_scatter_log2.png"), width=8, height=8, unit="in", res=150)
  par(cex=.2, mar=c(2.1,2.1,2.1,1.1))
  pairs(log2(rawdata), main="Log2 Raw Intensity Values", pch=".",  gap=.5, cex.labels=.5)
dev.off()


# Data Normalization
# ==== =============

# Load functions for normalization
library(affy)
library(preprocessCore)

normdata=normalize.quantiles(rawdata) 
colnames(normdata)=colnames(rawdata)
rownames(normdata)=rownames(rawdata)


# Probe Filtering
# ===== =========
# This step aims to removing probes that did not detect a transcript
# in any of the experimental groups. Note that this step can be optional.
#
# Create a vector or P/A calls for each probe using detection probabilities calculate by BeadStudio
present=which(apply(apply(Data.Raw[,detection]<.04, 1, tapply, design$Group, mean)>=.5, 2, any))


# Testing for differential expression
# ======= === ============ ==========

# Load the MAanova package
library(maanova)

# Create a madata object which only includes present detected transcripts
madata=read.madata(normdata[present,], design, log.trans=T)

# Some basis statistics on each experimental group
Means = t(apply(madata$data, 1, tapply, design$Group, mean)) 
colnames(Means)=paste("Mean", colnames(Means), sep=":")
SEs = t(apply(madata$data, 1, tapply, design$Group, function(x) sqrt(var(x)/length(x))))
colnames(SEs)=paste("SE", colnames(SEs), sep=":")

# Fit the model
fit.fix=fitmaanova(madata, formula=~Group)

# Construct a matrix of contrastas of interest
# --------------------------------------------
# In factorial design like this one, one can ask different questions from the data. Each
# question can be tested by a comparison between some set of experimental groups. These
# comparisons are called contrasts. The matest function from MAanova can take a matrix 
# of contrasts and test whether those comparisons explain a significant proportion of variance
# in the expression levels measured by each probe.
#
# But, before we create our matrix of contrasts, lets define some terms to simplify nomenclature.
# Let,
# I      : intact (not castrated) treatment
# C      : castrated treatment
# B      : C57BL/6J genotype
# B.Y    : C57BL/6J-chrY<A> genotype (chromosome Y congenic strain on a C57BL/6J genomic background)
# Geno   : genotype 
# Trt    : treatment
# Int    : genotype x treatment interaction
# Geno_I : genotype effect in I animals
# Geno_C : genotype effect in C animals
# Trt_B  : treatment effect in the B genotype
# Trt_BY : treatment effect in the B.Y genotype
#
# Then the four experimental groups can be denoted by B.C, B.I, BY.C, and BY.I.
# 
# Contrasts are vectors of coefficients that when multiplied to the vector of means by
# experimental group, create comparisons that relevant and that can be tested statistically.
# To create the vector of contrast coefficients, assume the experimental groups are sorted
# alphabetically.
# 
# Now, construct a matrix of contrasts where rows are contrasts and columns are experimental
# groups:
#                          B.C  B.I BY.C BY.I
cmat = rbind( Geno     =  c( 1,   1,  -1,  -1 )*.5,
              Trt      =  c( 1,  -1,   1,  -1 )*.5,
              Int      =  c( 1,  -1,  -1,   1 ),
              Geno_I   =  c( 0,   1,   0,  -1 ),
              Geno_C   =  c( 1,   0,  -1,   0 ),
              Trt_B    =  c( 1,  -1,   0,   0 ),
              Trt_BY   =  c( 0,   0,   1,  -1 ),
              B.C_BY.I =  c( 1,   0,   0,  -1 ),
              B.I_BY.C =  c( 0,   1,  -1,   0 ))
              
# We can use these contrasts for calculate some ratios (fold changes)
# that can be of interest.

# Log differences for main effects and the interaction
Geno = apply(Means, 1, function(x) sum(x*cmat[1,]))
Trt  = apply(Means, 1, function(x) sum(x*cmat[2,]))
Int  = apply(Means, 1, function(x) sum(x*cmat[3,]))

# Log differences for comparisons within factor levels
Geno_I   = apply(Means, 1, function(x) sum(x*cmat[4,]))
Geno_C   = apply(Means, 1, function(x) sum(x*cmat[5,]))
Trt_B    = apply(Means, 1, function(x) sum(x*cmat[6,]))
Trt_BY   = apply(Means, 1, function(x) sum(x*cmat[7,]))
B.C_BY.I = apply(Means, 1, function(x) sum(x*cmat[8,]))
B.I_BY.C = apply(Means, 1, function(x) sum(x*cmat[9,]))

# Bind columns into a matrix
logDiffs=cbind(Geno, Trt, Int, Geno_I, Geno_C, Trt_B, Trt_BY, B.C_BY.I, B.I_BY.C)

# Transform log differences to FC scale
FC = apply(logDiffs, 2, logdiff2FC)

# Test test each contrasts using 300 permutations of sample labels
test.cmat=matest(madata, fit.fix, term="Group", Contrast=cmat, n.perm=300, test.type = "ttest",
                 shuffle.method="sample", verbose=TRUE)

# Contrasts names are not kept in the matrix of permutation results, so lets 
# copy them from the matrix of tabular p-values
colnames(test.cmat$Fs$Pvalperm) = colnames(test.cmat$Fs$Ptab)

# Plot p-values comparing different ways of calculating them (see ?matest)
png(file.path(outdir,"P-values Hist.png"), width=6, height=6, unit="in", res=150)
  par(mfrow=c(2,2), oma=c(2,0,2,0), cex=.8, xpd=NA)
  palette(rainbow(3))
  plot(density(test.cmat$F1$Ptab[,1]), col=1, main="F1:Ptab", lwd=2)
  lines(density(test.cmat$F1$Ptab[,2]), col=2, lwd=2)
  lines(density(test.cmat$F1$Ptab[,3]), col=3, lwd=2)
  
  plot(density(test.cmat$F1$Pvalperm[,1]), col=1, main="F1:Pvalperm", lwd=2)
  lines(density(test.cmat$F1$Pvalperm[,2]), col=2, lwd=2)
  lines(density(test.cmat$F1$Pvalperm[,3]), col=3, lwd=2)
  
  plot(density(test.cmat$Fs$Ptab[,1]), col=1, main="Fs:Ptab", lwd=2)
  lines(density(test.cmat$Fs$Ptab[,2]), col=2, lwd=2)
  lines(density(test.cmat$Fs$Ptab[,3]), col=3, lwd=2)
  
  plot(density(test.cmat$Fs$Pvalperm[,1]), col=1, main="Fs:Pvalperm", lwd=2)
  lines(density(test.cmat$Fs$Pvalperm[,2]), col=2, lwd=2)
  lines(density(test.cmat$Fs$Pvalperm[,3]), col=3, lwd=2)
  
  legend(-.5, -1.6, legend=c("Geno", "Trt", "Int"), col=1:3, lwd=2, xjust=.5, ncol=3, xpd=NA)
dev.off()


# Multiple comparison control (FDR transformation, see ?adjPval)
test.cmat=adjPval(test.cmat, method="adaptive")

# Plot distribution of FDR values
png(file.path(outdir,"FDR Hist.png"), width=6, height=6, unit="in", res=150)
  par(mfrow=c(2,2), oma=c(2,0,2,0), cex=.8, xpd=NA)
  palette(rainbow(3))
  plot(density(test.cmat$F1$adjPtab[,1]), col=1, main="F1:Ptab", lwd=2, xlim=c(-.1,1.1))
  lines(density(test.cmat$F1$adjPtab[,2]), col=2, lwd=2)
  lines(density(test.cmat$F1$adjPtab[,3]), col=3, lwd=2)
  
  plot(density(test.cmat$F1$adjPvalperm[,1]), col=1, main="F1:Pvalperm", lwd=2, xlim=c(-.1,1.1))
  lines(density(test.cmat$F1$adjPvalperm[,2]), col=2, lwd=2)
  lines(density(test.cmat$F1$adjPvalperm[,3]), col=3, lwd=2)
  
  plot(density(test.cmat$Fs$adjPtab[,1]), col=1, main="Fs:Ptab", lwd=2, xlim=c(-.1,1.1))
  lines(density(test.cmat$Fs$adjPtab[,2]), col=2, lwd=2)
  lines(density(test.cmat$Fs$adjPtab[,3]), col=3, lwd=2)
  
  plot(density(test.cmat$Fs$adjPvalperm[,1]), col=1, main="Fs:Pvalperm", lwd=2, xlim=c(-.1,1.1))
  lines(density(test.cmat$Fs$adjPvalperm[,2]), col=2, lwd=2)
  lines(density(test.cmat$Fs$adjPvalperm[,3]), col=3, lwd=2)
   
  legend(-.5, -1.5, legend=c("Geno", "Trt", "Int"), col=1:3, lwd=2, xjust=.5, ncol=3)
dev.off()

# Summarize into a table of results for all present transcripts
# Here we are exporting only the tests from the permutations on the F values with shrinkage variance estimates,
# but other are available (see ?matest)
annot.out = c("SYMBOL", "ACCESSION", "REFSEQ_ID", "ENTREZ_GENE_ID", "PROTEIN_PRODUCT", "ONTOLOGY_COMPONENT", "ONTOLOGY_PROCESS", "ONTOLOGY_FUNCTION")
out = data.frame(Probe=rownames(rawdata)[present], Data.Raw[present, annot.out], Means, SEs, F_val=test.cmat$Fs$Fobs, P_val=test.cmat$Fs$Pvalperm, 
               FDR=test.cmat$Fs$adjPvalperm, FC=FC)

# Select probes with FDR <= 0.2 for any of the 9 tests
selected = rownames(rawdata)[present][apply(test.cmat$Fs$adjPvalperm <= fdr_th, 1, any)]

# Read probe annotations 
annot = read.delim("annot.txt")

# Merge to table of test results
out.annot = merge(out, annot, by.x="Probe", by.y="probe_name", all.x=T, all.y=F, sort=F)

# Export all results
write.table(out.annot , file=file.path(outdir,"YChrom_results_allpresent.csv"), sep=",", row.names=F)

# Export only selected
write.table(out[selected,], file=file.path(outdir, "YChrom_results_onlyselected.csv"), sep=",", row.names=F)

# Count genes differentially expressed
# ----- ----- -------------- ---------
# In this experiment, a question of interest was how many genes respond differently to the
# treatment of castration in the two genotypes. In other words, how important is the
# interaction between genotype and treatment. Secondly, it was interesting to assess the
# nature of the interaction. Are both genotypes responding to the treatment but in opposite
# directions? Or does the treatment have an effect in one of genotype and not in the other?
# To answer the first question, one could count the number of probes that show a significant
# effect for the Int contrast, i.e. how many probes have an FDR below a threshold for the
# Int.Pvalperm test? If you don;t know how to calculate this in R yet, open the exported
# YChrom_results_onlyselected.csv file in a spread sheet editor such as Calc (OpenOffice)
# and use filters to calculate this number. Then try to do this in R and compare the
# results.
# 
# To answer the second question, we can use Venn Diagrams. We want to count the number if
# genes that are selected in each genotype among those that show a significant interaction
# effect.
# 
# One caveat is that genes are represented by multiple probes in this microarray platform.
# It is a good idea to count genes only once, but there is a question of how to count a gene
# when different probes are giving different signals, i.e. one is saying that the gene is
# selected whereas the other is saying that is not. This platform was based RefSeq
# (http://www.ncbi.nlm.nih.gov/RefSeq/), which is a human-curated database of reference
# transcripts for known genes and it was designed to avoid redundancy
# (http://www.illumina.com/products/mouseref-8_expression_beadchip_kits_v2.ilmn). Therefore,
# we will assume that each probe is testing different biological signals, and not a repeated
# measure of the same transcript. We will count a gene as selected if any of it transcripts
# (probes) is selected. For other platforms where probes provide repeated measurements for
# the same transcript, one may want to use a voting or averaging approach to summarize
# results at the gene-level.

# First, eliminate probes without a gene annotation
out.annot=out.annot[!is.na(out.annot$gene_id),]

# Now count significant genes for each comparison. All comparison are done only for
# probes that are also significant for the Interaction term.
# Because were are interested in genes that show an interaction, we cannot the the Geno
# contrast to test for significance, since that contrast tests the marginal effect of 
# the genotype across treatments. In other words, it tests for effects that are consistent in both
# treatments. These can be zero even if the genotype has an effect in both treatments but those effects 
# have opposite signs. The same situation is true for testing treatment effects in probes with
# significant interaction. Therefore, we need to use the contrasts that tested genotype 
# differences within each level of treatment, and viceversa.
Genes.Int_Geno_I = cbind(Gene=with(out.annot, unique(gene_id[FDR.Int <= fdr_th & FDR.Geno_I <= fdr_th])), I=1)
Genes.Int_Geno_C = cbind(Gene=with(out.annot, unique(gene_id[FDR.Int <= fdr_th & FDR.Geno_C <= fdr_th])), C=1)
Genes.Int_Trt_B  = cbind(Gene=with(out.annot, unique(gene_id[FDR.Int <= fdr_th & FDR.Trt_B <= fdr_th])), B=1)
Genes.Int_Trt_BY = cbind(Gene=with(out.annot, unique(gene_id[FDR.Int <= fdr_th & FDR.Trt_BY <= fdr_th])), BY=1)

# Merge counts of genotype effects in both treatments into one table
Genes.Int_Geno   = merge(Genes.Int_Geno_I, Genes.Int_Geno_C, by="Gene", all=T)
Genes.Int_Geno   = cbind(Genes.Int_Geno$Gene,apply(Genes.Int_Geno[,2:3], 2, function(x) as.numeric(!is.na(x))))

# Merge counts of treatment effects in both genotypes into one table
Genes.Int_Trt    = merge(Genes.Int_Trt_B, Genes.Int_Trt_BY, by="Gene", all=T)
Genes.Int_Trt    = cbind(Genes.Int_Trt$Gene,apply(Genes.Int_Trt[,2:3], 2, function(x) as.numeric(!is.na(x))))

# Load the limma library for creating Venn diagrams
library(limma)

# Count genes for each compartment of the Venn diagram
Counts.Int_Geno=vennCounts(Genes.Int_Geno[,2:3])
print(Counts.Int_Geno)

Counts.Int_Trt=vennCounts(Genes.Int_Trt[,2:3]) 
print(Counts.Int_Trt)

# Plot genes responding to genotype in a treatment dependent manner
png(file.path(outdir, "vennDiagra_Int_Geno.png"), width=6, height=6, unit="in", res=150)
  vennDiagram(Counts.Int_Geno, main="\n\n\nGenes Responding to Genotype\nin a Treatment Dependent Manner")
dev.off()

# Plot genes responding to treatment in a genotype dependent manner
png(file.path(outdir, "vennDiagra_Int_Trt.png"), width=6, height=6, unit="in", res=150)
  vennDiagram(Counts.Int_Trt, main="\n\n\nGenes Responding to Treatment\nin a Genotype Dependent Manner")
dev.off()

# Interpretation of the Venn diagrams:
#
# In theory, both plots are showing the same test in two different ways, i.e. the number of
# genes with interaction effects, but partitioned either by treatment or by genotype.
# Because in practice we are showing results from four different (but related) tests, the
# total number of selected genes in each diagram is not exactly the same, but they should
# largely agree.
# 
# Although the numbers here are small because we used only a small sample or probes, you
# will see that more genes are responding to the treatment in the BY genotype. Also, you
# should see more differences between genotypes in the castrated animals. This was the
# pattern observed in the full dataset. See Figure 4 of Llamas et al 2008 (reference above).
# 
# Congratulations! You have completed this tutorial. Now you can use this script as a
# template and modify it for your own datasets. But first, I recommend you standing up and
# stretching a little. Data analyses can be addictive. Use with care.
