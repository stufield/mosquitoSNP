###########################################################################
# LIMMA NORMALISATION OF RAW DATA_ 4x44K ARRAY LOOPED DESIGN
###########################################################################
require(limma)
require(maanova)
require(qvalue)
require(lattice)
require(qvalue)
require(calibrate)
require(lattice)
#####################
# set working directory – in R using LIMMA & MAANOVA 
setwd("~/Documents/Dropbox/CSU/Black & Donnelly/LIMMA") # <- set path here!

#read in target file
targets <- readTargets("Targets.txt", row.names="arraynumber")
targets

#weight spots i.e. remove spots which signal is below background using 
#boolean agilent column "rIsWellAboveBG" 1=yes 0=no, R&G sig=1 R+G>1
myfun <- function (x) {
  a=x[,"rIsWellAboveBG"] == 1
  b=x[,"gIsWellAboveBG"] == 1
  as.numeric(a+b >= 1)
}

#read in array files using weight function, take default Limma setting for 
#Agilent arrays i.e. foreground mean signal, background median signal used 
#in analysis. 

RG <- read.maimages(targets, source="agilent", wt.fun= myfun)
#show(RG)
summary(RG)

#summary(RG$R)
#check all genes and arrays read correctly, this will list number of genes 
#and number of arrays
#dim(RG)
#names(RG$genes)

#order genes so duplicates are next to each other
# therefore order according to 'RG$genes$ProbeName'
index <- order(RG$genes$ProbeName)
RG <- RG[index,]

#read in SpotType file providing info on different spots i.e. control type, 
#targets etc
spottypes <- readSpotTypes("SpotTypes.txt")
spottypes

#match spot types to genes in the RG (read array)list ('status'=character 
#vector giving the control status of each spot on the array)
RG$genes$Status <- controlStatus(spottypes,RG)

#Create an MA-plot with color coding for control spots
#plotMA(RG) or to create png images to file
plotMA3by2(RG, prefix="prenormalisationMA")

# prenorm dye density plot
plotDensities(RG) # or command to produce the plot

########################################################################
#background correction none offset 50 to increase low intensity values
#implements the background correction methods reviewed                   
#or developed in Ritchie et al (2007) and Silver at al (2009).  If 
#method="none" then no correction is done, i.e., the background intensities 
#are treated as zero. Other methods are all designed to produce positive 
#corrected intensities. If method="half" then any intensity which is less 
#than 0.5 after background subtraction is reset to be equal to 0.5. The 
#offset can be used to add a constant to the intensities before log-
#transforming, so that the log-ratios are shrunk towards zero at the lower 
#intensities. This may eliminate or reverse the usual 'fanning' of log-
#ratios at low intensities associated with local background subtraction. 
#Background correction (background subtraction) is also performed by the 
#normalizeWithinArrays method for RGList objects, so it is not necessary to 
#call backgroundCorrect directly unless one wants to use a method other 
#than simple subtraction. Calling backgroundCorrect before 
#normalizeWithinArrays will over-ride the default background correction
##########################################################################

RG.b <- backgroundCorrect(RG, method="none", offset=50)
plotMA(RG)
plotMA(RG.b)

# plotMA3by2(RG.b, prefix="backgroundnoneoffset50MA")
# plot image of background corrected densities
plotDensities(RG.b)


#weighting different spots to be more/less important during normalisation 
#i.e. weight differentially expressed genes 0 and genes not meant to be 
#differentially expressed 1? “weighting simply multiplies the existing weights with either 1 
#or 0. c(rep(1,4),rep(0,4))”

RG$weights <- modifyWeights(RG$weights, RG$genes$Status, spottypes$SpotType, c (rep(1,4), rep(0,4)))

#within array normalisation method loess suggested
#DEF: smoothing function LOESS in the R open-source software package, we 
#then fit a regression line to this data plot. LOESS is derived from the S 
#statistical function LOWESS1, which uses a locally weighted least squares 
#estimate of a regression fit. In simple terms, we take continuous sections 
#of the data and fit a linear (or quadratic) regression line to that area 
#of data. The function is then applied in continuity to the rest of the 
#data set, using a moving window of local data points to derive a fit line, 
#the result being a smoothed curvilinear regression line. The amount of 
#fitting and smoothing that takes place is governed by the span parameter 
#of the LOESS function which sets the proportion of the total data set to 
#be used in each window for local fitting. 
#http://www.obgyn.cam.ac.uk/genearray/loess-normalisation.htm
#Cleveland WS, Robust Locally Weighted Regression and Smoothing 
#Scatterplots, Journal of the American Statistical Association 
#1979;74:829-836
#Workman C, Jensen LJ et al. A new non-linear normalization method for 
#reducing variability in DNA microarray experiments. Genome Biol 
#2002;3(9):research0048

MA <- normalizeWithinArrays(RG.b, method="loess")

plotMA3by2(MA,prefix="bendio_bkgrdnone_offset50_loess_MA", col=spottypes$colour, zero.weights=TRUE, common.lim=TRUE)


#BOXPLOT – check the distribution of the red and green ratios for each 
#array to check they are similar in distribution; i.e. the median is 
#centred around zero and the distribution is of a similar size.  
#If arrays #are not similar then a BETWEEN array normalisation may be required…. 

boxplot(MA$M[MA$weights==1]~col(MA$M)[MA$weights==1], names=colnames(MA$M), las=2, main="bendio_bkgrdnone_offset50_loess_BOX", ylab="M", xaxt="none", ylim=c(-10,10))
axis(1, at=c(1:length(targets[,1])), label= targets$Name, las=2)


#PLOT DENSITIES – check that the red and green dye distribution are similar.  If they are not the dye correction within MAANOVA should correct for this but will need to check this later on...

plotDensities(MA)



#between array normalisation (remove the variations resulted from multiple 
#microarray experiments) Pie uses method scale? QUANTILE
#DEF: The scale normalization method was proposed by Yang et al (2001, 
#2002) and is further explained by Smyth and Speed (2003). The idea is 
#simply to scale the log-ratios to have the same median-absolute-deviation 
#(MAD) across arrays. This idea has also been implemented by the 
#maNormScale function in the marrayNorm package. The implementation here is 
#slightly different in that the MAD scale estimator is replaced with the 
#median-absolute-value and the A-values are normalized as well as the M-
#values. 
#Quantile normalization was proposed by Bolstad et al (2003) for 
#Affymetrix-style single-channel arrays and by Yang and Thorne (2003) for 
#two-color cDNA arrays. method="quantile" ensures that the intensities have 
#the same empirical distribution across arrays and across channels. 
#method="Aquantile" ensures that the A-values (average intensities) have 
#the same empirical distribution across arrays leaving the M-values (log-
#ratios) unchanged. These two methods are called "q" and "Aq" respectively 
#in Yang and Thorne (2003). method="Tquantile" performs quantile 
#normalization separately for the groups indicated by targets. targets may 
#be a target matrix such as read by readTargets or can be a vector 
#indicating green channel groups followed by red channel groups. 

#NO between array norm as recommended by Agilent and Zahurak et al. 2007





###########################################################################
# CONVERSION TO MAANOVA FORMAT AND RUNNING MAANOVA ########################
###########################################################################
RG.corrected <- RG.MA(MA)
setwd("~/Documents/Dropbox/CSU/Black & Donnelly/LIMMA")


###########################################################################
source(paste(getwd(),"array.frame().R",sep="/"))  ### <- set path here
###########################################################################
#### Call array.frame() ###############
#######################################
### Determine ColNames to be included #
colNames <- c("ProbeName","GeneName","SystematicName","Description", "Status","Sequence","Row","Col")
Tororo.raw <- array.frame(RG.corrected, names= colNames)


###################
### Save Result ###
###################
write.csv(Tororo.raw, file="Tororo.csv")  # not needed for maanova

### Can be exported and then imported in read.madata() below or 
### can use Tororo.raw object directly in read.madata()
write.table(Tororo.raw, file="Tororo.txt", row.names=FALSE, sep="\t", quote=FALSE)

# after here data files can go either to Maanova or J/Maanova


#####################
#####################
#####################
#####################
### MAANOVA Analysis
#####################
#####################
#####################
#####################
require(maanova)
####################

# import array experimental design
design <- read.table("bendio1designfile.txt", header=TRUE)

write.table(ghanabendio1.raw, "ghanabendio1_bkgrd_none_offset50_loess.data.txt", row.names=FALSE, sep="\t", quote=FALSE)

#import data.txt and design file j/MAANOVA or continue in R and MAANOVA
#Note – can remove control spots at this point manually from the data file 
#are they are no longer required in analysis...
#Note – ensure order of arrays and dyes in the design file EXACTLY match the ordering in the data file
#
#JMAANOVA COMMANDS:
#read in data and Log transform (as recommended).  Specify no. of reps if 
#want to collapse them(4x44K has n rep=2) and where each bit of info is 
#stored in the data file as below.
#
#FIT MIXED EFFECT MODEL: (restricted maximum likelihood method used to fit 
#linear model as default setting) Given the data and formula, this function 
#fits the regression model for each gene and calculates the ANOVA 
#estimates, variance components for random terms, fitted values, etc. For a 
#mixed effect models, the output estimates will be BLUE and BLUP.



# The function to import data into J/maanova was:
##########################
Tororo.raw <- read.table("Tororo.txt", header=TRUE, sep="\t") # shortcut to above
##########################
Tororo <- read.madata(
  datafile= as.matrix(Tororo.raw),          # must be a matrix -> convert
  #datafile= "~/Desktop/LIMMA/Tororo.txt",  # or from imported from *.txt
  designfile= "~/Desktop/LIMMA/Designfile.txt", 
  arrayType="twoColor",
  header=TRUE, 
  spotflag=TRUE, 
  n.rep=1, 
  probeid=1,
  GeneName=2,
  Description=4,
  Status=5,
  row=7,
  col=8,
  metarow=9,
  metacol=10,
  intensity=11, 
  log.trans=TRUE)

names(Tororo); class(Tororo)



# Note: the change in metarow and metacol positions relative to previous versions.
# No prelication of probes in 8x15K array (replication within genes is present) 
# The script to fit a model was.....
Tororo.Fit <- fitmaanova(
   madata= Tororo, 
   formula= ~Array+Dye+Sample+Group, 
   random= ~Array+Sample, 
   method= "REML", 
   verbose= TRUE, 
   subCol= FALSE)

names(Tororo.Fit)


######################################################################
# Model test was based on a single permutation without replication
# F-TEST UNIVERSAL – tests all combinations of groups. As the experimental 
# size is small (<5 pools per group) permutation testing in not viable or 
# recommended therefore we will rely on tabulated pvalues hence n.perm=1.
############################################################################
Tororo.Fit.F <- matest(
   data= Tororo, 
   anovaobj= Tororo.Fit,
   term= "Group", 
   test.type= "ftest", 
   MME.method= "REML", 
   n.perm= 1, 
   verbose=TRUE)

names(Tororo.Fit.F); class(Tororo.Fit.F)
volcano(Tororo.Fit.F)

FoldChange <- calVolcanoXval(Tororo.Fit.F) ### vector 1 col

###########################
# pairwise T test: defined
###########################
Tororo.Fit.T <- matest(
   data= Tororo, 
   anovaobj= Tororo.Fit, 
   term= "Group", 
   test.type= "ttest", 
   Contrast= new("matrix", nrow=3, ncol=3, byrow=TRUE, 
                           #   K     R    T
                           c(-1.0,  1.0, 0.0,   ### Comp #1
                              0.0, -1.0, 1.0,   ### Comp #2
                             -1.0,  0.0, 1.0)), ### Comp #3
   MME.method= "REML", 
   n.perm= 1, 
   verbose= FALSE)

names(Tororo.Fit.FT; class(Tororo.Fit.T)
volcano(Tororo.Fit.T, highlight=FALSE)

FoldChange <- calVolcanoXval(Tororo.Fit.T)   # matrix 3 cols



########################
####### Shortcuts ######
load("Tororo.Fit.Rdata")    # Tororo.Fit
load("Tororo.Fit.T.Rdata")  # Tororo.Fit.T
load("Tororo.Fit.F.Rdata")  # Tororo.Fit.F
#####################################################
########################################################################
# CREATE A RESULTS TABLE from pairwise or F test as there are only two 
# treatment groups here -- create one before and after FDR. 
# ENTER EACH LINE ONE AT A TIME!
######################################################################
index <- 1:length(Tororo$probeid)
Tororo.tab <- data.frame(
	ProbeID = Tororo.Fit.T$obsAnova$probeid[index], 
	Status = Tororo$Status[index],
	Gene = Tororo$GeneName[index], 
	Description = Tororo$Description[index],
	logFC = FoldChangeT,
	Pval = Tororo.Fit.T$Fs$Ptab,
	Flags = Tororo.Fit.T$obsAnova$flag)

head(Tororo.tab)
dim(Tororo.tab)


###########################################################################
#ADJUST P-VALUE - Benjamini and Hochberg False discovery rate (FDR. The 
#FDR is the expected false positive rate; for example, if 1000 observations 
#were experimentally predicted to be different, and a maximum FDR for these 
#observations was 0.10, then 100 of these observations would be expected to 
#be false positives. The q value is defined to be the FDR analogue of the 
#p-value. The q-value of an individual hypothesis test is the minimum FDR 
#at which the test may be called significant. One approach is to directly 
#estimate q-values rather than fixing a level at which to control the FDR.
#load required packages for p-value adjustment
##################################################
Qobj <- qvalue(as.matrix(Tororo.tab[,c("Pval.1","Pval.2","Pval.3")]))
names(Qobj)
qplot(Qobj)


logQvals <- -log10(Qobj$qvalues)   # pull out qvalues and log
colnames(logQvals) <- c("logQval.1","logQval.2","logQval.3")  # rename

head(logQvals)

Tororo.tab <- cbind(Tororo.tab, logQvals) # bind qvalues onto data frame
head(Tororo.tab)                          # take a look; 
dim(Tororo.tab)                           # make sure added


####################################################
### Subset #########################################
### Remove: CV, negativecontrol, & positivecontrol #
####################################################
Tororo.tab <- subset(Tororo.tab, Tororo.tab["Status"]=="detoxtarget" | Tororo.tab["Status"]=="WGtarget")


##################
### Quick Plots ##
##################
par(mfrow=c(1,3))
plot(logQval.1 ~ logFC.1, data=Tororo.tab, cex=0.5, col="navy", main="Comparison #1")
plot(logQval.2 ~ logFC.2, data=Tororo.tab, cex=0.5, col="darkgreen", main="Comparison #1")
plot(logQval.3 ~ logFC.3, data=Tororo.tab, cex=0.5, col="darkred", main="Comparison #1")


###############################################
### Subset ###########
### By critical Qval
### and Fold Change
###############################################
critP = -log10(0.0001)     # set critical P value = 0.0001

### Subset based on Pvals of comparison #2 (group 3 - group 2)
Left <- subset(Tororo.tab, logQval.2 > critP & logFC.2 < 0)
Right <- subset(Tororo.tab, logQval.2 > critP & logFC.2 > 0)
Left
Right


###########################################################
### Plot highlighting top Pvals
###########################################################
par(mfrow=c(1,3))
### Compare Sig. points of focal comparison (#2) in comparison #1
plot(logQval.1 ~ logFC.1, data=Tororo.tab, cex=0.5, col="navy")
points(logQval.1 ~ logFC.1, data=Left, col=2, pch=19, cex=0.5)
points(logQval.1 ~ logFC.1, data=Right, col=3, pch=19, cex=0.5)

### Compare Sig. points of focal comparison (#2) in comparison #2
plot(logQval.2 ~ logFC.2, data=Tororo.tab, cex=0.5, col="navy")
points(logQval.2 ~ logFC.2, data=Left, col=2, pch=19, cex=0.5)
points(logQval.2 ~ logFC.2, data=Right, col=3, pch=19, cex=0.5)
title("Comparison #2 on ...")

### Compare Sig. points of focal comparison (#2) in comparison #3
plot(logQval.3 ~ logFC.3, data=Tororo.tab, cex=0.5, col="navy")
points(logQval.3 ~ logFC.3, data=Left, col=2, pch=19, cex=0.5)
points(logQval.3 ~ logFC.3, data=Right, col=3, pch=19, cex=0.5)



#################################
### Labelling ###################
### Subset points of interest ##########
### On Comparison/Graph #3 = K vs T ####
########################################
s <- c("CUST_7725_PI4","CUST_8478_PI4","CUST_12524_PI")   # search for
r <- c("stu1","stu2","stu3")                              # replace with

Tororo.tabX <- searchReplace(Tororo.tab, S=s, R=r)
head(Tororo.tabX)

########################################
### Label points by clicking on them ###
### Labels are those in Tororo.tabX ####
########################################
identify(Tororo.tab$logFC.3, Tororo.tab$logQval.3, labels=Tororo.tabX$ProbeID, cex=0.5) ### click on points for ProbeID




####################################
### Heatmap from a Matrix ##########
####################################
x <- matrix(runif(225,1,20),ncol=15)  ### dummy matrix data
colnames(x) = as.character(1:ncol(x))
rownames(x) = as.character(1:nrow(x))

stu.heatmap <- function(x, key="right", grayscale=FALSE, ...){
   grid <- cbind(expand.grid(Rows=rownames(x), Cols=colnames(x)), x= as.vector(x))
   grid$Rows <- ordered(grid$Row, levels=rev(levels(grid$Cols))) ## reorder axis
   n <- range(grid$x); L <- length(x)
   if (grayscale) type <- rev(grey(1:L/L)) else type <- rev(heat.colors(L))
   levelplot(x ~ Cols*Rows, 
      data= grid, 
      at= seq(n[1], n[2], length=L), 
      col.regions= type,
      colorkey= list(space=key), 
      pretty= TRUE, ...)
}

stu.heatmap(x, grayscale=FALSE, main="Matrix Heatmap")



##################
### Plotting junk
##################
require(rgl)
require(emdbook)
z <- rnorm(15000) + atan2(x,y)
plot3d(x, y, z, col=rainbow(1000))
plot3d(Tororo.tab$logFC, -log10(Tororo.tab$adjPval[,1]), z, col=rainbow(15000), type="s", size=0.3)

