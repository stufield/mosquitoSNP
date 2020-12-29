########################################################
LIMMA NORMALISATION OF RAW DATA  - DDT 8x15k MICROARRAYS
########################################################

#set wd
setwd ("c:/Documents and Settings/mjames/My Documents/SaraMicroarray/8x15K DDT MICROARRAY ANALYSIS")
getwd()
#open limma package
library(limma)
targets<- readTargets ("DDTTargets.txt", row.names="arraynumber")
targets
#weight spots i.e. remove spots which signal is below background using #boolean agilent column "rIsWellAboveBG" 1=yes 0=no, R&G sig=1 R+G>1
myfun<-function (x) {
      a=x[,"rIsWellAboveBG"] == 1
         b=x[,"gIsWellAboveBG"] == 1
           as.numeric(a+b >= 1)
}
#read in array files using weight function, take default Limma setting for #Agilent arrays i.e. foreground mean signal, background median signal used #in analysis. 
RG<- read.maimages(targets, source="agilent", wt.fun=myfun)
show(RG)
summary(RG)
summary(RG$R)
#check all genes and arrays read correctly, this will list number of genes #and number of arrays
dim(RG)
names(RG$genes)

#order genes so duplicates are next to eachother therefore order according #to 'RG$genes$ProbeName'
index=order(RG$genes$ProbeName)
RG=RG[index,]

#read in SpotType file providing info on different spots i.e. control type, #targets etc
spottypes<- readSpotTypes("DDTSpotTypes.txt")
spottypes
#match spot types to genes in the RG (read array)list ('status'=character #vector giving the control status of each spot on the array)
RG$genes$Status<- controlStatus(spottypes,RG)
#Create an MA-plot with color coding for control spots
plotMA(RG)
#or to create png images to file
plotMA3by2(RG)

#background correction half offset 50 recommended by limma
RG.b<- backgroundCorrect(RG, method="half", offset=50)
plotMA3by2(RG.b, prefix="rawintensitiesbkgdcorrectedHALF")

#SAVING PLOT IMAGE
#jpeg('name of plot.jpg')
#plotDensities(RG) or command to produce the plot
#dev.off()

#saving plot image of background corrected densities
jpeg('densitiesHALF.jpg')
plotDensities(RG.b)
dev.off()

#weighting different spots to be more/less important during normalisation #i.e. weight differentially expressed genes 0 and genes not meant to be #differentially expressed 1? “weighting simply multiplies the existing weights with either 1 #or 0. c(rep(1,4),rep(0,4))”

RG$weights=modifyWeights(RG$weights, RG$genes$Status, spottypes$SpotType, c (rep(1,3), rep(0,4)))
#within array normalisation method loess suggested
MA=normalizeWithinArrays(RG.b, method="loess")
 
 #between array normalisation Pie uses method scale? QUANTILE
 #MA<- normalizeBetweenArrays (RG.b, method="vsn")
 MA=normalizeBetweenArrays(MA, method="quantile")
 #below plot command produces blank plot???
 #plotMA3by2(MA, prefix="normalisedintensitiesQUANTILE",col=spottypes$color, #zero.weights=TRUE, common.lim=TRUE)
 #removed the col=spottype$color command and plots worked

 plotMA3by2(MA, prefix="normalisedintensitiesQUANTILE", zero.weights=TRUE, common.lim=TRUE)
 #pies script additional plots
 boxplot(MA$M[MA$weights==1]~col(MA$M)[MA$weights==1], names=colnames(MA$M), las=2, main="Normalised intensities QUANTILE between arrays", ylab="M", xaxt="none", ylim=c(-10,10))
 axis(1, at=c(1:length(targets[,1])), label=targets$Name, las=2)
 plotDensities(MA)

################################################
CONVERSION TO MAANOVA FORMAT AND RUNNING MAANOVA
################################################

#select detox only probes
#MA.detox=MA[MA$genes$Status=="detoxtarget",]
#index=order(MA.detox$genes$ProbeName)
#MA.detox=MA.detox[index,]
#MA.detox
RG.corrected<- RG.MA(MA)
MetaRow=1
MetaCol=1
#bind all the arrays together horizontally i.e. probe id running down #column 1 and each array red signal-green signal-weighting for probe #running horizontally for each of the 18 arrays..
ghanaDDT1.raw=as.data.frame(cbind(RG.corrected[,1]$genes[c(7:11, 4)], MetaRow, MetaCol,RG.corrected[,1]$genes[1:2], RG.corrected[,1]$R, RG.corrected[,1]$G, 1-RG.corrected[,1]$weights))
ghanaDDT1.raw.new=as.data.frame(cbind(RG.corrected[,2]$R, RG.corrected[,2]$G, 1-RG.corrected[,2]$weights))
ghanaDDT1.raw=cbind(ghanaDDT1.raw, ghanaDDT1.raw.new)
ghanaDDT1.raw.new=as.data.frame(cbind(RG.corrected[,3]$R, RG.corrected[,3]$G, 1-RG.corrected[,3]$weights))
ghanaDDT1.raw=cbind(ghanaDDT1.raw, ghanaDDT1.raw.new)
ghanaDDT1.raw.new=as.data.frame(cbind(RG.corrected[,4]$R, RG.corrected[,4]$G, 1-RG.corrected[,4]$weights))
ghanaDDT1.raw=cbind(ghanaDDT1.raw, ghanaDDT1.raw.new)
ghanaDDT1.raw.new=as.data.frame(cbind(RG.corrected[,5]$R, RG.corrected[,5]$G, 1-RG.corrected[,5]$weights))
ghanaDDT1.raw=cbind(ghanaDDT1.raw, ghanaDDT1.raw.new)
ghanaDDT1.raw.new=as.data.frame(cbind(RG.corrected[,6]$R, RG.corrected[,6]$G, 1-RG.corrected[,6]$weights))
ghanaDDT1.raw=cbind(ghanaDDT1.raw, ghanaDDT1.raw.new)
ghanaDDT1.raw.new=as.data.frame(cbind(RG.corrected[,7]$R, RG.corrected[,7]$G, 1-RG.corrected[,7]$weights))
ghanaDDT1.raw=cbind(ghanaDDT1.raw, ghanaDDT1.raw.new)
ghanaDDT1.raw.new=as.data.frame(cbind(RG.corrected[,8]$R, RG.corrected[,8]$G, 1-RG.corrected[,8]$weights))
ghanaDDT1.raw=cbind(ghanaDDT1.raw, ghanaDDT1.raw.new)
ghanaDDT1.raw.new=as.data.frame(cbind(RG.corrected[,9]$R, RG.corrected[,9]$G, 1-RG.corrected[,9]$weights))
ghanaDDT1.raw=cbind(ghanaDDT1.raw, ghanaDDT1.raw.new)
ghanaDDT1.raw.new=as.data.frame(cbind(RG.corrected[,10]$R, RG.corrected[,10]$G, 1-RG.corrected[,10]$weights))
ghanaDDT1.raw=cbind(ghanaDDT1.raw, ghanaDDT1.raw.new)
ghanaDDT1.raw.new=as.data.frame(cbind(RG.corrected[,11]$R, RG.corrected[,11]$G, 1-RG.corrected[,11]$weights))
ghanaDDT1.raw=cbind(ghanaDDT1.raw, ghanaDDT1.raw.new)
ghanaDDT1.raw.new=as.data.frame(cbind(RG.corrected[,12]$R, RG.corrected[,12]$G, 1-RG.corrected[,12]$weights))
ghanaDDT1.raw=cbind(ghanaDDT1.raw, ghanaDDT1.raw.new)
ghanaDDT1.raw.new=as.data.frame(cbind(RG.corrected[,13]$R, RG.corrected[,13]$G, 1-RG.corrected[,13]$weights))
ghanaDDT1.raw=cbind(ghanaDDT1.raw, ghanaDDT1.raw.new)
ghanaDDT1.raw.new=as.data.frame(cbind(RG.corrected[,14]$R, RG.corrected[,14]$G, 1-RG.corrected[,14]$weights))
ghanaDDT1.raw=cbind(ghanaDDT1.raw, ghanaDDT1.raw.new)
ghanaDDT1.raw.new=as.data.frame(cbind(RG.corrected[,15]$R, RG.corrected[,15]$G, 1-RG.corrected[,15]$weights))
ghanaDDT1.raw=cbind(ghanaDDT1.raw, ghanaDDT1.raw.new)
ghanaDDT1.raw.new=as.data.frame(cbind(RG.corrected[,16]$R, RG.corrected[,16]$G, 1-RG.corrected[,16]$weights))
ghanaDDT1.raw=cbind(ghanaDDT1.raw, ghanaDDT1.raw.new)
ghanaDDT1.raw.new=as.data.frame(cbind(RG.corrected[,17]$R, RG.corrected[,17]$G, 1-RG.corrected[,17]$weights))
ghanaDDT1.raw=cbind(ghanaDDT1.raw, ghanaDDT1.raw.new)
ghanaDDT1.raw.new=as.data.frame(cbind(RG.corrected[,18]$R, RG.corrected[,18]$G, 1-RG.corrected[,18]$weights))
ghanaDDT1.raw=cbind(ghanaDDT1.raw, ghanaDDT1.raw.new)
#add suffix to the array columns bound in the above command to indicate #which values are red/green and what the weighting for each spot is i.e. .R #.G .W at the end of array identifier
names(ghanaDDT1.raw)[seq(11,length(names(ghanaDDT1.raw)),3)]=paste(names(ghanaDDT1.raw)[seq(11,length(names(ghanaDDT1.raw)),3)],"R",sep=".")
names(ghanaDDT1.raw)[seq(12,length(names(ghanaDDT1.raw)),3)]=paste(names(ghanaDDT1.raw)[seq(12,length(names(ghanaDDT1.raw)),3)],"G",sep=".")
names(ghanaDDT1.raw)[seq(13,length(names(ghanaDDT1.raw)),3)]=paste(names(ghanaDDT1.raw)[seq(13,length(names(ghanaDDT1.raw)),3)],"W",sep=".")
ghanaDDT1.raw$ProbeName=substr(ghanaDDT1.raw$ProbeName, 1, 13)
#write a csv table of combined array info
write.csv(ghanaDDT1.raw, file="ghanaDDT1raw.csv")
#import array experimental design
design=read.table("DDTdesignfile.txt", header=TRUE)
write.table(ghanaDDT1.raw, "ghanaDDT.data.txt", row.names=FALSE, sep="\t", quote=FALSE)
library(maanova)
ghanaDDT<-read.madata(datafile="ghanaDDT.data.txt", designfile="DDTdesignfile.txt", arrayType="twoColor", header=TRUE, probeid=1, genename=2, systematicname=3, description=4, status=5, sequence=6, metarow=7, metacol=8, row=9, col=10, intensity=11, spotflag=TRUE, avgreps=0, log.trans=TRUE)
#log2 transformation of the data, 8x15 K = no reps
ghanaDDT
(ghanaDDT)
gridcheck(ghanaDDT)

#analysis according to Pie
#LOAD PACKAGE
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")

#FIT MIXED MODEL ANOVA (removed spot as no reps therefore no spot effect#
ghanaDDT.fit<-fitmaanova(ghanaDDT, formula=~Group+Array+Dye+Sample, random=~Array+Sample, covariate=~1)
save (ghanaDDT.fit, file="ghanaDDTfit")
#RESIDUAL PLOT
jpeg('DDTresidualplot.jpeg')
resiplot(ghanaDDT, ghanaDDT.fit)
dev.off()
#TEST STATISTIC (need to alter this command to use tabulated p-values)
test.unadj=matest(ghanaDDT, ghanaDDT.fit, test.type=c("ftest"), test.method=c(1,1), term="Group")
test=adjPval(test.unadj, method="jsFDR")
volcano(test, method="fdrperm")
#NOTE THE FDRPERM METHOD CAUSES TRUNCATION OF THE VOLCANO PLOT PERHAPS THE #PERMUTATION TEST IS NOT RUNNING CORRECTLY…?

# create a top table
ind=seq(1, length(ghanaDDT$probeid), 2)
DDT.tab=as.data.frame(cbind(ProbeID=test$obsAnova$probeid[ind], Gene=ghanaDDT$genename[ind], Description=ghanaDDT$description[ind], Status=ghanaDDT$status[ind]))
DDT.tab$logFC=test$obsAnova$Group[,1]-test$obsAnova$Group[,2]
DDT.tab$Pval=test$Fs$Ptab
DDT.tab$adjPval=test$Fs$adjPtab
DDT.tab$Flag=test$obsAnova$flag
write.table(DDT1.tab, row.names=F, file="ghanaDDT resultshalfquantile perm2.txt", sep="\t", quote=F)

# plot volcano and save to file
jpeg('volcano DDT1 perm 50')
plot(DDT1.tab$logFC, -log10(DDT1.tab$adjPval))

#add labels to volcano plot from values in results file i.e. foldchange #value and gene name pos= 1,2,3,4 1= below, 2=left, 3=above, 4=right
#i.e.
#plot(DDT1.tab$logFC, -log10(DDT1.tab$adjPval))
#text(x=1.343289, y=1.113238, labels="GSTS1_1", pos=2)
#text(x=1.670643, y=1.194637, labels="GSTS1_1", pos=2)
#text(x=1.784621, y=1.336079, labels="COEJHE2E", pos=2)
#text(x=1.202384, y=1.396278, labels="CYP325C2", pos=4)
dev.off()

save.image()
save.image("DDTanalysis250310.Rhistory")
savehistory()
savehistory("DDTanalysis250310.Rhistory")
