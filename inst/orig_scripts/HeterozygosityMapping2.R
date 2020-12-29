###############################
### Stu Field #################
### & #########################
### Martin Donnelly ###########
### Colorado State University #
### Dec 2, 2010 ##############
###### sgf@colostate.edu ######
##############################
require(pegas)
require(mosquito)
##########################
rm(list=ls()) ############
##########################
datafile <- "Harr Craig.csv"
##############################################################
path2file(datafile, setdir=TRUE)
###############################################################
#############################################
### Calculate Heterozygosities with Pegas ###
#############################################
###################################################
Data <- read.csv(datafile)
###############################################################
Locusnames <- c("L3RIH59", "L3RIND30", "L3RIR01", "L3RIR02", 
			    "L3RIR03", "L3RIR04", "L3RIR05", "L3RIR07", 
			    "L3RIR08", "L3RIR09", "L3RIR10", "L3RIR12")

### Convert ###
MicroData <- formatMicrosat(Data, col.loci=4:27, locus.names=Locusnames)
head(MicroData)

### Remove Alajo & Malawi because don't contain both phenotypes ####
MicroData <- subset(MicroData, Population!="Alajo" & Population!="Malawi")


### Convert to "loci" object for pegas functions ###
LociCols <- 4:15
MicroData <- as.loci(MicroData, row.names=1, col.pop=2, col.loci= LociCols)
MicroData



#################
### Quick Look
#################
print(MicroData, details=TRUE)
attr(MicroData, "locicol")
getAlleles(MicroData)
str(MicroData)

S <- summary(MicroData)
sapply(S, function(x) H(x$allele))        # Aggregate Heterozygosities (pegas)
plot(MicroData, las=2)


###############
### Fst #######
###############
Fst(MicroData)


##############################
# Calculate Heterozygosities #
# using HeFUN()
###########################################
Pops <- levels(MicroData$population)
types <- c("R","S")
L <- length(types)

### Summarize MicroData by Popn, Phenotype, & locus
### This is a list of each combination or Popn & Phenotype
### with each element of the list containing a list of each Locus
S <- by(MicroData, INDICES= c(MicroData["population"], MicroData["Phenotype"]), FUN= summary)
Ls <- length(S)

###############################
### Calculate He with HeFUN() ###
###############################
for (i in 1:Ls) {
   if (i==1) { He <- sapply(S[[i]], function(x) HeFUN(x$allele))   
   } else He <- cbind(He, He=sapply(S[[i]], function(x) HeFUN(x$allele)))
   ###################################
   # name columns by pop & phenotype #
   ###################################
   if (i==Ls) {
   	  A <- dimnames(S)
   	  Colnames <- NA; count <- 1
   	  Lj <- length(A$Phenotype)
   	  Lk <- length(A$population)
   	  for (j in 1:Lj) {
   	     for (k in 1:Lk) {
   	        Colnames[count] <- paste(A[[1]][k], A[[2]][j], sep="|")
   	        count <- count + 1
   	     }
   	  }
   colnames(He) <- Colnames
   He <- He[,order(colnames(He))]
   }
}
He
#read.csv("FStatHe.csv")

############################################
### Function to calculate 1 / Homozygosity ###
############################################
HoFUN <- function(x) { y <- (1 / (1 - x) )^2 - 1; y }

Ho <- HoFUN(He)
Ho

##################################
### Calculate log-ratio of R/S ###
##################################
index <- seq(1, ncol(Ho), by=2)
for (i in index) {
  logRatio <- log(Ho[,i] / Ho[,i+1])
  if (i==1) PlotData <- logRatio
  else PlotData <- cbind(PlotData, logRatio)
  if (i==last(index)) colnames(PlotData) <- Pops
}
PlotData


### Base pair positions ###
bp <- c(4282425,13629366,6773660,6816889,6834197,6887232, 
        6900390,7029097,7034734,7074313,7083144,7192866)

PlotData <- as.data.frame(cbind(PlotData, bp))
PlotData


### Reorder by "bp" for plotting ###
PlotData <- PlotData[order(PlotData["bp"]),]
PlotData

##################
### Plot Data ####
##################
plot(PlotData[,"bp"], PlotData[,1],
	ylim= c(min(PlotData[,-ncol(PlotData)]-0.5), ### remove final column (bp) in 
	       max(PlotData[,-ncol(PlotData)]+0.5)), ### determination of limits
	xlab= "Chromosome base pair",
	ylab= "log(Rh/Sh)",
	type= "b",
	pch= 19,
	cex= 0.75,
	lty= 1)
for (l in 2:ncol(PlotData)) {
	lines(PlotData[,"bp"], 
	      PlotData[,l], 
	      type="b", 
	      pch=19,
	      cex= 0.75, 
	      col=l)
}
abline(h=0, lty=2)


