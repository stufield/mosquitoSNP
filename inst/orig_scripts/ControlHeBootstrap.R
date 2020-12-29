###############################
### Stu Field #################
### & #########################
### Martin Donnelly ###########
### Colorado State University #
### Dec 2, 2010 ##############
###### sgf@colostate.edu ######
##############################
require(pegas)
require(gplots)
require(mosquito)
###########################################################
datafile <- "Harr Craig.csv"
path2file(datafile, setdir=TRUE)
#########################################
### mini order function for use with lapply
order2 <- function(x,n) { 
					if( !is.data.frame(x) ) 
						x <- as.data.frame(x)
					x[ order(x[n]), ]
}

#############################################
### Import Raw Data & manipulate ############
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

###############################################################
### Convert to "loci" object for pegas functions ###
###############################################################
MicroData <- as.loci(MicroData, row.names=1, col.pop=2, col.loci= 4:15)
MicroData


#################
### Quick Look ##
#################
print(MicroData, details=TRUE)
attr(MicroData, "locicol")
getAlleles(MicroData)
str(MicroData)
S <- summary(MicroData)
sapply(S, function(x) H(x$allele, variance=TRUE))   # Aggregate Heterozygosities
plot(MicroData, las=2)

###############
### Fst #######
###############
Fst(MicroData)


##############################
# Calculate Heterozygosities
#################################
# function for Heterozygosities #
# Corrected from pegas function #
###########################################
### Summarize MicroData by Popn, Phenotype, & locus
### This is a list of each combination or Popn & Phenotype
### with each element of the list containing a list of each Locus
S <- by(MicroData, INDICES=c(MicroData["population"], MicroData["Phenotype"]), FUN=summary)
Ls <- 1:length(S)



###############################
### Calculate He with HeFUN()
### Using variance calculated in HeFUN()
### NOT Bootstrap method
######################################
source(paste(FD,"HeFUN().R", sep="/"))
######################################
for (i in seq(length(S))) {
   X <- sapply(S[[i]], function(x) HeFUN(x$allele, var=TRUE))
   if (i==1) { 
		He <- X[1,]
		HeVar <- X[2,]
		He.n <- X[3,]
	}
   else {
		He <- cbind(He, He=X[1,])
		HeVar <- cbind(HeVar, HeVar=X[2,]) 
		He.n <- cbind(He.n, He.n=X[3,]) 
	}
	###################################
	# name columns by pop & phenotype #
	###################################
	if ( i==length(S) ) {
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
		colnames(He) = colnames(HeVar) = colnames(He.n) <- Colnames
		He <- He[,order(colnames(He))]
		HeVar <- HeVar[,order(colnames(HeVar))]
		He.n <- He.n[,order(colnames(He.n))]
	}
	HeSE <- 1.96 * (sqrt(HeVar) / sqrt(He.n))
	HeList <- list(He=He, variance=HeVar, n=He.n, HeCI95=HeSE)
}

## Errors will be symmetrical ##
HeList


#####################
####### Plotting ####
#####################
bp <- c(4282425,13629366,6773660,6816889,6834197,6887232, 
        6900390,7029097,7034734,7074313,7083144,7192866)

### Add bp to each list entry data frame ###
HeList <- lapply(HeList, cbind, bp)

### Reorder by "bp" and by list entry ###
HeList <- lapply(HeList, order2, "bp")

### Create x-axis ###
xaxis <- HeList$He[,"bp"]

### Resistants ###
PlotPop <- 1           # used to index the population to plot; 1 = 1st/2nd col
### To plot all 4 pops, simply put a loop around the two plotting functions 
### below & index PlotPop = i, where i becomes: for (i in seq(1,7,2)) { etc.
plotCI(x= xaxis/1000, 
		y= HeList$He[,PlotPop],
		uiw= HeList$HeCI95[,PlotPop],
		col= "navy",
		barcol= "navy",
		pch= 19, 
		gap= 0,
		type= "b",
		cex= 0.5,
		xlab= "Base pair along chromosome 3R (x1000)",
		ylab= "Expected Heterozygosity"); box()

legend("right", legend=c("R","S"), col=c("navy","darkred"), lty=1, pch=19, bg="gray80", pt.cex=0.7, box.lty=0)

### Susceptibles ###
plotCI(x= xaxis/1000, 
		y= HeList$He[,PlotPop+1],
		uiw= HeList$HeCI95[,PlotPop+1],
		col = "darkred",
		barcol= "darkred",
		pch= 19,
		gap= 0,
		type= "b",
		cex= 0.5,
		xlab= "",
		add= TRUE)













#######################
### Error CI95 of He ##
### via Bootstrapping ##
###########################
### Set Working Directory #################################
setwd("~/Documents/Dropbox/CSU/Black & Donnelly/R code/Data")
Data <- read.csv("Harr Craig.csv")
###############################################################
Locusnames <- c("L3RIH59", "L3RIND30", "L3RIR01", "L3RIR02", "L3RIR03", "L3RIR04",
                "L3RIR05", "L3RIR07", "L3RIR08", "L3RIR09", "L3RIR10", "L3RIR12")

MicroData <- formatMicrosat(Data, col.loci= 4:27, locus.names= Locusnames)
MicroData <- subset(MicroData, Population != "Alajo" & Population != "Malawi")
head(MicroData)

PlotData <- HeBootstrapFUN(MicroData, nboot=5, LociCols=4:15)
PlotData

# load(file="PlotDataBoot1000.rda") # shortcut








#####################
####### Plotting ####
#####################
bp <- c(4282425,13629366,6773660,6816889,6834197,6887232, 
        6900390,7029097,7034734,7074313,7083144,7192866)

PlotData <- lapply(PlotData, cbind, bp)

### Reorder by "bp" and by list entry ###
PlotData <- lapply(PlotData, order2, "bp")

xaxis <- PlotData$estimate[,"bp"]

### Resistants ###
PlotPop <- 1           # used to index the population to plot; 1 = 1st/2nd col
plotCI(x= xaxis/1000,
	   y= PlotData$estimate[,PlotPop], 
	   li= PlotData$lowerCI95[,PlotPop], 
	   ui= PlotData$upperCI95[,PlotPop],
	   col= "navy",
	   barcol= "navy",
	   pch= 19, 
	   gap= 0,
	   type= "b",
	   cex= 0.7,
	   xlab= "Base pair along chromosome 3R (x1000)",
	   ylab= "Expected Heterozygosity"); box()

legend("right", legend=c("R","S"), col=c("navy","darkred"), lty=1, pch=19, bg="gray80", pt.cex=0.7, box.lty=0)

### Susceptibles ###
plotCI(x= xaxis/1000,
	   y= PlotData$estimate[,PlotPop+1], 
	   li= PlotData$lowerCI95[,PlotPop+1], 
	   ui= PlotData$upperCI95[,PlotPop+1],
	   col = "darkred",
	   barcol= "darkred",
	   pch= 19,
	   gap= 0,
	   type= "b",
	   cex= 0.7,
	   xlab= "",
	   add= TRUE)


