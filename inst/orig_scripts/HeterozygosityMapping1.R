###############################
### Stu Field #################
### & #########################
### Martin Donnelly ###########
### Colorado State University #
### Nov 29, 2010 ##############
###### sgf@colostate.edu ######
### Heterozygosities Calculated
### via External Program ######
##############################
rm(list=ls()) ##############
##############################
datafile <- "Log RH Craig.csv"
##################################################################
path2file(datafile, setdir=TRUE)
##################################################################
data <- read.csv(datafile)
head(data)


Locusnames <- c("L3RIH59", "L3RIND30", "L3RIR01", "L3RIR02", 
			    "L3RIR03", "L3RIR04", "L3RIR05", "L3RIR07", 
			    "L3RIR08", "L3RIR09", "L3RIR10", "L3RIR12")
L <- length(Locusnames)
Pops <- levels(data$Population)


### Function to calculate 1 / Homozygosity ###
HomoFUN <- function(x) {
   y <- (1 / (1 - x) )^2 - 1
   y
}


### add measure to data frame ###
data$Homozygosity <- HomoFUN(data[,"Heterozygosity"])

### Subset the Data by pop & phenotype
Y <- by(data, INDICES=c(data["Phenotype"], data["Population"]), FUN=list)

##################################
### Calculate log-ratio of R/S ###
##################################
index <- seq(1, length(Y), by=2)
for (i in index) {
  logRatio <- log(Y[[i]][,"Homozygosity"] / Y[[i+1]][,"Homozygosity"])
  if (i==1) PlotData <- logRatio
  else PlotData <- cbind(PlotData, logRatio)
}

colnames(PlotData) <- Pops
bp <- c(4282425,13629366, 6773660,6816889, 
        6834197,6887232, 6900390,7029097, 
        7034734,7074313, 7083144,7192866)

PlotData <- as.data.frame(cbind(PlotData,bp), row.names= as.character(Y[[1]][,"Locus"]))
PlotData


### Reorder by "bp" for plotting ###
PlotData <- PlotData[order(PlotData["bp"]),]


##################
### Plot Data ####
##################
plot(PlotData[,"bp"], PlotData[,1], 
	xlab= "Chromosome base pair", 
	ylab= "log(Rh/Sh)",
	type= "b", 
	pch= 19, 
	lty= 1)
for (l in 2:ncol(PlotData)) {
	lines(PlotData[,"bp"], PlotData[,l], type="b", pch=19, col=l)
}
abline(h=0, lty=2)

