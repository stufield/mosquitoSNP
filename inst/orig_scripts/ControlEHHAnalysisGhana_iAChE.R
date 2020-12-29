###############################
### Stu Field & Martin Donnelly
### Sept 27, 2010
### Colorado State University
### sgf@colostate.edu
###############################
#####################
#####################
#########################
### Ghana iAChE data set
#########################
#####################
#####################
rm(list=ls()) #######
##########################################
require(mosquito)
##########################################
datafile <- "Phase input data minus F2a M-form.txt"
datafile2 <- "iACHEphase.txt"
path2file(datafile, setdir=TRUE)

#################################################
# Load & orgainize data set #####################
#################################################
data <- read.table(datafile, header=FALSE, skip=4, row.names=1) # core = "1085"
ColNames <- as.character(read.table(datafile, header=FALSE, skip=2, nrows=1, row.names=1))

data <- read.table(datafile2, header=FALSE, skip=4, row.names=1) # core = "1094"
ColNames <- as.character(read.table(datafile2, header=FALSE, skip=2, nrows=1, row.names=1)) 

colnames(data) <- ColNames
head(data)

SNPcore <- "1085"  # first data set
SNPcore <- "1094"  # second data set

######################
### Remove InDels ####
### using rm.col() ###
######################
new.data <- rm.col(data, index=0)

###########################
# Convert nucleotide codes
# using search.replace()
######################################
new.data <- search.replace(new.data, S= 1:4 , R= c("A","C","G","T"))
head(new.data)

###########################
# Separate into 2 datasets
# Resistant & Susceptible
# R = "A" * S = "G"
############################
R.data <- subset(new.data, subset= new.data[,SNPcore] == "A")
S.data <- subset(new.data, subset= new.data[,SNPcore] == "G")


################################
### lineage.map function call ##
################################
data.core <- which(names(R.data) == SNPcore)
lineage.map(R.data, core=data.core)
lineage.map(S.data, core=data.core)

### Save lineage map ###
S.data.map <- lineage.map(S.data, core=data.core, csv=FALSE, fileOut="Ghana iAChE_S_LineageMap_")$Haplo.Lineages
R.data.map <- lineage.map(R.data, core=data.core)$Haplo.Lineages

################################
### SNP.boot call ##############
################################
### First column MUST be character/factor; i.e. AA, pop#, etc.
S.data.map <- as.data.frame(cbind(rep(LETTERS[1],nrow(S.data.map)), S.data.map)) 
SNP.boot(S.data.map, nBoot=100, core= data.core, csv=FALSE, fileOut="Ghana iAChE_S_")

R.data.map <- as.data.frame(cbind(rep(LETTERS[1],nrow(R.data.map)), R.data.map)) 
SNP.boot(R.data.map, nBoot=25, core= data.core, csv=FALSE, fileOut="Ghana iAChE_R_")


