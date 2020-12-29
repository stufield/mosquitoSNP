###############################
### Haplotype Lineage for Amy Data
###############################
### Stu Field & Martin Donnelly
### Sept 27, 2010
### Colorado State University
### sgf@colostate.edu
###############################
rm(list=ls())
#############################################
require(mosquito)
###################################################
#path2file("HapData.csv", setdir=TRUE)
#data <- read.csv("HapData.csv") # test data set
#####################
#####################
#####################
#####################
## Ghana KDR Data ###
#####################
#####################
#####################
####################################
# Load & orgainize data set
########################################
datafile <- "AmyMformGhanaStu.txt"
path2file(datafile, setdir=TRUE)
#############################################
data <- read.table(datafile, header=FALSE)
names <- paste(data[,1], data[,2], sep="")
new.data <- data[-c(1,2)]
colnames(new.data) <- 1:ncol(new.data); rownames(new.data) <- names
new.data <- new.data[-59,]    # remove resistant individual
head(new.data)

###########################
# Convert nucleotide codes
# using search.replace()
######################################
new.data <- search.replace(new.data, S=1:4 , R=c("A","C","G","T"))

################################
### lineage.map() function call 
################################
lineage.map(new.data, core= 25, csv=FALSE)
Ghana.KDR <- lineage.map(new.data, core= 25, csv=FALSE)$Haplo.Lineages

#########################################
### SNP.boot function call ##############
#########################################
### First column MUST be character/factor; i.e. AA, pop#, etc.
Ghana.KDR <- as.data.frame(cbind(rep(LETTERS[1],nrow(Ghana.KDR)), Ghana.KDR)) 
SNP.boot(Ghana.KDR, nBoot=50, core= 25, csv=FALSE, fileOut="Ghana KDR_")


