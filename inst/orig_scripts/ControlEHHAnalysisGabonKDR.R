###############################
### Stu Field & Martin Donnelly
### Sept 27, 2010
### Colorado State University
### sgf@colostate.edu
###############################
########################
########################
########################
########################
### Gabon KDR Data Set #
########################
########################
########################
########################
rm(list=ls()) ##########
################################################
require(mosquito)
###################################################
datafile <- "AmyGabonJuly.txt"
path2file(datafile, setdir=TRUE)
###################################################
# Load & orgainize data set
####################################################
data <- read.table(datafile, header=FALSE)
names <- paste(data[,1], data[,2], sep="")
new.data <- data[-c(1,2)]
colnames(new.data) <- 1:ncol(new.data); rownames(new.data) <- names
head(new.data)


###########################
# Convert nucleotide codes
# using search.replace()
######################################
new.data <- search.replace(new.data, S= 1:4 , R= c("A","C","G","T"))
head(new.data)


###########################
# Separate into 3 datasets
############################
AA.index <- paste(new.data[,24], new.data[,25], sep="")
Leu.data <- subset(new.data, subset= AA.index == "TT")
Ser.data <- subset(new.data, subset= AA.index == "CT")
Phe.data <- subset(new.data, subset= AA.index == "TA")

################################
### lineage.map function call ##
################################
lineage.map(Leu.data, core= 25, csv=FALSE)   # Wildtype
lineage.map(Ser.data, core= 25, csv=FALSE)
lineage.map(Phe.data, core= 25, csv=FALSE)
Leu.boot <- lineage.map(Leu.data, core= 25, csv=FALSE)$Haplo.Lineages

################################
### SNP.boot function call #####
################################
### First col MUST be character/factor; i.e. AA, pop#, etc
Leu.boot <- as.data.frame(cbind(rep(LETTERS[1],nrow(Leu.boot)), Leu.boot))
SNP.boot(Leu.boot, nBoot=50, core= 25, csv=FALSE, fileOut="Gabon KDR_")

