###############################
### Field & Donnelly ##########
### Sept 15, 2010 #############
### Bootstrap for mosquito SNP data
#####################################
data <- read.csv("Kenyan SNP EHH Data.csv", header=TRUE)    # read into R memory
head(data)

#   Subset Data sets by Allele   #
data.Leu <- subset(data, data == "Leu")
data.Ser <- subset(data, data == "Ser")
SNP.boot(data.Leu, 50, core=21)
SNP.boot(data.Ser, 25, core=21)

