########################
#     Pegas Example    #
#########################################
#         Stu Field                     #
#         Department of Biology         #
#         Colorado State University     #
#         Fort Collins, CO              #
#         80523-1878                    #
#         sgf@colostate.edu             #
#########################################
#         June 23, 2011                 #
#########################################
rm(list=ls())
##################
require(pegas)
#require(adegenet)
###############################################
#       Import as loci (pegas) object         #
###############################################
data <- read.loci("MicroSatData.csv", 
				loci.sep = ",", 
				row.names = 1,
				col.pop = 2, 
				col.loci = 3:11)
print(data, details=TRUE)
attributes(data)
data

###########################
#      Calculate Fst      #
###########################
Fst(data)


##########################################
#      Calculate Heterozygosities        #
#      using H() function (pegas)        #
##########################################
data(nancycats)    # import nancycats dataset
## convert the data and compute frequencies: 
S <- summary(as.loci(nancycats))
class(S)
S$allele[4]

## compute H for all loci: 
sapply(S, function(x) H(x$allele))

## ... include variance (row2)
sapply(S, function(x) H(x$allele, variance = TRUE))

?sapply

