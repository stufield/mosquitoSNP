###############################
### Stu Field #################
### Martin Donnelly ###########
### Colorado State University #
### Nov 23, 2010 ##############
###############################
datafile <- "HarrCraigData.csv"
path2file(datafile, setdir=TRUE)
data <- read.csv(datafile)
######################################################
lociCols <- 4:27                                     # define Loci columns
lociBind <- seq(lociCols[1], tail(lociCols,1), 2)    # Define Loci to bind
nLoci    <- length(lociBind)                         # No. of Loci


Locusnames <- c("L3RIH59", "L3RIND30", "L3RIR01",
                "L3RIR02", "L3RIR03", "L3RIR04", 
                "L3RIR05", "L3RIR07", "L3RIR08", 
                "L3RIR09", "L3RIR10", "L3RIR12")

### If you want to convert format to allele1/allele2 use formatMicrosat()
data2 <- formatMicrosat(data, col.loci=4:27, locus.names=Locusnames)
head(data2)

#################################
# split up again if necessary
# should be same as original data
#################################
#for ( i in 1:12 ) {
#   split <- strsplit(as.character(data[, i + 3]), "/")
#   a <- as.numeric(sapply(split, "[", 1))
#   b <- as.numeric(sapply(split, "[", 2))
#   Y <- cbind(a, b)
#   if ( i==1 ) {
#      X <- Y
#   } else {
#      X <- cbind(X, Y)
#   }
#}
#splitdata <- cbind(data[,1:3],X)
#colnames(splitdata) <- colnames(data)
#head(splitdata)





########################################################
# Calculate Means per locus per phenotype per popn ###
########################################################
for ( i in lociBind ) {
  tableA <- tapply(data[,i],
                   INDEX=c(data["Phenotype"], data["Population"]),
                   FUN=mean, na.rm=TRUE)
  tableB <- tapply(data[,i+1],
                   INDEX=c(data["Phenotype"], data["Population"]),
                   FUN=mean, na.rm=TRUE)
  table <- (tableA + tableB) / 2
  if ( i==lociBind[1] ) {
     Means <- list(table)
     counter <- 1
  } else {
     Means[[counter]] <- table
  } 
  counter <- counter + 1
}
names(Means) <- Locusnames
Means
print(length(Means)==nLoci)    # double-check



###########################################################
### Calculate Variance per locus per phenotype per popn ###
###########################################################
for ( i in lociBind ) {
  tableA <- tapply(data[,i],
                   INDEX=c(data["Phenotype"], data["Population"]),
                   FUN=var, na.rm=TRUE)
  tableB <- tapply(data[,i+1],
                   INDEX=c(data["Phenotype"], data["Population"]), 
                   FUN=var, na.rm=TRUE)
  table <- (tableA + tableB) / 2
  if ( i==lociBind[1] ) {
     Sigma <- list(table)
     counter <- 1
  } else {
      Sigma[[counter]] <- table
  }
  counter <- counter + 1
}
names(Sigma) <- Locusnames
Sigma


####################################################
### Calculate N per locus per phenotype per popn ###
####################################################
for ( i in lociBind ) {
  y <- data[ !is.na(data[,i]), ]    # remove NAs -> na.omit()???
  tableA <- tapply(y[,i], INDEX=c(y["Phenotype"], y["Population"]), FUN=length)
  tableB <- tapply(y[,i+1], INDEX=c(y["Phenotype"], y["Population"]), FUN=length)
  table  <- (tableA + tableB) / 2
  if ( i==lociBind[1] ) {
     N <- list(table)
     counter <- 1
  } else {
      N[[counter]] <- table
  }
  counter <- counter + 1
}
names(N) <- Locusnames
N


#####################################################################
# Calculate the 'unbiased variance' per locus per phenotype per popn
#####################################################################
for ( i in 1:nLoci ) {
  Vartable <- (Sigma[[i]] * N[[i]]) / (Means[[i]] * (N[[i]]-1))
  if ( i==1 ) {
     CorrectedVar <- list(Vartable)
  } else {
     CorrectedVar[[i]] <- Vartable
  }
}
names(CorrectedVar) <- Locusnames
CorrectedVar
print(length(CorrectedVar)==nLoci)   # double check



##############################################################
# Calculate the MEAN 'unbiased variance' per locus per phenotype
##############################################################
VR <- VS <- NA
for ( j in 1:nLoci ) {
  VR[j] <- mean(CorrectedVar[[j]]["R",], na.rm=TRUE)
  VS[j] <- mean(CorrectedVar[[j]]["S",], na.rm=TRUE)
}
names(VR) <- names(VS) <- Locusnames
VR
VS

print(length(VR)==nLoci)




###############################
# Calculate log-Variance Ratio
###############################
logVarRS <- log(VR/VS)
logVarRS



#############
# Plotting
#############
plot(1:nLoci, logVarRS, type="b",
     main="Selective Sweep Trough",
     xlab="Nucleotide (bp)", 
     ylab="log-Variance Ratio (VR / VS)",
     pch=19, col="navy", cex=0.75, lwd=1.5)
abline(h=0, lty=2)
box()




###################################
### SAS code                      #
### Estimating the 'trough'       #
### Needs mutation rate estimates #
###################################
/*data het2; set het1;

/*proc sort; by lnRS;
run;
data xi; set x;
if locus='2'  then vR1l=VR; 
proc print; 
run;
data xii; set x;
if locus='5'  then vS1l=VR; 
proc print; 
run;               *<Need to work out how to get the program to select the top twop values without inputting the locus numbers manually;
if locus='5' then vs2l=lnrs;
X=((-4*(8.3*10**-6)*(LOG((-VS1L)/(VR1L-VS1L))))-(k*n*(1.2*10**-8)*(LOG((-VS1L)/(VR1L-VS1L))))+(4*(8.3*10**-6)*(LOG((-VS2L)/(VR2L-VS2L)))))/((n*(1.2*10**-8)*(LOG((-VS1L)/(VR1L-VS1L))))-(n*(1.2*10**-8)*(LOG((-VS2L)/(VR2L-VS2L)))));
* Eq 1 Harr et al should calculate the distance between the microsatellite;
Proc print; run;

/* The following equation Eq 1 Harr et al should calculate the distance between the microsatellite and 

X=((-4*(8.3*10**-6)*(LOG((-VS1L)/(*VR1L-VS1L))))-(k*n*(1.2*10**-8)*(LOG((-VS1L)/(*VR1L-VS1L))))+(4*(8.3*10**-6)*(LOG((-VS2L)/(*VR2L-VS2L)))))/((n*(1.2*10**-8)*(LOG((-VS1L)/(*VR1L-VS1L))))-(n*(1.2*10**-8)*(LOG((-VS2L)/(*VR2L-VS2L)))));  *by phenotype locus ascending;
*run;
*DATA Distance ;                   
*Phenotype='AAAAAAAA'; *Locus=111111111111; *V=44444; 

**************************************************************
*DATA a ; *SET final ;
*IF FIRST.locus;                 *< takes the first locus for each population from above ordering;
*   proc print ;
*   PROC APPEND BASE=BOOT; * must think of a way of selecting the same locus in
   the susceptible populaiton 
run;
*
********************88

