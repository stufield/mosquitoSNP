##############################
### Field & Donnelly #########
### Oct. 8th 2010 ############
### Sliding Window ###########
### SNP Analysis #############
##############################
rm(list=ls()) ################
#####################################
### Set Dir & load required functions 
###############################################
require(mosquito)
###############################################
datafile <- "Sliding Window Data.csv"
path2file(datafile, setdir=TRUE)
######################################################################
mydata <- read.csv(datafile, header=TRUE, skip=4)
head(mydata)



####################
### Count data #####
### Conduct Fisher #
####################
mydata <- read.csv("Ghana MvS SNP data.csv", header=TRUE, skip=0)
mydata <- na.omit(mydata) # remove NAs
head(mydata)




###################################
### Calculate Fisher Exact test ###
###################################
p.value <- rep(NA,nrow(mydata))
for (i in 1:nrow(mydata)){
  p.value[i] <- fisher.test(matrix(as.numeric(mydata[i,4:7]), nrow=2, byrow=TRUE))$p.value
}

negLogP <- -log10(p.value)


###########################
### Define Data to plot ###
### Window in SNPs ########
###########################
Win <- 5
skip <- 1
negLogP <- mydata$negLogP
data.out <- slide.window(negLogP, window= Win, frame.skip= skip, FUN= max)

plot(data.out$values ~ data.out$midpoints, 
     type="s", 
     col= "navy", 
     ylab= "-logP (max)", 
     xlab= "midpoints", 
     lwd= 1.5)
legend("topleft", 
       legend= format(paste("WindowSize =",Win)), 
       bg= "gray75", 
       cex= 0.75, 
       box.lty= 0)






###########################
###########################
###########################
### Define Data to plot ###
### Window in bp ##########
###########################
###########################
###########################
BPwin <- 7500000
negLogP <- mydata$negLogP
data.out <- slide.window.bp(negLogP, Pos= mydata$ScaledPosition, bp= BPwin, FUN= max)



#####################
### Plot SNP data ###
#####################
chr.breaks <- c(60707344, 62812924, 110916950, 162910923, 206568235)
chr.cols <- c(2,2,3,2,3)
chr.lines <- c(2,2,1,2,1)
##############################
plotSlideSNP(x=data.out$midpoints, y=data.out$values, 
   type= "s", 
   col= "navy", 
   lwd= 1.25, 
   ylabel= "-logP (max)", 
   breaks= chr.breaks, 
   font= 4,
   window= Win,
   line.cols= chr.cols, 
   line.lty= chr.lines,
   save.plot= FALSE)

