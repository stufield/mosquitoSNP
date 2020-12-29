############################
### Translating a ##########
### diploid DNA sequence ###
### into an AA sequence ####
############################
### Stu Field & Bill Black #
### Colorado State Univ ####
### sgf@colostate.edu ######
### Sept 30th, 2010 ########
############################
### function to calculate
### possible AAs from
### sequences containing
### ambiguity codons
############################
rm(list=ls()) ##############
###############################################
### Set Dir & load required package
###############################################
require(mosquito)
datafile <- "CodonTable.csv"
path2file(datafile, setdir=TRUE)
################################################
### Sample DNA sequence ####
############################
DNAseq <- c("YYYAAGATCGTCGGTGGCGATGAGGCCGAAGCGCACGAATTTCCCTACCAAATCTCGCTGCAGTGGAACTTCAACGATGGACAAACGGAGACCATGCACTTCTGYGGAGCTTCGGTGTTGAACGAAAACTTYGTCCTGACGGCTGCTCACTGCAAGACCGCATACTCCAATACCGGGTWCATCGAAGTGGTTGCCGCTGAACATGATGTGGCYGTTGCGGAAGGATCCGAACAGCGTCGYTTGGTTGCGGAGTTCATCGTCCACGAGGACTATCAAGGRGGAGTCAGTCCCGATGAGATTGCCGTCA")


###################
### Function call
###################
AAtranslate("YRY")
DNAtranslate(DNAseq, csv=FALSE)

