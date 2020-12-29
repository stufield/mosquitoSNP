###############################################
# Universal translation of all diploid codons
###############################################
codon.p <- c("A","C","G","T","R","M","W","S","Y","K")
L <- 1:length(codon.p)
codons <- NA
count  <- 0
clist  <- list()

for(a in L) {
  for(b in L) {
    for(c in L) {
      count <- count + 1
      codons[count] <- paste0(codon.p[a], codon.p[b], codon.p[c])
    }
  }
  cat(count,"\n")
  clist[[a]] <- codons[(count-(b*c-1)):count]
}
codon.p
codons
clist
