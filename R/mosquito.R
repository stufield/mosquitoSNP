#' Translate Codons into Amino Acids
#'
#' A function to calculate possible Amino Acid codes
#' from a codon containing  ambiguity codes.
#'
#' Translating a diploid DNA sequence into an AA sequence.
#'
#' @param x A character string containing ambiguity codes for a codon.
#' @return A list containing:
#' \item{Codons}{A list of possible codons resulting from the ambiguous codon}
#' \item{AAs}{A character string of the possible amino acids that could result,
#'   written as single letter AA codes}
#' @author Stu Field and Bill Black IV
#' @seealso [DNAtranslate()], [ambiguityCodons()], [codons()]
#' @references Bill Black IV
#' @examples
#' AAtranslate("WYR")
#' @export AAtranslate
AAtranslate <- function(x) {

  # Separate variable codon into 3 chunks
  cut_codon <- substring(x, 1:3, 1:3)
  L1 <- ifelse(any(cut_codon[1] == c("A", "G", "C", "T")), 1, 2) # A,G,T,C = 1
  L2 <- ifelse(any(cut_codon[2] == c("A", "G", "C", "T")), 1, 2) # R,S,M,Y,W,K = 2
  L3 <- ifelse(any(cut_codon[3] == c("A", "G", "C", "T")), 1, 2)
  poss_codons <- rep(NA, L1 * L2 * L3)   # Setup storage variable
  count <- 0

  # Calculate all possible codons
  for (i in 1:L1) {
    for (j in 1:L2) {
      for (k in 1:L3) {
        count <- count + 1L
        poss_codons[count] <- paste0(
          mosquitoSNP::ambiguityCodons[cut_codon[1L]][i, 1L],
          mosquitoSNP::ambiguityCodons[cut_codon[2L]][j, 1L],
          mosquitoSNP::ambiguityCodons[cut_codon[3L]][k, 1L])
      }
    }
  }

  # Translate all possible codons
  L4  <- length(poss_codons)
  AAs <- character(L4)
  nr  <- nrow(mosquitoSNP::codons)
  for ( i in 1:L4 ) {
    AAs[i] <- colnames(
      mosquitoSNP::codons[ceiling(which(mosquitoSNP::codons == poss_codons[i]) / nr)]
      )
  }

  # Condense duplicate AAs and create single entry output
  list(Codons = poss_codons,
       AAs    = paste(names(table(AAs)), collapse = ""))
}


#' DNA Translate
#'
#' Translate a sequence of DNA nucleotides to AA sequence.
#' Ambiguity codes are allowed.
#'
#' @param x A vector of DNA sequence for translation into AA sequence.
#' @param csv Logical. Create a csv output file?
#' @param file Character. Name of the csv output file if desired
#' @return A list containing:
#' \item{DNAsequence}{Original DNA sequence}
#' \item{AAsequence}{Translated AA sequence based on DNA sequence}
#' \item{TranslationTable}{Table combining in nice format}
#' @author Stu Field, Martin Donnelly, Bill Black
#' @seealso \code{\link{AAtranslate}}, \code{\link{searchReplace}}
#' @examples
#' chunk1 <- c("YYYAAGATCGTCGGTGGCGATGAGGCCGAAGCGCACGAATTTCCCTACCAAATCTCG")
#' chunk2 <- c("CTGCAGTGGAACTTCAACGATGGACAAACGGAGACCATGCACTTCTGYGGAGCTTCGG")
#' chunk3 <- c("TGTTGAACGAAAACTTYGTCCTGACGGCTGCTCACTGCAAGACCGCATACTCCAATA")
#' chunk4 <- c("CCGGGTWCATCGAAGTGGTTGCCGCTGAACATGATGTGGCYGTTGCGGAAGGATCCGA")
#' chunk5 <- c("ACAGCGTCGYTTGGTTGCGGAGTTCATCGTCCACGAGGACTATCAAGGRGGAGTCAGT")
#' chunk6 <- c("CCCGATGAGATTGCCGTCA")
#' DNAseq <- paste0(chunk1, chunk2, chunk3, chunk4, chunk5, chunk6)
#' DNAseq
#'
#' AAtranslate("YRY")
#' DNAtranslate(DNAseq, csv = FALSE)
#' @importFrom utils write.csv
#' @export DNAtranslate
DNAtranslate <- function(x, csv = FALSE, file) {

  ntds      <- nchar(x)
  codon.vec <- substring(x, seq(1, ntds, 3), seq(3, ntds, 3))
  L         <- floor(ntds / 3)
  AASeq     <- character(L)
  for ( i in 1:L ) {
    AASeq[i] <- AAtranslate(codon.vec[i])$AAs
  }

  tabl <- as.data.frame(rbind(AASeq, codon.vec[1:L]))
  names(tabl) <- 1:L
  rownames(tabl) <- c("AAsequence", "DNAsequence")

  if ( csv ) {
    write.csv(tabl, file = sprintf("%s_%s.csv", file, Sys.Date()))
  }

  list(DNAsequence      = x,
       AAsequence       = paste(AASeq, collapse = "-"),
       TranslationTable = tabl)
}


#' Are All Entries Identical?
#'
#' Determine if all entries of a vector are identical. Can work on a matrix
#' also because matrix could be converted to a vector.
#'
#' @param x A numeric or character vector or a matrix that could be converted
#' to a vector
#' @return A logical value indicating if all values in the vector are identical
#' @author Stu Field
#' @seealso \code{\link{identical}}
#' @examples
#' allSame(1:15)
#' allSame(rep("DOG", 10))
#' allSame(diag(5))
#' @export allSame
allSame <- function(x) {
  all(isTRUE(x == x[1L]))
}


#' Array Data to Data Frame
#'
#' Reorganize array data into a data frame. Uses \code{\link{cbind}} to combine
#' arrays horizontally with Probe ID running down column 1. Each array
#' red-signal green-signal weighting for probe running horizontally for each of
#' the arrays
#'
#' Bind all the arrays together horizontally i.e. probe id running down
#' column 1 and each array red signal-green signal-weighting for probe
#' running horizontally for each of the arrays
#'
#' Takes as input x, the RG.corrected limma class object
#' Meta = 1) MetaRow, 2) MetaCol as included in original R-code
#'
#' @param x The "RGList" class object that is the result of
#' a call to [limma::RG.MA()].
#' @param names The column names of `x` that you wish to
#' include in the order you want them.
#' @param Meta The row (1) and col (2) of MetaRow and MetaCol
#' as seen in original R-code.
#' Not always used so default is set to c(1, 1).
#' @param probechar The cut off for the length of the probe name (default = 13)
#' @return A data frame containing all the desired array data.
#' @author Stu Field
#' @seealso [limma()]
#' @examples
#' RG.corrected$targets
#' colNames <- c("ProbeName", "GeneName", "SystematicName",
#'               "Description", "Status", "Sequence", "Row", "Col")
#' Tororo.raw <- array2df(RG.corrected, names = colNames)
#' head(Tororo.raw, 25)
#' @export array2df
array2df <- function(x, names, Meta = c(1, 1), probechar = 13) {

  arrays <- dim(x)[2L]
  for ( i in 1:arrays ) {
    if ( i == 1 ) {
      X <- cbind(x[, 1]$genes[names], MetaRow = Meta[1L], MetaCol = Meta[2L])
      L <- ncol(X)
    }
    X <- cbind(X, x[, i]$R, x[, i]$G, 1 - x[, i]$weights)
  }
  X <- as.data.frame(X)
  n <- ncol(X)
  append <- c("R", "G", "W")    # append colnames with R, G, or W
  for ( j in 1:length(append) ) {
    names(X)[seq(L+j, n, 3)] <- paste(names(X)[seq(L+j, n, 3)], append[j], sep = ".")
  }
  cat("* Combining", i, "arrays\n")
  X$ProbeName <- substr(X$ProbeName, 1, probechar) # cut off ProbeName to max 13 char (default)
  return(X)
}


#' Bootstrap Phenotypes
#'
#' Create bootstrap population of phenotypes
#'
#' @param x The original vector of population data
#' @param loci The loci present in the population
#' @param pops A factor variable of the populations
#' @param phenotypes Phenotypes of the population. Usually R/S
#' @param rseed Optional for reproducibility
#' @return A bootstrap population based on the original population data
#' @author Stu Field, Martin Donnelly
#' @examples
#' # add example
#' \dontrun{}
#' @export boot.pop
boot.pop <- function(x, loci, pops, phenotypes, rseed = sample(1000, 1)) {

  set.seed(rseed)
  I <- length(loci)         # column with locus as header
  J <- length(pops)         # no. of populations
  K <- length(phenotypes)   # no. of phenotypes
  # Resample with replacement
  for ( i in 1:I ) {
    for ( j in 1:J ) {
      for ( k in 1:K ) {
        boot.sample <- x[ x$Population == pops[j] & x$Phenotype == phenotypes[k], loci[i] ]
        x[x$Population == pops[j] & x$Phenotype == phenotypes[k], loci[i]] <- sample(boot.sample, replace = TRUE)
       }
    }
  }
  return(x)
}



#' Find and Calculate SNPs
#'
#' Determine the relative positions of SNPs for variable SNP data.
#' Individuals are rows and nucleotide positions are cols.
#' Must be a data frame or a matrix. See Examples.
#'
#' From Diploid genetic sequences.
#'
#' @param x A matrix of raw SNP data, individuals as rows & loci as cols.
#' @return A list containing:
#' \item{Data }{The raw data reorganized as a matrix. Individuals as rows.
#' Loci as cols}
#' \item{Individuals }{Number of individuals in the data set. Equal to
#' nrows}
#' \item{Length}{Length of the SNP sequence. Equal to ncol}
#' \item{nSNP}{Number of SNPs in the data set}
#' \item{SNPs}{Position of the SNP loci, identified by the col number}
#' @author Stu Field, Bill Black
#' @seealso [allSame()]
#' @keywords SNP
#' @examples
#' Ago2
#' findSNPs(Ago2)
#' @export findSNPs
findSNPs <- function(x) {

  ntd <- nchar(x)
  X <- matrix(NA, nrow = nrow(x), ncol = ntd,
              dimnames = list(c(rownames(x)), c(1:ntd)))

  # separate characters into vector entries
  for ( i in 1:nrow(x) ) {
    X[i,] <- substring(x[i,], 1:ntd, 1:ntd)
  }

  # Remove dashed individuals
  X <- X[-which(X[, 1L] == "-"), ]

  # which sites are variable
  SNPs <- unname(which(!apply(X, 2, allSame)))

  list(Data = X,
       Individuals = nrow(X),
       Length = ntd,
       nSNP = length(SNPs),
       SNPs = SNPs)
}



#' Convet Microsatellite Data Format
#'
#' Convert adjacent microsatellite data from separate columns for each of the
#' diploid loci to a "/" separating the alleles (i.e. allele1/allele2).
#'
#' @param rawdata Raw microsattellite data containing separate columns for each
#' of the diploid loci. Converts to a "/" separating the diploid alleles.
#' @param col.loci Defining the columns which contain microsat data.
#' @param locus.names Names of the microsat loci
#' @return Matrix of microsatellite data with "/" as the allele separator,
#' rather than separate cols for alleles.
#' @author Stu Field
#' @seealso [HarrCraigData]
#' @references "inst/HarrCraig.csv"
#' @examples
#' head(HarrCraigData)
#' Locusnames <- c("L3RIH59", "L3RIND30", "L3RIR01", "L3RIR02",
#'                 "L3RIR03", "L3RIR04", "L3RIR05", "L3RIR07",
#'                 "L3RIR08", "L3RIR09", "L3RIR10", "L3RIR12")
#'
#' # Reformat
#' FormattedData <- formatMicrosat(HarrCraigData, col.loci = 4:27,
#'                                 locus.names = Locusnames)
#' head(FormattedData)
#'
#' @importFrom utils tail
#' @export formatMicrosat
formatMicrosat <- function(rawdata, col.loci, locus.names) {

  col.1    <- col.loci[1L]              # first locus col
  col.x    <- tail(col.loci, 1)        # final locus col
  Y        <- rawdata[, 1:(col.1 - 1)]
  locibind <- seq(col.1, col.x - 1, 2)   # loci to bind

  if ( length(locibind) %% 2 != 0 ) {
    stop("Odd number of Loci to bind")
  }

  # paste & bind #
  for ( i in locibind ) {
    if ( i == locibind[1L] ) {
      X <- paste(rawdata[, i], rawdata[, i + 1], sep = "/")
    } else {
      X <- cbind(X, paste(rawdata[, i], rawdata[, i + 1], sep = "/"))
    }
  }
  colnames(X) <- locus.names
  cbind(Y, X)
}



#' Heterozygosity Bootstrapping
#'
#' Conduct bootstrap analysis on heterozygosity microsatellite data for
#' resistance/susceptible phenotypes. Originally written for the analysis of
#' Harr Craig data set.
#'
#' Bootstrapping routine on original MicroData: HarrCraigData
#'
#' @param x Microsatellite data object (usually data frame). See example.
#' @param nboot Number of bootstrap replicates.
#' @param LociCols Which cols are microsatellite loci.
#' @return A list containing:
#' \item{estimate }{The empirical estimate of heterozygosity by population
#' and by phenotype}
#' \item{upperCI9 }{upper confidence limit of the heterozygosity estimate}
#' \item{lowerCI95 }{lower confidence limit of the heterozygosity estimate}
#' @note Bootstrapping routine on original MicroData
#' @author Stu Field, Martin Donnelly
#' @seealso \code{\link{HeFUN}}, \code{\link{boot.pop}}
#' @references Harr Craig microsatellite data
#' @examples
#' Locusnames <- c("L3RIH59", "L3RIND30", "L3RIR01", "L3RIR02",
#'                 "L3RIR03", "L3RIR04", "L3RIR05", "L3RIR07",
#'                 "L3RIR08", "L3RIR09", "L3RIR10", "L3RIR12")
#'
#' MicroData <- formatMicrosat(HarrCraigData, col.loci = 4:27, locus.names = Locusnames)
#'
#' # exclude Alajo and Malawi because:
#' # no Alajo|S and no Malawi|R
#' # cannot estimate heterozygosity with 1 phenotype
#' MicroData <- subset(MicroData, Population != "Alajo" & Population != "Malawi")
#' head(MicroData)
#'
#' HeBootstrapFUN(MicroData, nboot = 100, LociCols = 4:15)
#'
#' @importFrom stats quantile
#' @export HeBootstrapFUN
HeBootstrapFUN <- function(x, nboot, LociCols) {  # x = Microsat object
   # general housekeeping #
  Popns     <- levels(factor(x$Population))  # used in boot.pop() below
  types     <- levels(factor(x$Phenotype))   # used in boot.pop()
  lociNames <- colnames(x[LociCols])         # used in boot.pop()
  nPops     <- length(Popns)
  ntypes    <- length(types)
  L         <- length(lociNames)

  # Storage Array for Bootstrap Estimates of He ###
  HeArray <- array(NA, dim = c(L, nPops * ntypes, nboot))
  #print(dim(HeArray))

  # convert to loci; necessary for summary method nested within by()
  y     <- convert2loci(x, col.pop = 2, col.loci = LociCols)
  Sy    <- by(y, INDICES = list(y$population, y$Phenotype), FUN = summaryLoci)
  EstHe <- sapply(Sy, function(i) sapply(i, function(j) HeFUN(j$allele)))
  colnames(EstHe) <- paste(rep(Popns, ntypes), rep(types, each = nPops), sep = "|")
  EstHe <- EstHe[, order(colnames(EstHe)) ]  # reorder columns alpha by types

  # Calculate He for each bootstrap popn
  for ( n in 1:nboot ) {
    if ( n%%10 == 0 ) {
      cat("Bootstrap ...", n, "\n")
    }

    # Create bootstrapped dataset with replacement
    # and convert to "loci" object
    Boot.x <- boot.pop(x, loci = lociNames, pops = Popns, phenotypes = types, rseed = n)
    Boot.y <- convert2loci(Boot.x, col.pop = 2, col.loci = LociCols)

    # Summarize MicroData by Popn, Phenotype, & locus
    # This is a list of each combination or Popn & Phenotype
    # with each element of the list containing a list of each Locus
    S  <- by(Boot.y, INDICES = list(Boot.y$population, Boot.y$Phenotype), FUN = summaryLoci)
    He <- sapply(S, function(i) sapply(i, function(j) HeFUN(j$allele)))
    colnames(He) <- paste(rep(Popns, ntypes), rep(types, each = nPops), sep = "|")
    He <- He[, order(colnames(He)) ]        # reorder columns alpha by types
    HeArray[,, n]     <- He                 # Store He for boot popn
    dimnames(HeArray) <- dimnames(EstHe)    # Conserve col/row names among EstHe and HeArray for indexing
  }

  # Calculate CI95s
  lower <- apply(HeArray, 1:2, quantile, probs = 0.025) # c(1, 2) = by row & by col
  upper <- apply(HeArray, 1:2, quantile, probs = 0.975)

  list(estimate  = EstHe,
       upperCI95 = upper,
       lowerCI95 = lower)
}



#' Calculate Heterozygosity
#'
#' Calculate a heterozygosity estimate given a set of microsatellite data.
#'
#' This function is based on the Heterozygosity estimate function within the
#' pegas package by Emmanuel Paradis. However Martin believed there was an
#' error in the variance calculation so it was modified.
#'
#' @param x Microsatellite data.
#' @param var Should variance also be computed according to: 2 * (2 * (n - 2) *
#' (sp3 - sp2^2) + sp2 - sp2^2)/(n * (n - 1))
#' @return A list containing:
#' \item{H }{Point estimate of population heterozygosity}
#' \item{var.H }{Variance associated with the point estimate}
#' \item{n }{Number of samples used in the heterozygosity point estimate}
#' @author Stu Field, Martin Donnelly
#' @seealso [HeBootstrapFUN()]
#' @references Based on the Heterozygosity function in the pegas package by
#' Emmanuel Paradis [pegas::pegas()].
#' @examples
#' Locusnames <- c("L3RIH59", "L3RIND30", "L3RIR01", "L3RIR02",
#'                 "L3RIR03", "L3RIR04", "L3RIR05", "L3RIR07",
#'                 "L3RIR08", "L3RIR09", "L3RIR10", "L3RIR12")
#'
#' MicroData <- formatMicrosat(HarrCraigData, col.loci = 4:27, locus.names = Locusnames)
#'
#' apply(MicroData[, 4:ncol(MicroData)], 2, HeFUN)
#' apply(MicroData[, 4:ncol(MicroData)], 2, HeFUN, var = TRUE)
#' @export HeFUN
HeFUN <- function(x, var = FALSE) {

  if ( !is.factor(x) ) {
    if ( is.numeric(x) ) {
      x <- x[ names(x)!="NA" ] # remove NAs
      n <- sum(x)
      k <- length(x)
      freq <- x / n
    } else {
      x <- factor(x)
    }
   }

  if ( is.factor(x) ) {
    n    <- length(x)
    k    <- nlevels(x)
    freq <- table(x) / n
  }
  sp2  <- sum(freq^2)
  H    <- n * (1 - sp2) / (n - 1)

  if ( var ) {
    sp3   <- sum(freq^3)
    var.H <- 2 * (2 * (n - 2) * (sp3 - sp2^2) + sp2 - sp2^2) / (n * (n - 1))
    return(c(H, var.H, n))
  } else {
    return(H)
  }
}



#' Lineage Map & EHH Analysis
#'
#' Creates a lineage map of a given data set as part of the Extended Haplotype
#' Heterozygosity (EHH) analysis.
#'
#' @param data Raw data of nucleotide changes.
#' @param core Set the nucleotide position of the core.
#' @param csv Write a csv output file of the resulting lineage map?
#' @param file Character. If `csv = TRUE`, the name of the output file.
#' Otherwise ignored.
#' @return A list containing:
#' \item{SNP.Data}{Original SNP data.}
#' \item{Haplotypes}{The number of haplotypes lineages left & right of the
#' core}
#' \item{Lineages}{Table showing the number of haplotypes in each
#' lineage}
#' \item{Haplo.Lineages}{The actual lineage map; ordered from top
#' to bottom}
#' @author Stu Field, Martin Donnelly
#' @seealso [allSame()], [SNPboot()]
#' @examples
#' names <- paste(GhanaKDR[, 1], GhanaKDR[, 2], sep = "")
#' new.data <- GhanaKDR[-c(1, 2)]
#' colnames(new.data) <- 1:ncol(new.data); rownames(new.data) <- names
#' new.data <- new.data[-59,]    # remove resistant individual
#' head(new.data)
#'
#' # replace numerics with nucleotides
#' new.data <- searchReplace(new.data, s = 1:4 , r = c("A","C","G","T"))
#'
#' # the call
#' lineageMap(new.data, core = 25)
#'
#' \dontrun{
#' lineageMap(new.data, core = 25, csv =TRUE, file = "lineageOutput.csv")
#' }
#' @export lineageMap
lineageMap <- function(data, core, csv = FALSE, file) {

  sites <- ncol(data)
  n     <- nrow(data)
  SNPnames <- paste0("SNP", seq(1-core, ncol(data) - core))
  SNPnames[core] <- "Core"
  Map   <- matrix(0, n, sites, dimnames = list(rownames(data), SNPnames))
  Map[, core] <- rep(1, n)

  # Centromeric # Left <- ----
  for ( i in (core-1):1 ) {
    if ( allSame(data[, i]) ) {
      Map[, i] <- Map[, i + 1]
      next
    } else {
      lineages  <- names(sort(table(Map[, i + 1]), decreasing = TRUE))
      nlineages <- length(lineages)

      for ( k in 1:nlineages ) {
        lineage.index <- which(Map[, i + 1] == lineages[k])
        nextSNP <- data[lineage.index, i]

        if ( allSame(nextSNP) ) {
          Map[lineage.index, i] <- Map[lineage.index, i + 1]
          next
        } else {
          htypes <- names(sort(table(nextSNP), decreasing = TRUE))
          old.v  <- which(data[lineage.index, i] == htypes[1L])
          new.v  <- which(data[lineage.index, i] == htypes[2L])
          Map[lineage.index[old.v], i] <- Map[lineage.index[old.v], i + 1]
          Map[lineage.index[new.v], i] <- max(Map[, i:(i +1 )]) + 1
          if ( length(htypes) == 3 ) {
            third <- which(data[lineage.index, i] == htypes[3L])
            Map[lineage.index[third], i] <- max(Map[, i:(i + 1)]) + 1
          }
        }
      }
    }
  }

  # Telomeric # Right -> ----
  for ( i in (core+1):sites ) {

    if ( allSame(data[, i]) ) {
      Map[, i] <- Map[, i-1]
      next
    } else {
      lineages  <- names(sort(table(Map[, i - 1]), decreasing = TRUE) )
      nlineages <- length(lineages)
      for ( k in 1:nlineages ) {
        lineage.index <- which(Map[, i - 1] == lineages[k])
        nextSNP <- data[lineage.index, i]

        if ( allSame(nextSNP) ) {
          Map[lineage.index, i] <- Map[lineage.index, i-1]
          next
        } else {
          htypes <- names(sort(table(nextSNP), decreasing = TRUE) )
          old.v <- which(data[lineage.index, i] == htypes[1L])
          new.v <- which(data[lineage.index, i] == htypes[2L])
          Map[lineage.index[old.v], i] <- Map[lineage.index[old.v], i - 1]
          Map[lineage.index[new.v], i] <- max(Map[, i:(i - 1)]) + 1
          if ( length(htypes) == 3 ) {
            third <- which(data[lineage.index, i] == htypes[3L])
            Map[lineage.index[third], i] <- max(Map[, i:(i - 1)]) + 1
          }
        }
      }
    }
  }

  if ( csv ) {
    write.csv(Map[order(Map[, 1]),], file = sprintf("%s_%s.csv", file, Sys.Date()),
              row.names = TRUE)
  }

  list(SNP.Data       = data,
       Haplotypes     = c(Left = max(Map[, 1]), Right = max(Map[, sites])),
       Lineages       = list(Left = table(Map[, 1]), Right = table(Map[, sites])),
       #Haplo.Lineage = Map,
       Haplo.Lineages = Map[order(Map[, 1]), ])
}


#' Sliding SNP Plot
#'
#' General wrapper for plotting sliding window analyses on SNP
#' chromosome base pair positions. Primary utility in [slide.window.bp()].
#'
#' @param x Midpoints. Typically the output of `slide.window.bp()$midpoints`.
#' @param y Values. Typically the output of `slide.window.bp()$values`.
#' @param breaks Base pairs where the chromosome breaks occur from one
#' chromosome to the next.
#' @param ylabel Label for the y-axis.
#' @param font Font for the plot, passed to [par()].
#' @param line.cols Line colours for plotting.
#' @param line.lty Line type for plotting
#' @param save Logical.
#' @param ... Additional arguments passed to [plot()].
#' @return A step "sliding window" plot.
#' @note Length of `breaks =`, `line.cols =`, and `line.lty =`
#' *must* all be equal.
#' @author Stu Field, Martin Donnelly
#' @seealso [slide.window.bp()]
#' @examples
#' head(Ghana_MvS_SNPdata)
#'
#' BPwin    <- 7500000
#' negLogP  <- Ghana_MvS_SNPdata$negLogP
#' data.out <- slide.window.bp(negLogP, Pos = Ghana_MvS_SNPdata$ScaledPosition,
#'                             bp = BPwin, FUN = max)
#' data.out
#'
#' chr.breaks <- c(60707344, 62812924, 110916950, 162910923, 206568235)
#' chr.cols   <- c(2, 2, 3, 2, 3)
#' chr.lines  <- c(2, 2, 1, 2, 1)
#'
#' plotSlideSNP(x = data.out$midpoints, y = data.out$values, type = "s",
#'              col = "navy",
#'              lwd = 1.25,
#'              ylabel = "-logP (max)",
#'              breaks = chr.breaks,
#'              font = 4,
#'              line.cols = chr.cols,
#'              line.lty = chr.lines,
#'              save = FALSE)
#' @importFrom graphics plot abline par box axis
#' @importFrom grDevices pdf dev.off
#' @export plotSlideSNP
plotSlideSNP <- function(x, y, breaks, ylabel, font, line.cols,
                         line.lty, save = FALSE, ...) {

  if ( length(breaks) != length(line.cols) )
    stop("Line Colours Missing!")
  if ( length(breaks) != length(line.lty) )
    stop("Line Types Missing!")
  if ( save ) {
    file <- sprintf("WindowPlot_%s.pdf", Sys.Date())
    pdf(file, height = 8, width = 12, title = "SlidingSNPwindow")
    on.exit(dev.off())
  }
  plot(y ~ x, ylab = ylabel, xlab = "Chromosome Position", ...)
  box()
  par(font = font)
  abline(v = breaks, col = line.cols, lty = line.lty, lwd = 0.75)
}





#' Sliding Window SNP Analysis (bp)
#'
#' Sliding window analysis for SNPs along chromosome based on base pair
#' positions rather than SNP positions. Could be any vector of values.
#'
#' @param X Vector object containing data appropriate to be applied by
#' `FUN`.
#' @param Pos Vector of base pair (SNP) positions to evaluate.
#' @param bp Size of the window over which the `FUN` should be evaluated
#' in number of base pairs. Equivalent to `window =` in
#' [slide.window()].
#' @param FUN The summary function that evaluates the numbers in the window,
#' typically `mean` or `max`.
#' @return A list containing:
#' \item{values}{The result of the application of `FUN` in each of the windows.}
#' \item{midpoints}{The midpoint of the window, used primarily for plotting.}
#' \item{n.SNPs}{How many SNPs evaluated in each window based on constant bp
#' jumps along the chromosome.}
#' @author Stu Field, Martin Donnelly
#' @seealso [plotSlideSNP()], [slide.window()]
#' @examples
#' r.vec <- sample(1:20, 100, replace = TRUE)  # random vector
#' slide.window.bp(X = r.vec, Pos = 1:length(r.vec), bp = 5, FUN = mean)
#'
#' head(Ghana_MvS_SNPdata)
#'
#' BPwin   <- 7500000
#' negLogP <- Ghana_MvS_SNPdata$negLogP
#' plot.data <- slide.window.bp(negLogP, Pos = Ghana_MvS_SNPdata$ScaledPosition,
#'                              bp = BPwin, FUN = max)
#'
#' chr.breaks <- c(60707344, 62812924, 110916950, 162910923, 206568235)
#' chr.cols <- c(2, 2, 3, 2, 3)
#' chr.lines <- c(2, 2, 1, 2, 1)
#'
#' plotSlideSNP(x = plot.data$midpoints, y = plot.data$values,
#'              type = "s", col = "navy",
#'              lwd = 1.25, ylabel = "-logP (max)",
#'              breaks = chr.breaks, font = 4,
#'              line.cols = chr.cols,
#'              line.lty = chr.lines,
#'              save = FALSE)
#' @export slide.window.bp
slide.window.bp <- function(X, Pos, bp, FUN) {

  M <- max(Pos)
  L <- M %/% bp
  if ( bp > M ) {
    stop("Error: bp window too large!")
  }
  out1 <- numeric(L)
  out2 <- numeric(L)
  out3 <- numeric(L)
  hi <- seq(bp, M, bp)
  lo <- hi - (bp - 1)
  for ( i in 1:L ) {
    win <- which(Pos >= lo[i] & Pos <= hi[i])
    out1[i] <- FUN(X[win])             # summary function
    out2[i] <- mean(range(Pos[win]))   # calculate position midpoint
    out3[i] <- length(win)             # calc # SNPs in window
  }
  list(values    = out1,
       midpoints = out2,
       n.SNPs    = out3)
}


###############################
# Function for Calculating
# value of a sliding window
# for simple functions
################################
# The slide.window() function
# window = size of window for analysis at each step
# frame.skip = how far to jump to next window frame
# X = vector for analysis
# FUN = summary function over window values
###########################################

#' Sliding Window SNP Analysis
#'
#' Sliding window analysis for SNPs along chromosome.
#' Can be any vector of values.
#'
#' @param X Vector object containing data appropriate to be applied by
#' \code{FUN}.
#' @param window Size of the window over which the \code{FUN} should be
#' evaluated.
#' @param frame.skip How far the window should jump from one iteration to the
#' next. Default = 1.
#' @param FUN The summary function that evaluates the numbers in the window,
#' typically \code{mean} or \code{max}.
#' @return A list containing:
#' \item{values}{The result of the application of \code{FUN} in each of the windows.}
#' \item{midpoints}{The midpoint of the window, used primarily for plotting.}
#' \item{length}{How many 'windows' were evaluated during the analysis.}
#' @author Stu Field, Martin Donnelly
#' @seealso [lotSlideSNP()], [slide.window.bp()]
#' @examples
#' # Example 1
#' r.vec <- sample(1:20, 100, replace = TRUE)       # create random vector
#' slide.window(X = r.vec, window = 5, frame.skip = 2, FUN = max)
#'
#' # Example 2
#' head(Ghana_MvS_SNPdata)
#'
#' Win <- 5
#' negLogP <- Ghana_MvS_SNPdata$negLogP
#' data.out <- slide.window(negLogP, window = Win, frame.skip = 2, FUN = max)
#'
#' plot(data.out$values ~ data.out$midpoints,
#'      type = "s", col = "navy",
#'      ylab = "-logP (max)",
#'      xlab = "midpoints", lwd = 1.5)
#' legend("topleft",
#'        legend = format(paste("WindowSize =", Win)),
#'        bg = "gray75", cex = 0.75, box.lty = 0)
#' @export slide.window
slide.window <- function(X, window, frame.skip = 1, FUN) {

  L <- length(X)

  if ( window >= L ) {
    out <- FUN(X)
    cat("Warning: window >= length(X)", "\n")
    return(out)
  }

  P  <- 1:L                    # SNP positions
  hi <- seq(window, L, frame.skip)
  lo <- hi - (window - 1)
  out1 <- numeric(length(hi))
  out2 <- numeric(length(hi))

  for (i in 1:length(hi)) {
    win <- lo[i]:hi[i]
    if (win[window] > L)
      break     # break if window hangs over end
    out1[i] <- FUN(X[win])
    out2[i] <- mean(range(P[win]))   # calculate position midpoint
  }

  list(values    = out1,
       midpoints = out2,
       length    = length(out1))
}



#' SNP Bootstrap of EHH data
#'
#' Perform EHH analysis on SNP data and calculate CI95
#' via bootstrapping method.
#'
#' 1st col *must* be character class!
#'
#' @param x A data frame or matrix with haplotype lineage data, typically
#' output from [lineageMap()]. Column 1 MUST be a character vector,
#' usually with an Amino Acid string describing the locus.
#' @param nBoot Number of bootstrap samples to perform for each SNP locus.
#' @param core Set which column (SNP) should be considered the `core`.
#' Typically determined from [lineageMap()].
#' @param lo Lower confidence limit passed to [quantile()].
#' @param up Upper confidence limit passed to [quantile()].
#' @param pts Size of the points for the EHH plot.
#' @param hist Logical. Should a histogram also be plotted showing the
#' distribution of the bootstraps of EHH?
#' @param plot Logical. Should the Expected Haplotype Homozygosity be plotted?
#' @param csv Logical. Should the data frame containing EHH point estimate &
#' CI95 bootstrap errors be saved to a file?
#' @param file Character. If `csv = TRUE`, the name of the *.csv file to be
#' produced. Otherwise ignored.
#' @return A data frame (with `nrow = nSNP`), consisting of:
#' \item{Est H-H}{the EHH point estimate}
#' \item{2.5\%}{the lower CI95}
#' \item{97.5\%}{the upper CI95}
#'
#' This data frame is then used by the [gplots::plotCI()] function
#' plot the EHH plot with bootstrap CI95 errors.
#' @note Requires the \code{gplots} package.
#' @note If the bootstrap confidence interval is very short
#' relative to the point estimate, [gplots::plotCI()] will give the
#' warning: `In arrows(...): zero-length arrow is of indeterminate
#' angle and so skipped`.
#' @author Stu Field, Martin Donnelly
#' @seealso \code{\link{boot.pop}}, [HeBootstrapFUN()],
#' gplots::plotCI()], [lineageMap()], [quantile()]
#' @references "Kenyan_SNP_EHH_Data.csv"
#' @examples
#'
#' head(KenyaSNP_EHHdata)   # pkg data
#'
#' # Subset data sets by Leu/Ser
#' data.Leu <- subset(KenyaSNP_EHHdata, KenyaSNP_EHHdata == "Leu")
#' data.Ser <- subset(KenyaSNP_EHHdata, KenyaSNP_EHHdata == "Ser")
#'
#' SNPboot(data.Leu, nBoot = 50, core = 21)  # nBoot >500 recommended
#' SNPboot(data.Ser, nBoot = 25, core = 21)
#'
#' @importFrom gplots plotCI
#' @export SNPboot
SNPboot <- function(x, nBoot, core, lo = 0.025, up = 0.975, pts = 19,
                    hist = FALSE, plot = TRUE, csv = FALSE, file) {

  cols <- ncol(x)   # number of SNPs
  n    <- nrow(x)   # pop size

  # Calculate Point Estimate Haplotype Homozygosity
  EstHH <- NA
  for ( z in 2:cols ) {
    EstHH[z-1] <- ( (sum(as.vector(table(x[, z]))^2) / n) - 1) / (n - 1)
  }

  # Create Bootstrap Samples
  for ( bt in 1:nBoot ) {
    if ( bt%%10 == 0 ) {
      cat("Bootstrap ...", bt, "\n")
    }

    boot_vec <- sample(1:nrow(x), replace = TRUE)   # select bootstraps indices

    if ( hist ) {
      plot(table(boot_vec))
    }

    BootSamplePop <- x[boot_vec, ]

    # Calculate Haplotype Homozygosity estimate
    EH <- numeric(cols - 1)      # placeholder for Estimated Homozygosity

    for ( j in 2:cols ) {        # 1st col is factor (skip)
      counts <- as.vector(table(BootSamplePop[, j]))
      EH[j-1] <- ( (sum(counts^2) / n) - 1) / (n - 1)
    }

    if ( bt == 1 ) {
      BootHH <- EH
    } else {
      BootHH <- rbind(BootHH, EH)
    }
  }

  # Calculate upper & lower percentiles
  BootHH <- unname(BootHH)  # Remove names from table
  lower  <- apply(BootHH, 2, quantile, probs = lo)
  upper  <- apply(BootHH, 2, quantile, probs = up)

  ret <- data.frame(cbind(EstHH, lower, upper))
  names(ret) <- c("Est H-H","2.5%","97.5%")
  rownames(ret) <- 1:(cols - 1)

  if ( csv ) {
    write.csv(ret, file = sprintf("%_%s.csv", file, Sys.Date()), row.names = TRUE)
  }

  if ( plot ) {
    plotCI(x = 1:(cols-1), y = ret[, "Est H-H"],
           li = ret[,"2.5%"], ui = ret[, "97.5%"],
           pch = pts, gap = 0, cex = 0.85, xaxt = "n", xlab = "",
           ylab = "Extended Haplotype Homozygosity")
    box()
    axis(1, at = seq(1, cols-1, by = 2), cex.axis = 0.85,
         labels = seq(-core+1, cols-1-core, by = 2))
  }
  ret
}

