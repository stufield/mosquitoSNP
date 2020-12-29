
#' Borrowed functions from the pegas package

#' @noRd
convert2loci <- function(x, allele.sep = "/|", col.pop = NULL,
                         col.loci = NULL) {

  if ( is.null(col.pop) ) {
    ipop <- which(tolower(names(x)) == "population")
    if ( length(ipop) ) {
      col.pop <- ipop
    }
  }
  if ( is.character(col.pop) ) {
    col.pop <- which(names(x) == col.pop)
  }
  if ( is.numeric(col.pop) ) {
    names(x)[col.pop] <- "population"
    x[, col.pop] <- factor(x[, col.pop])
  }
  if ( is.null(col.loci) ) {
    col.loci <- 1:ncol(x)
    if ( is.numeric(col.pop) ) {
      col.loci <- col.loci[-col.pop]
    }
  }
  if ( is.character(col.loci) ) {
    col.loci <- match(col.loci, names(x))
  }
  if ( allele.sep != "/|" ) {
    if ( allele.sep == "" ) {
      stop("alleles within a genotype must be separated")
    }
    for ( i in col.loci ) {
      levels(x[, i]) <- gsub(allele.sep, "/", levels(x[, i]))
    }
  }
  class(x) <- c("loci", "data.frame")
  attributes(x)$locicol <- col.loci
  checkOrderAlleles(x)
}

#' @noRd
checkOrderAlleles <- function (x) {

  reorderAlleles <- function(x) {
    for ( i in seq_along(x) ) {
      y <- x[i]
      if (!length(grep("/", y))) 
        next
      y <- unlist(strsplit(y, "/"))
      y <- paste(.sortAlleles(y), collapse="/")
      x[i] <- y
    }
    return(x)
  }

  for ( k in attributes(x)$locicol ) {
    y <- x[, k]
    if ( is.numeric(y) ) {
      x[, k] <- factor(y)
      next
    }
    lv <- levels(y)
    if ( !length(grep("/", lv)) ) {
      next
    }
    a <- reorderAlleles(lv)
    if ( !identical(a, lv) ) {
      levels(x[, k]) <- a
    }
  }
  return(x)
}

#' @noRd
.sortAlleles <- function (x) {
  locale <- Sys.getlocale("LC_COLLATE")
  if ( !identical(locale, "C") ) {
    Sys.setlocale("LC_COLLATE", "C")
    on.exit(Sys.setlocale("LC_COLLATE", locale))
  }
  sort(x)
}

#' @noRd
summaryLoci <- function(object) {
  L <- attributes(object)$locicol
  ans <- vector("list", length(L))
  names(ans) <- names(object[L])
  ii <- 1L
  for ( i in L ) {
    geno <- levels(object[, i])
    alle <- strsplit(geno, "[/|]") 
    unialle <- sort(unique(unlist(alle))) 
    l <- tabulate(object[, i], length(geno))
    names(l) <- geno
    tab <- matrix(0, length(unialle), length(geno),
                  dimnames=list(unialle, geno))
    for ( j in seq_along(alle) ) {
      for ( k in alle[[j]] ) {
        tab[k, j] <- tab[k, j] + 1
      }
    }
    ans[[ii]] <- list(genotype=l, allele=drop(tab %*% l))
    ii <- ii + 1L
  }
  class(ans) <- c("summaryLoci", class(ans))
  return(ans)
}

#' @noRd
print.summaryLoci <- function (x, ...) {
  nms <- names(x)
  for (i in 1:length(x)) {
    cat("Locus", nms[i], ":\n")
    cat("-- Genotype frequencies:\n")
    print(x[[i]][[1]])
    cat("-- Allele frequencies:\n")
    print(x[[i]][[2]])
    cat("\n")
  }
  invisible(x)
}

