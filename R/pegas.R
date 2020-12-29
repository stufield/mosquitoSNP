
#' Borrowed functions from the pegas package
#' But heavily modified/simplified


#' @importFrom tibble as_tibble
#' @noRd
as.loci2 <- function(x, allele.sep = "/|", col.pop = NULL, col.loci = NULL) {

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
    x$population %<>% factor()
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
  x <- checkOrderAlleles(x)
  x <- tibble::as_tibble(x)
  structure(x, class = c("loci2", class(x)), locicol = col.loci)
}

#' @noRd
checkOrderAlleles <- function (x) {

  # internal closure
  reorderAlleles <- function(x) {
    for ( i in seq_along(x) ) {
      y <- x[i]
      if (!length(grep("/", y)))
        next
      y <- unlist(strsplit(y, "/"))
      y <- paste(.sortAlleles(y), collapse = "/")
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
#' @importFrom purrr map set_names
#' @export
summary.loci2 <- function(object, ...) {
  names(object)[attr(object, "locicol")] %>%
    purrr::set_names() %>%
    purrr::map(~ {
      .vec <- as.character(object[[.x]])
      list(genotype = c(table(.vec)),
           allele   = strsplit(.vec, "[/|]") %>% unlist() %>% table() %>% c()
      )
    }) %>%
    structure(class = c("summaryLoci", "list"))
}

#' @noRd
#' @importFrom usethis ui_value ui_todo ui_line
#' @importFrom purrr walk
#' @export
print.summaryLoci <- function (x, ...) {
  purrr::walk(names(x), ~ {
    usethis::ui_todo("Locus {ui_value(.x)}:")
    usethis::ui_line("-- Genotype frequencies:")
    print(x[[.x]]$genotype)
    usethis::ui_line("-- Allele frequencies:")
    print(x[[.x]]$allele)
    cat("\n")
  })
  invisible(x)
}
