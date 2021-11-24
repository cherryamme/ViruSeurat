#' Title Used to check the proportion of viruses in SCRNA data
#'
#' @param vct vclist input
#' @param sct sclist input
#'
#' @description {in: vclist} {out: print(table)}
#'
#' @author Jesse
#' @export

CompareVirusCB <- function(vct = vclist, sct = sclist) {
  UseMethod("CompareVirusCB")
}


#' Title single vclist compare
#'
#' @param vct single vclist input
#' @param sct sclist input
#'
#' @return table
#' @export
#'
CompareVirusCB.vcsingle <- function(vct, sct = sclist) {
  if (missing(vct)) {
    stop("Please select the read data set for the virus, otherwise automatically exits")
  } else{
    message("Input dataset is  ", substitute(vct), "\n")
  }
  if (missing(sct)) {
    message("Use default sclist")
    if (!exists("sclist")) {
      stop("Check the input SCLIST, check the transposon SCLIST data")
    }
  }
  if (!inherits(sct, "sclist")) {
    stop("SCLIST is not SCLIST format, you need to import by Input2Seurat, please check the transposon SCLIST data")
  }
  CompareV <- sclistbase(sct)
  CompareV(vct)
}




#' Title vclist compare
#'
#' @param vct vclist input
#' @param sct sclist input
#'
#' @return table
#' @export
#'
#'
CompareVirusCB.vclist <- function(vct, sct = sclist) {
  if (missing(vct)) {
    stop("Please select the read data set for the virus, otherwise automatically exits")
  } else{
    message("Input dataset is ", substitute(vct), "\n")
  }
  if (missing(sct)) {
    message("Use default sclist")
    if (!exists("sclist")) {
      stop("Check the input SCLIST, check the transposon SCLIST data")
    }
  }
  if (!inherits(sct, "sclist")) {
    stop("SCLIST is not SCLIST format, you need to import by Input2Seurat, please check the transposon SCLIST data")
  }
  CompareV <- sclistbase(sct)
  for (i in seq_along(vct)) {
    cat("\n \n \nNo ", i, "  element----  ", names(vct)[[i]], " Compare:\n \n")
    CompareV(vct[[i]])
  }
}


sclistbase <- function(y) {
  function(x) {
    for (i in 1:length(y)) {
      cat("No", i)
      cat("Sample-name:", as.character(unique(y[[1]]@active.ident)))
      cat(" ", as.character(unique(y[[i]]$orig.ident)))
      cat("    All-cell-num:")
      cat(length(colnames(y[[i]])))
      cat("    Has-virus-cell:")
      cat(length(names(x)))
      cat("    Compared-cell:")
      cat(length(intersect(names(x), colnames(y[[i]]))), "\n")
    }
  }
}
