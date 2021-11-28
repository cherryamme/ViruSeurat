#' It will generat a Test Data
#'
#' @param pos enviroment position
#'
#' @export
#'
#'
#'
Input_Testdata <- function(pos = 1) {
  Testdata <- paste0(system.file(package = "ViruSeurat"), "/testdata/")
  assign("Testmeta", read.csv(paste0(Testdata, "sample_info.tsv"), sep = "\t", header = F, col.names = c("RAW_data", "RAW_meta"), stringsAsFactors = F), envir = as.environment(pos))
  assign("TestSclist", Input_Sclist(meta = Testmeta, inputfiledir = Testdata), envir = as.environment(pos))
  assign("TestVclist", Input_Vclist(inputfiledir = Testdata), envir = as.environment(pos))
}