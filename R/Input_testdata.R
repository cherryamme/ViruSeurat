#' It will generat a Test Data
#'
#' @export
#'
#'
Input_Testdata <- function() {
  Testdata <- paste0(system.file(package = "ViruSeurat"), "/testdata")
  assign(Testmeta,read.csv("sample_info.tsv", sep = "\t", header = F, col.names = c("RAW_data", "RAW_meta"), stringsAsFactors = F),envir = .GlobalEnv)
  assign(TestSclist,Input_Sclist(meta = Testmeta, inputfiledir = Testdata),envir = .GlobalEnv)
  assign(TestVclist,Input_Vclist(inputfiledir = Testdata),envir = .GlobalEnv)
}
