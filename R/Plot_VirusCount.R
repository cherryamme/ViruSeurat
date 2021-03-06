#' Title Read number distribution histogram for drawing the imported viral data
#'
#' @param x vclist you want to count from
#' @param meta metadata
#'
#' @description defaultly,it will automatic detect Data (vclist/metadata) in global env.
#' This function is used to plot virusCount in your dataset, you can visualize your VirusCounts distribution .
#' p1<-Plot_VirusCount()  # it will automatic detect Data (vclist/metadata) in global env.
#' Plot_VirusCount(x=vclist,meta=metadata) # set your own parameter .
#' @return wrap-plot-list
#' @export
#' @importFrom utils read.csv setTxtProgressBar str txtProgressBar
#'
Plot_VirusCount <- function(x = vclist, meta = metadata) {
  if (missing(x)) {
    message("\nMissing X , Use default vclist if there is")
    if (!exists("vclist")) {
      stop(
        "No VCLIST is found, retrieve whether Vclist exists.VCLIST should be generated by Input_vclist, a compliant list object with Virus feature."
      )
    } else {
      message("\nFind vclist, check if it is valid ....")
      if (class(vclist) == "vclist") {
        message("\nThis vclist is generated by Input_Vclist, it is legal and valid.")
      } else {
        stop("This vclist is not generated by Input_Vclist, Please check the input X (vclist)")
      }
    }
  }
  if (missing(meta)) {
    message("\nMissing metadata, Detecting metadata in the global environment")
    meta <- Check_meta()
  }


  message("\n \nCheck if Raw_data in VCLIST and META correspond...")
  if (!all(is.element(names(x), meta$RAW_data))) {
    stop("The elements in the vclist are not corresponding to the Raw_Data in Meta, please check the input")
  }
  cat("\nPass Check!, start the figure")

  p <- list()
  idm <- match(names(x), meta$RAW_data)
  pd <- txtProgressBar(style = 3)
  for (i in 1:length(x)) {
    cat("\n Starting process No.", i, " .....", i, "/", length(x), "\n")
    setTxtProgressBar(pd, i / length(x))
    p[i] <-
      DataExplorer::plot_histogram(x[[i]], title = paste0(names(x)[i], "__", meta$RAW_meta[i]))
  }
  close(pd)
  q <- patchwork::wrap_plots(p)
  Check_figdir()
  message("\nSave the generated map as Viruscounts_histogram.jpg, saved under the fig/ folder.")
  ggplot2::ggsave(
    "fig/VirusCounts_histogram.jpg",
    q,
    width = 18,
    height = 12,
    limitsize = F
  )
  message("\n \nComplete All !!")
  return(q)
}
