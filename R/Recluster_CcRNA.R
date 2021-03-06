#' Title ReCalculate Cluster information
#'
#' @param x ccRNA
#' @param res new resolution
#'
#' @return new Recluster ccRNA
#' @export
#' @importFrom ggraph guide_edge_colourbar
#' @importFrom stringr str_detect
#' @importFrom magrittr %>%
#' @import clustree
#'
Recluster_CcRNA <- function(x = ccRNA, res = c(seq(.1, 1.6, .2))) {
  if (missing(x)) {
    message("Missing x, Use default scRNA")
    if (!exists("scRNA")) {
      stop("didn't find scRNA, Please check scRNA is exists or not, scRNA is generated by Normalize_IntegrateSclist, it is a legal integrated Seurat object. ")
    } else {
      checkmate::assert_class(scRNA, "Seurat")
    }
  } else {
    checkmate::assert_class(scRNA, "Seurat")
  }
  if (missing(res)) {
    message("Missing res, Use defalut res=c(0.1,0.3,0.5,0.7,0.9,1.1,1.3,1.5)")
  } else {
    checkmate::assertDouble(res)
  }

  message("\n \n \nClean previous find_Cluster record...")
  idn <- names(x@meta.data) %>% str_detect("integrated_snn_res.")
  x@meta.data[idn] <- NULL
  message("\nClean successfully !!  Use new res to calculate Cluster")
  x <- FindClusters(
    object = x,
    resolution = res
  )
  message("\n \n \nUse clustree Visualize every node in reduction, Result save as Clustree_res.jpg in fig/")
  Check_figdir()
  p1 <- clustree(x@meta.data, prefix = "integrated_snn_res.")
  ggsave("fig/Clustree_Reres.jpg", p1, width = 18, height = 8, limitsize = F)

  message("\n \nSave ccRNA-ReCluster object as ccRNA_Cluster.Rdata (Rdata)")
  save(x, file = "ccRNA_ReCluster.Rdata")
  message("\n \n \n       --------------ReCluster ccRNA Successfully------------")
  return(x)
}
