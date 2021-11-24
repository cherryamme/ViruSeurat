#' Title Normalize and Integrate Svclist, generate a scRNA : an integrated Seurat object
#'
#' @param x svclist which has virus counts, it is a sclist object
#' @param features use to find High variable genes
#' @param fmethod method use to calculate VariableFeatures
#'
#' @return scRNA
#' @export
#'
Normalize_IntegrateSclist<-function(x=svclist,fmethod="vst",features = 2000){
  # check param x
  if (missing(x)) {
    message('Missing x, Use default svclist')
    if (!exists("svclist") || !inherits(svclist, "sclist")){
      stop("Not find a svclist, Please check if svclist is exists. svclist is generat by viralinSeurat, it is a Seurat-list object and has Virus counts feature")
    }
  }
  checkmate::assertCount(features)
  cat("Starting Normalize, High variable gene discovery method is ",fmethod,", gene features = ",features,"\n")
  scNlist<-list()
  pd<-txtProgressBar(style = 3)
  for (i in 1:length(x)) {
    cat("\n    Process No.",i," element.........",i,"/",length(x),"\n")
    setTxtProgressBar(pd,i/length(x))
    cat("\n \n")
    scNlist[[i]] <- NormalizeData(x[[i]])
    scNlist[[i]] <- FindVariableFeatures(x[[i]], selection.method = fmethod,nfeatures = features)
    message("\n\n ~~~~ No.",i," element is completed ~~~~~\n \n")
  }
  close(pd)
  nfeatures <- SelectIntegrationFeatures(object.list = scNlist)
  message("Select Integration Features, Use feature:  ",nfeatures,"....\n")
  scRNA.anchors <- FindIntegrationAnchors(object.list = scNlist,anchor.features = nfeatures)
  message("Find anchor Complete, keep the anchor feature as scRNA.anchors.Rdata (Rdata)")
  save(scRNA.anchors,file="scRNA.anchors.Rdata")
  message("Use anchors to integrate Normalized scList ....\n")
  scRNA <- IntegrateData(anchorset = scRNA.anchors)
  message("\n    Integrate Data successfully, save integrate data as scRNA_intergra.Rdata (Rdata)  \n ")
  save(scRNA,file="scRNA_intergra.Rdata")
  return(scRNA)
}
