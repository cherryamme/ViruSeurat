#' Title Input transposon data as Seurat-list Format
#'
#' @param cellmin the min cell you want to keep for
#' @param fmin the min feature you want to keep for
#' @param meta metadata that needs input
#' @param inputfiledir transposon data dir
#'
#' @return Seurat-list
#' @export
#' @description defaultly,it will automatic detect Data in data/ directory.
#' The metadata type is like "CRR703  H03P" generate from sample_info.tsv. This tsv-file must be in your workdirectory.
#' cellmin is the min cell you want to keep. below this value, the observation will be dropped.
#' fmin is the min feature you want to keep .below this value, the variable will be dropped.
#' sclist<-Input_Sclist() #it will automatic detect Data in data/ directory.
#' vclist<-Input_Sclist(data=metadata,cellmin =3,fmin=20) #you can use your own metadata dataset and specific cellmin or fmin.
#' @import Seurat
#' @importFrom stringr str_sub
#' @importFrom utils read.csv setTxtProgressBar str txtProgressBar

Input_Sclist <- function(meta = metadata,
                         inputfiledir = "data/",
                         cellmin = 3,
                         fmin = 20) {
  # Judge must be param
  if (missing(meta)) {
    message(
      "Missing input meta, you need to give imported data, save it in the working directory to be named ** _ info.tsv, import to metadata, try to use the default metadata"
    )
  meta <- Check_meta()
  }
  if (all(colnames(meta) == c("RAW_data", "RAW_meta"))) {
    message(
      "Detected Meta two columns of 'Raw_Data', 'Raw_meta', continued to perform import \n\n"
    )
  } else {
    stop(
      "Detected that Meta does not contain two columns 'Raw_Data', 'Raw_meta', not a predetermined format, abort process."
    )
  }
  x <- meta[["RAW_data"]]
  meta <- meta[["RAW_meta"]]
  # Judge optinal param
  if (missing(cellmin)) {
    message("Missing cellmin, Use default min.cells = 3")
  }
  if (missing(fmin)) {
    message("Missing fmin, Use default min.features = 20")
  }
  # Input data to X
  message(
    "\n \n x has ",
    length(x),
    " elements, Will be imported sequentially, use the read.csv parameter: header = T,row.names = 1 \n"
  )
  pb <- txtProgressBar(style = 3)
  for (i in 1:length(x)) {
    cat(
      "\n Inporting No.",
      i,
      " element............",
      i,
      "/",
      length(x),
      "\n"
    )
    Sys.sleep(0.5)
    setTxtProgressBar(pb, i / length(x))
    if (!file.exists(paste0(inputfiledir, x[i], ".csv"))) {
      stop(
        "\n \n haven't detected  ",
        x[i],
        "  file, Please check data/ directory if it has target files."
      )
    }
    assign(paste0("data", i), as.sparse(t(
      read.csv(
        paste0(inputfiledir, x[i], ".csv"),
        header = T,
        row.names = 1
      )
    )))
    message("\n ", x[i], "  Complete Input ")
  }
  close(pb)
  message("Input x Successful!! \n \n")
  message("Prepare to generate Seurat-list Object..\n")


  # Input to Seurat-list object, Use param cellmin,fmin
  message("Seurat-list Configuration : \n")
  message("          min.cells : ", cellmin, "\n")
  message("          min.features : ", fmin, "\n")
  message("        split meta-data to patient:", str_sub(meta, -1), "\n")
  message("        split meta-data to sample:", str_sub(meta, 1, 3), "\n")
  message("        Initialization viral counts as 0 \n")

  sclist <- list()
  pd <- txtProgressBar(style = 3)
  for (i in 1:length(x)) {
    cat(
      "\n Inputting No.",
      i,
      "element....",
      x[i],
      "........",
      i,
      "/",
      length(x),
      "\n"
    )
    Sys.sleep(0.5)
    setTxtProgressBar(pd, i / length(x))
    sclist[[i]] <-
      CreateSeuratObject(
        get(paste0("data", i)),
        project = x[i],
        min.cells = cellmin,
        min.features = fmin
      )
    sclist[[i]]@meta.data$orig.ident <- factor(meta[i])
    sclist[[i]]@meta.data$patient <- factor(str_sub(meta[i], -1))
    sclist[[i]]@meta.data$sample <- factor(str_sub(meta[i], 1, 3))
    sclist[[i]]@meta.data$virus <- 0
    sclist[[i]]@meta.data$hasvirus <- F
    message(
      "\n No.",
      i,
      " sclist element already compelete! ",
      "   Total: ",
      length(x),
      "\n"
    )
  }
  close(pd)
  message("Save inputted sclist as 'Sclist_input.Rdata' Rdata ")
  class(sclist) <- "sclist"
  save(sclist, file = "Sclist_input.Rdata")
  return(sclist)
}




#' Title Use to input Virus counts data to Vclist,{in: $sample viralname } {out: $(sample)sum}
#'
#' @param x give dataset that contail viruscount
#' @param inputfiledir input Virus counts data dir
#' @param viralname give a viralname to filter and input to vclist
#'
#' @return vclist
#' @export
#' @description defaultly,it will automatic detect Data in data/ directory.
#' The Virus-data type is like "CRR703.Exp.tsv" generate from Viral-track.
#' viralname is which virus you want to count for.
#' vclist<-Input_Vclist()  #it will automatic detect Data in data/ directory.
#' vclist<-Input_Vclist(x=viruscount,viralname="NC_003977") #you can use your own viruscount dataset and specific viralname.
#' @importFrom utils read.csv setTxtProgressBar str txtProgressBar

Input_Vclist <-
  function(x = viruscount,
           inputfiledir = "data/",
           viralname = "NC_003977") {
    if (missing(x)) {
      message(
        'Missing X, the Viruscount format is:x=c("CRR115711","CRR115712")    --- ',
        "Try Use default viruscount replace "
      )
    }
    if (!exists("viruscount")) {
      message("Missing viruscount, Try import **Exp.tsv in data/ ")
      if (checkmate::testFileExists(paste0(
        inputfiledir,
        list.files(path = inputfiledir, pattern = "_Exp.tsv")
      ))) {
        viruscount <- list.files(path = inputfiledir, pattern = "_Exp.tsv") %>% gsub(pattern = "_Exp.tsv", replacement = "")
        cat("\n data/ has **_Exp.tsv file:  ", viruscount, "\n")
      } else {
        stop("data/ hasn't detect **_Exp.tsv file: Interrupt import")
      }
    } else {
      message("!! Detect Viruscount variables in global variables, try to apply")
      cat("\n viruscount in global variables has content: ", viruscount)
      checkmate::assertFileExists(paste0(inputfiledir, viruscount, "_Exp.tsv"))
    }
    sum <- 0
    pb <- txtProgressBar(style = 3)
    vclist <- list()
    for (i in 1:length(x)) {
      cat(
        "\n Importing No.",
        i,
        " element............",
        i,
        "/",
        length(x),
        "\n"
      )
      Sys.sleep(0.5)
      setTxtProgressBar(pb, i / length(x), title = paste0("Importing ", x[1], "....."))
      exp <-
        read.csv(paste0(inputfiledir, x[i], "_Exp.tsv"),
          sep = "\t",
          row.names = 1
        )
      idx <- str_detect(rownames(exp), viralname)
      if (any(idx)) {
        assign(x[i], exp[idx, ])
        vclist[[x[i]]] <- colSums(get(x[i]))[colSums(get(x[i])) > 0]
        class(vclist[[i]]) <- "vcsingle"
        message(
          "\n Import  ",
          paste0(inputfiledir, x[i], "_Exp.tsv"),
          " Successfuly! "
        )
        sum <- sum + 1
      }
    }
    message(
      "\n \n#####    Finish import  Virus-Count files! ",
      " Summary: Total ",
      length(x),
      " files, ",
      "There are ",
      sum,
      "files have ",
      viralname,
      "!!"
    )
    close(pb)
    class(vclist) <- "vclist"
    message("Save imported Vclist as 'Vclist_input.Rdata' Rdata")
    save(vclist, file = "Vclist_input.Rdata")
    return(vclist)
  }










#' Title Input metadata
#' @return meta
#' @export
#' @description meta<-Input_meta()
#' @importFrom checkmate assertFileExists
#' @importFrom utils read.csv setTxtProgressBar str txtProgressBar
Input_meta <- function() {
  checkmate::assertFileExists("sample_info.tsv")
  meta <-
    read.csv(
      "sample_info.tsv",
      sep = "\t",
      header = F,
      col.names = c("RAW_data", "RAW_meta"),
      stringsAsFactors = F
    )
  message("\nImported sample_info.tsv")
  return(meta)
}




#' Title Check metadata exists, ifnot inputmetadata
#' @return meta
#' @description  meta<-Check_meta
#' @importFrom utils read.csv setTxtProgressBar str txtProgressBar
#'
Check_meta <- function() {
  if (!exists("metadata")) {
    message(
      "\nMissing metadata!! Try looking for Meta file import: The default META file is sample_info.tsv"
    )
    meta <- Input_meta()
    message("\nImport meta successfully")
  } else {
    message("\n!! Metadata variables were found in the Global Environment, trying to apply")
    if (colnames(metadata) != c("RAW_data", "RAW_meta")) {
      stop(
        "The metadata variable is invalid variable!Terminate the process, check the input parameters"
      )
    }
    cat("\n \nThe metadata has content:")
    str(metadata)
  }
  return(meta)
}

#' Title check fig/ dir, create it
#' @export
#' @description  Check_meta() ,create fig/
#'
Check_figdir <- function() {
  if (!dir.exists("fig")) {
    message("\nThe current directory does not exist fig/, create a new fig/ folder ...")
    dir.create("fig")
  } else {
    message("\n The current directory exist fig/ folder.")
  }
}





globalVariables(c("metadata", "vclist", "sclist", "scRNA", "ccRNA", "svclist","Testmeta","TestSclist","TestVclist"))
