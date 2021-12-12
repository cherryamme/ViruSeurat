


# loaded the package what todo 
.onAttach<-function(x) {
    suppressMessages({
        suppressWarnings({
            loaded=sapply(pkgs, require, character.only=TRUE)
        })
    })
    if(all(loaded)){
        Print_logo()
        Print_intro()
        Print_attach()
        Print_dir()
    }
}












#' @title Print_logo
#' @description Print_logo when attach

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export 
#' @rdname Print_logo
#' @seealso 
#'  
#' @import glue
Print_logo <- function() {
    logo<-c("\033[35m@\033[39m_   _  _     \033[32m*\033[39m         _____             \033[31mby Jesse\033[39m      _   ",
"| | | |(_)  \033[33m.\033[39m          /  ___|        \033[37mo\033[39m            \033[33m.\033[39m   | |  ",
"| | | | _  _ __  _   _ \ `--.   ___  _   _  _ __  __ _ | |_ ",
"| | | || || '__|| | | | `--. \\ / _ \\| | | || '__|/ _` || __|",
"\\ \\_/ /| || |   | |_| |/\\__/ /|  __/| |_| || |  | (_| || |_\033[35m@\033[39m",
" \\___/ |_||_|  \033[37mo\033[39m \\__,_|\\____/  \\___| \\__,_||_| \033[32m*\033[39m \\__,_| \\__|")
   glue::glue_col("{blue {logo}}")
}

#' @title Print_intro
#' @description print_intro when attach

#' @return intro
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export 
#' @rdname Print_intro

Print_intro <- function() {
    LOGO=c(bell="\ud83d\udd14",
       bulb="\ud83d\udca1",
       gift="\ud83c\udf81",
       bolt="\u26a1",
       star="\u2b50")
    cat(sample(LOGO,1),magenta(" ViruSeurat: Virus in Seurat Object with tracposon expression \n"))
}

#' @title Print_attach
#' @description Print_attach when attach

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export 
#' @rdname Print_attach
#' @seealso 
#'  
#' @import glue
Print_attach <- function() {
    pkgs=c("Seurat", "clustree", "ggplot2", "stringr", "patchwork")
    cat(blue(bold("Packages also been loaded:\n")))
    pkgs1= paste0(pkgs,collapse = " / ")
    glue::glue_col("{green {pkgs1}} \n")
}

#' @title Print_dir
#' @description Print_dir when attach

#' @return OUTPUT_DESCRIPTION
#' @details DETAILS
#' @examples 
#' \dontrun{
#' if(interactive()){
#'  #EXAMPLE1
#'  }
#' }
#' @export 
#' @rdname Print_dir
#' @seealso 
#'  
#' @import glue
Print_dir <- function() {
    dir=getwd()
    cat(blue(bold("Directory now in :")),yellow(dir),"\n")
    indir=paste0(list.dirs(recursive = F,full.names = F),collapse = " / ")
    glue::glue_col("{blue {bold Has contents dir :} {green {indir}} \n}")
}





