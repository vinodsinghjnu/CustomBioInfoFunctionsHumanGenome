#' @title  create a directory and delete if it already exists
#'
#' @description create a directory and delete if it already exists
#' @param dir Name of the dir to be created
#' @return dir created in the current diretory
#' @export
#' @examples
#' createDir_delIfExists(dir='testDir')
#'
#'
#'
createDir_delIfExists=function(dir){
  if (dir.exists(dir)) {
    unlink(dir,recursive = TRUE)
    cat("Old Directory has been deleted")
    dir.create(dir)
    cat("New Directory has been created")
  }else{
    dir.create(dir)
    cat("New Directory has been created")
  }
}
