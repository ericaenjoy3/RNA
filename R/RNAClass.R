#####rlang .data prevents R CMD check from giving a NOTE about undefined global variables
#' @import data.table
#' @import limma
#' @import affy
#' @import ggplot2
#' @import scatterplot3d
#' @import RColorBrewer
#' @import grDevices
#' @import scales
#' @import GGally
#' @import GOtest
#' @import ape
#' @import pheatmap
#' @import ComplexHeatmap
#' @import circlize
#' @import multiplot
#' @import addgrids3d
#' @importFrom ModClusterR
###

tpm <- setClass(
  Class = "tpm",
  slots = c(tpm.value = "matrix", grps = "character"),
  validity = function(object){
    if (nrow(object@tpm.value) < 1) {
      return("Empty tpm.value matrix was given.")
    }
    if (length(grps) > 1) {
      return("Empty grps was given")
    }
    return(TRUE)
  }
)
