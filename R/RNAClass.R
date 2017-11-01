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
#' @import ModClusterR
#' @import methods
#' @import ggrepel
#' @importFrom utils read.table write.table combn
#' @importFrom stats cor model.matrix as.dist hclust cutree quantile lm var prcomp kmeans
#' @importFrom graphics plot text legend par
#' @importFrom cluster pam
###

tpm <- setClass(
  Class = "tpm",
  slots = c(tpm.value = "matrix", grps = "character"),
  validity = function(object){
    if (nrow(object@tpm.value) < 1) {
      return("Empty tpm.value matrix was given.")
    }
    if (length(object@grps) < 1) {
      return("Empty grps was given")
    }
    return(TRUE)
  }
)
