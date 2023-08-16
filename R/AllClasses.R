#' @title The SVPExperiment class
#' @docType class
#' @importFrom methods setClass
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @exportClass SVPExperiment
setClass('SVPExperiment', contains = 'SingleCellExperiment')

setClass('SCEByColumn', slot = c(sce = 'SingleCellExperiment'))
