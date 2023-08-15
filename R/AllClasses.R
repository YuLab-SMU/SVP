#' @name SVPExperiment
#' @title The SVPExperiment class
#' @importFrom methods setClass
#' @importClassesFrom SingleCellExperiment SingleCellExperiment
#' @exportClass SVPExperiment
setClass('SVPExperiment', contains = 'SingleCellExperiment')
