#' @title The Gene List of Cancer Single-cell State Atlas (CancerSEA)
#' 
#' @description CancerSEA is the first dedicated database that aims to 
#' comprehensively decode distinct functional states of cancer cells at 
#' single-cell resolution.
#' CancerSEASymbol is a gene symbol list, and CancerSEAEnsemble is a Ensemble 
#' gene list, they are a list contained gene signature names collected in the database.
#'
#' @name data_CancerSEA
#' @aliases CancerSEASymbol
#' @format list
#' a gene symbol list with gene signature names collected in CancerSEA:
#' \describe{
#'    \item{Angiogenesis}{Angiogenesis ensures that cancer cells receive continuous supplies 
#'       of oxygen and other nutrients.}
#'    \item{Apoptosis}{The inactivation of apoptosis in cancer cells lead to the persistence 
#'       of such grossly abnormal cells in the tissues.}
#'    \item{Cell Cycle}{Cell cycle,a critical process to ensure correct cell division,lies 
#'       at the heart of cancer.}
#'    \item{Differentiation}{The degree of cell differentiation can be used to measure the 
#'       progress of cancer,and dedifferentiated cells can lead to the formation of cancer.}
#'    \item{DNA damage}{DNA damage is an alteration in the chemical structure of DNA, and 
#'       un-repaired DNA damages accumulate in replicating cells possibly contribute to 
#'       progression to cancer.}
#'    \item{DNA repair}{DNA repair plays a fundamental role in the maintenance of genomic 
#'       integrity,it's deficits may lead to carcinogenesis.}
#'    \item{EMT}{EMT has been indicated to be involved in the initiation of metastasis in 
#'       cancer progression and in acquiring drug resistance.}
#'    \item{Hypoxia}{Tumor-hypoxia contributes to cell mobility,metastasis and therapy resistance.}
#'    \item{Inflammation}{Chronic inflammation can cause about 15% to 25% of human cancers.}
#'    \item{Invasion}{Invasion is a critical carcinogenic event in which cancer cells escape 
#'       from their primary sites and spread to blood or lymphatic vessels.}
#'    \item{Metastasis}{Metastasis promotes the malignant transformation of cancer and causes 
#'       most cancer deaths.}
#'    \item{Proliferation}{Proliferation,as one of the cancer hallmarks,is responsible for tumor 
#'       progression.}
#'    \item{Quiescence}{Quiescent cancer cells are resistant to chemotherapy.}
#'    \item{Stemness}{Cancer cells with high stemness fuel the growth of cancer.}
#'  }
#' @docType data
#' @return a list object
#' @keywords data
#' @references Yuan, H., Yan, M., Zhang, G., Liu, W., Deng, C., Liao, G., Xu, L., Luo, T., Yan, H., 
#' Long, Z., Shi, A., Zhao, T., Xiao, Y., & Li, X. (2019). CancerSEA: a cancer single-cell state atlas. 
#' Nucleic acids research, 47(D1), D900–D908. https://doi.org/10.1093/nar/gky939
#' @source
#' \url{http://biocc.hrbmu.edu.cn/CancerSEA/goDownload}
#' @examples
#' data(CancerSEASymbol)
NA

#' @rdname data_CancerSEA
#' @name data_CacerSEA
#' @aliases CancerSEAEnsemble
#' @format list
#' @docType data
#' @keywords data
#' @examples
#' data(CancerSEAEnsemble)
NA

#' @title A gene set identifies senescent cells and predicts senescence-associated 
#' pathways across tissues
#'
#' @description
#' SenMayoSymbol is a gene symbol list that can be used to identify senescent cells and predicts 
#' senescence-associated pathways across tissues
#'
#' @name data_SenMayo
#' @aliases SenMayoSymbol
#' @format list
#' @docType data
#' @keywords data
#' @return a list object
#' @references Saul, D., Kosinsky, R.L., Atkinson, E.J. et al. A new gene set identifies senescent 
#' cells and predicts senescence-associated pathways across tissues. Nat Commun 13, 4827 (2022). 
#' https://doi.org/10.1038/s41467-022-32552-1
#' @source
#' \url{https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-022-32552-1/MediaObjects/41467_2022_32552_MOESM4_ESM.xlsx}
#' @examples
#' data(SenMayoSymbol)
NA

#' @title a subset data of pbmck3 from SeuratData
#'
#' @description
#' a small SingleCellExperiment data set from pbmck3 which
#' contains 1304 genes and 800 cells (extract randomly)
#' 
#' @name data_sceSubPbmc
#' @format S4 class:SingleCellExperiment
#' @aliases sceSubPbmc
#' @docType data
#' @keywords data
#' @return a \linkS4class{SingleCellExperiment} object
#' @examples
#' data(sceSubPbmc)
NA

#' @title the Cell Cycle gene set
#' 
#' @description
#' the S and G2M gene list are from the Seurat which refer to this article (doi:10.1126/science.aad050),
#' the G1 gene list is from the G1_PHASE of Human Gene Set in MSigDB, but remove the duplicated records 
#' with S and G2M gene list.
#' 
#' @name CellCycle.Hs
#' @aliases data_CellCycle.Hs
#' @format list
#' @docType data
#' @keywords data
#' @return a list object
#' @examples
#' data(CellCycle.Hs)
NA


#' @title an example of result of runSGSA by extracting with gsvaExp
#'
#' @description
#' The result of runSGSA with HPDA A sample from (doi:10.1038/s41587-019-0392-8) 
#'
#' @name data_hpda_spe_cell_dec
#' @aliases hpda_spe_cell_dec
#' @format S4 class:SpatialExperiment
#' @docType data
#' @keywords data
#' @return a \linkS4class{SpatialExperiment} object
#' @examples
#' data(hpda_spe_cell_dec)
NA
