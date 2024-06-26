% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-lisaRes.R
\name{LISAResult}
\alias{LISAResult}
\title{LISAResult}
\usage{
LISAResult(x, type = NULL, features = NULL, ...)
}
\arguments{
\item{x}{object \linkS4class{SingleCellExperiment}.}

\item{type}{character, the name of \code{method} parameter of \code{runLISA()} combining
\code{.SVP},so it can be one of \code{localG.SVP} or \code{localmoran.SVP},
default is NULL.}

\item{features}{character or index which have been specified in \code{features} of
\verb{code\{runLISA()\}} and \code{action='add'}, default is NULL.}

\item{...}{additional parameter, meaningless now.}
}
\value{
a data.frame or SimpleList.
}
\description{
Extracting the result of \code{runLISA()}
}
\examples{
data(hpda_spe_cell_dec)
hpda_spe_cell_dec <- hpda_spe_cell_dec |> 
                     runLISA(features = 'Cancer clone A', 
                             assay.type = 'affi.score', 
                             method = 'localG',
                             action = 'add'
                     ) 
hpda_spe_cell_dec <- hpda_spe_cell_dec |>
                     runLISA(features = 'Cancer clone A', 
                             assay.type = 'affi.score', 
                             method = 'localmoran',
                             action = 'add'
                     )
local.G <- LISAResult(hpda_spe_cell_dec, 
                      type='localG.SVP', features='Cancer clone A'
           )
localmoran <- LISAResult(hpda_spe_cell_dec, 
                         type = 'logcalmoran.SVP', 
                         features = 'Cancer clone A'
           )
hpda_spe_cell_dec |> LISAResult() |> head()
}
