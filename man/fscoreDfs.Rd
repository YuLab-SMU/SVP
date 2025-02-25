% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods-fscore.R
\docType{methods}
\name{fscoreDfs}
\alias{fscoreDfs}
\alias{fscoreDf}
\alias{fscoreDfNames}
\alias{fscoreDf,SingleCellExperiment,missing-method}
\alias{fscoreDf,SingleCellExperiment,numeric-method}
\alias{fscoreDf,SingleCellExperiment,character-method}
\alias{fscoreDfs,SingleCellExperiment-method}
\alias{fscoreDfNames,SingleCellExperiment-method}
\alias{fscoreDf<-}
\alias{fscoreDfs<-}
\alias{fscoreDfNames<-}
\alias{fscoreDf<-,SingleCellExperiment,missing-method}
\alias{fscoreDf<-,SingleCellExperiment,numeric-method}
\alias{fscoreDf<-,SingleCellExperiment,character-method}
\alias{fscoreDfs<-,SingleCellExperiment-method}
\alias{fscoreDfNames<-,SingleCellExperiment,character-method}
\title{features score matrix extract method}
\value{
see \code{Getter} and \code{setter}
}
\description{
In some experiment, to calculated the contribution value of original features
(such as genes) in the new features (gene sets), if the result is stored with
the original object, which will simplify book-keeping in long workflows and
ensure that samples remain synchronised.
}
\section{Getters}{

In the following examples, \code{x} is a \linkS4class{SingleCellExperiment} object.
\describe{
\item{\code{fscoreDf(x, type)}:}{
Retrieves a \linkS4class{DataFrame} containing the new features (gene sets) (rows)
for the specified \code{type}.
\code{type} should either be a string specifying the name of the features scores matrix
in \code{x} to retrieve, or a numeric scalar specifying the index of the desired matrix,
defaulting to the first matrix is missing.
}

\item{\code{fscoreDfNames(x)}:}{
Retures a character vector containing the names of all features scores DataFrame Lists in
\code{x}. This is guaranteed to be of the same length as the number of results.
}

\item{\code{fscoreDfs(x)}:}{
Returns a named \linkS4class{List} of matrices containing one or more \linkS4class{DataFrame} objects.
Each object is guaranteed to have the same number of rows, in a 1:1 correspondence to those in \code{x}.
}
}
}

\section{Single-object setter}{

\code{fscoreDf(x, type) <- value} will add or replace an features scores matrix in a
\linkS4class{SingleCellExperiment} object \code{x}.
The value of \code{type} determines how the result is added or replaced:
\itemize{
\item If \code{type} is missing, \code{value} is assigned to the first result.
If the result already exists, its name is preserved; otherwise it is given a default name \code{"unnamed.fscore1"}.
\item If \code{type} is a numeric scalar, it must be within the range of existing results, and \code{value} will
be assigned to the result at that index.
\item If \code{type} is a string and a result exists with this name, \code{value} is assigned to to that result.
Otherwise a new result with this name is append to the existing list of results.
}
}

\section{Other setter}{

\describe{
\item{\code{fscoreDfs(x) <- value}:}{
Replaces all features score matrixs in \code{x} with those in \code{value}.
The latter should be a list-like object containing any number of \linkS4class{DataFrame} objects
with number of row equal to \code{nrow(x)}.

If \code{value} is named, those names will be used to name the features score matrixs in \code{x}.
Otherwise, unnamed results are assigned default names prefixed with \code{"unnamed.fscore"}.

If \code{value} is \code{NULL}, all features score matrixs in \code{x} are removed.
}

\item{\code{fscoreDfNames(x) <- value}:}{
Replaces all names for features score matrixs in \code{x} with a character vector \code{value}.
This should be of length equal to the number of results currently in \code{x}.
}
}
}

\examples{
# Using the class example
example(SVPExperiment, echo = FALSE)
dim(counts(svpe))
rownames(svpe) <- paste0("gene", seq(nrow(svpe)))
colnames(svpe) <- paste0("cell", seq(ncol(svpe)))
# Mocking up some GSVA Experiments
sce1 <- SingleCellExperiment(matrix(rpois(1000, 5), ncol=ncol(svpe)))
rownames(sce1) <- paste0("GO:",seq(nrow(sce1)))
colnames(sce1) <- colnames(svpe)
sce2 <- SingleCellExperiment(matrix(rpois(1000, 5), ncol=ncol(svpe)))
rownames(sce2) <- paste0("KEGG:", seq(nrow(sce2)))
colnames(sce2) <- colnames(svpe)

# Mocking up some relationship score between new feature and gene
fscore1 <- lapply(seq(nrow(sce1)), function(i) abs(rnorm(5, 0.5)) |> 
             setNames(sample(rownames(svpe),5))) |> 
             List() |>  
             DataFrame() |> setNames("rwr_score") 
rownames(fscore1) <- rownames(sce1)
fscore2 <- lapply(seq(nrow(sce2)), function(i) abs(rnorm(5, 0.8)) |> 
             setNames(sample(rownames(svpe),5))) |> 
             List() |> 
             DataFrame() |> setNames("hyper_test")
# Setting the score
fscoreDfs(sce1) <- list()
fscoreDfs(sce2) <- list()
fscoreDf(sce1, "rwr_score") <- fscore1
fscoreDf(sce2, "hyper_test") <-fscore2

# Setting the GSVA Experiments
gsvaExp(svpe, "GO1") <- sce1
gsvaExp(svpe, "KEGG1") <- sce2

# Getting the GSVA Experiment data
fscoreDf(gsvaExp(svpe), "rwr_score")
fscoreDf(gsvaExp(svpe, 'KEGG1'), "hyper_test")
fscoreDf(gsvaExp(svpe, 'KEGG1'), 1)
fscoreDfNames(gsvaExp(svpe))
fscoreDfs(gsvaExp(svpe))

# Setting the names of features score DataFrame
fscoreDfNames(gsvaExp(svpe, withColData=FALSE)) <- "rwr.score"
fscoreDfNames(gsvaExp(svpe, withColData=FALSE))[1] <- "Test"
}
