\name{matrix.coo-class}
\docType{class}
\alias{matrix.coo-class}
\title{Class "matrix.coo"}
\description{A new class for sparse matrices stored in coordinate format}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("matrix.coo", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{ra}:}{Object of class \code{numeric}, a real array of nnz elements containing the non-zero
        elements of A.} 
    \item{\code{ja}:}{Object of class \code{integer}, an integer array of nnz elements containing the column
        indices of the elements stored in `ra'.}
    \item{\code{ia}:}{Object of class \code{integer}, an integer array of nnz elements containing the row
        indices of the elements stored in `ra'.}
    \item{\code{dimension}:}{Object of class \code{integer}, dimension of the matrix }
  }
}
\section{Methods}{
  \describe{
    \item{as.matrix.csr}{\code{signature(x = "matrix.coo")}: ... }
    \item{as.matrix}{\code{signature(x = "matrix.coo")}: ... }
    \item{dim}{\code{signature(x = "matrix.coo")}: ... }
  }
}

\seealso{ \code{\link{matrix.csr-class}}}

\keyword{classes}

