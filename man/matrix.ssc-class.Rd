\name{matrix.ssc-class}
\docType{class}
\alias{matrix.ssc-class}
\title{Class "matrix.ssc"}
\description{A new class for sparse matrices stored in symmetric sparse column format}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("matrix.ssc", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{ra}:}{Object of class \code{numeric}, a real array of nnz elements containing the non-zero
        elements of the lower triangular part of A, stored in column order. Thus, if i<j, all elements
        of column i precede elements from column j. The order of elements
        within the column is immaterial.}
    \item{\code{ja}:}{Object of class \code{integer}, an integer array of nnz elements containing the row
        indices of the elements stored in `ra'.}
    \item{\code{ia}:}{Object of class \code{integer}, an integer array of n+1 elements containing pointers to
        the beginning of each column in the arrays `ra' and `ja'. Thus
        `ia[i]' indicates the position in the arrays `ra' and `ja'
        where the ith column begins. The last, (n+1)st, element of `ia'
        indicates where the n+1 column would start, if it existed.
 }
    \item{\code{dimension}:}{Object of class \code{integer}, dimension of the matrix }
  }
}
\section{Methods}{
  \describe{
    \item{as.matrix.csc}{\code{signature(x = "matrix.ssc")}: ... }
    \item{as.matrix.csr}{\code{signature(x = "matrix.ssc")}: ... }
    \item{as.matrix.ssr}{\code{signature(x = "matrix.ssc")}: ... }
    \item{as.matrix}{\code{signature(x = "matrix.ssc")}: ... }
    \item{dim}{\code{signature(x = "matrix.ssc")}: ... }
  }
}
\seealso{\code{\link{matrix.csr-class}}}
\keyword{classes}

