\name{matrix.csc-class}
\docType{class}
\alias{matrix.csc-class}
\title{Class "matrix.csc"}
\description{A new class for sparse matrices stored in compressed sparse column format}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("matrix.csc", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{ra}:}{Object of class \code{"numeric", from class "matrix.csr"} a real array of nnz elements containing the non-zero
        elements of A, stored in row order. Thus, if i<j, all elements
        of row i precede elements from row j. The order of elements
        within the rows is immaterial.}
    \item{\code{ja}:}{Object of class \code{"numeric", from class "matrix.csr"} an integer array of nnz elements containing the column
        indices of the elements stored in `ra'.}
    \item{\code{ia}:}{Object of class \code{"numeric", from class "matrix.csr"} an integer array of n+1 elements containing pointers to
        the beginning of each row in the arrays `ra' and `ja'. Thus
        `ia[i]' indicates the position in the arrays `ra' and `ja'
        where the ith row begins. The last, (n+1)st, element of `ia'
        indicates where the n+1 row would start, if it existed.
 }
    \item{\code{dimension}:}{Object of class \code{"numeric", from class "matrix.csr"} dimension of the matrix }
  }
}
\section{Extends}{
Class \code{"matrix.csr"}, directly.
}
\section{Methods}{
  \describe{
    \item{as.matrix.csr}{\code{signature(x = "matrix.csc")}: ... }
    \item{as.matrix.ssc}{\code{signature(x = "matrix.csc")}: ... }
    \item{as.matrix.ssr}{\code{signature(x = "matrix.csc")}: ... }
    \item{as.matrix}{\code{signature(x = "matrix.csc")}: ... }
    \item{chol}{\code{signature(x = "matrix.csc")}: ... }
    \item{dim}{\code{signature(x = "matrix.csc")}: ... }
    \item{t}{\code{signature(x = "matrix.csc")}: ... }
  }
}

\seealso{ \code{\link{matrix.csr-class}}}

\keyword{classes}

