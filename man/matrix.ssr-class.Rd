\name{matrix.ssr-class}
\docType{class}
\alias{matrix.ssr-class}
\title{Class "matrix.ssr"}
\description{A new class for sparse matrices stored in symmetric sparse row format}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("matrix.ssr", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{ra}:}{Object of class \code{numeric}, a real array of nnz elements containing the non-zero
        elements of the lower triangular part of A, stored in row order. Thus, if i<j, all elements
        of row i precede elements from row j. The order of elements
        within the rows is immaterial.}
    \item{\code{ja}:}{Object of class \code{integer}, an integer array of nnz elements containing the column
        indices of the elements stored in `ra'.}
    \item{\code{ia}:}{Object of class \code{integer}, an integer array of n+1 elements containing pointers to
        the beginning of each row in the arrays `ra' and `ja'. Thus
        `ia[i]' indicates the position in the arrays `ra' and `ja'
        where the ith row begins. The last, (n+1)st, element of `ia'
        indicates where the n+1 row would start, if it existed.
 }
    \item{\code{dimension}:}{Object of class \code{integer}, dimension of the matrix }
  }
}
\section{Methods}{
  \describe{
    \item{as.matrix.csc}{\code{signature(x = "matrix.ssr")}: ... }
    \item{as.matrix.csr}{\code{signature(x = "matrix.ssr")}: ... }
    \item{as.matrix.ssr}{\code{signature(x = "matrix.ssr")}: ... }
    \item{as.matrix}{\code{signature(x = "matrix.ssr")}: ... }
    \item{dim}{\code{signature(x = "matrix.ssr")}: ... }
  }
}
\seealso{\code{\link{matrix.csr-class}}}
\keyword{classes}

