\name{matrix.csr-class}
\docType{class}
\alias{matrix.csr-class}
\title{Class "matrix.csr"}
\description{A new class for sparse matrices stored in compressed sparse row format}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("matrix.csr", ...)}.
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

\section{Methods}{
  \describe{
    \item{\%*\%}{\code{signature(x = "matrix.csr", y = "matrix.csr")}: ... }
    \item{\%*\%}{\code{signature(x = "matrix.csr", y = "numeric")}: ... }
    \item{as.matrix.csc}{\code{signature(x = "matrix.csr")}: ... }
    \item{as.matrix.ssc}{\code{signature(x = "matrix.csr")}: ... }
    \item{as.matrix.ssr}{\code{signature(x = "matrix.csr")}: ... }
    \item{as.matrix}{\code{signature(x = "matrix.csr")}: ... }
    \item{chol}{\code{signature(x = "matrix.csr")}: ... }
    \item{diag}{\code{signature(x = "matrix.csr")}: ... }
    \item{diag<-}{\code{signature(x = "matrix.csr")}: ... }
    \item{dim}{\code{signature(x = "matrix.csr")}: ... }
    \item{image}{\code{signature(x = "matrix.csr")}: ... }
    \item{solve}{\code{signature(a = "matrix.csr")}: ... }
    \item{t}{\code{signature(x = "matrix.csr")}: ... }
  }
}
\seealso{\code{\link{matrix.csc-class}}}
\keyword{classes}

