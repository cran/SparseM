\name{matrix.csr.chol-class}
\docType{class}
\alias{matrix.csr.chol-class}
\title{Class "matrix.csr.chol"}
\description{A class of objects returned from Ng and Peyton's (1993) block
sparse Cholesky algorithm}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("matrix.csr.chol", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{nrow}:}{Object of class \code{"numeric"} number of rows in the linear system of equations }
    \item{\code{nnzlindx}:}{Object of class \code{"numeric"} number of non-zero elements in lindx}
    \item{\code{nsuper}:}{Object of class \code{"numeric"} number of supernodes}
    \item{\code{lindx}:}{Object of class \code{"numeric"} vector of integer
containing, in column major order, the row subscripts of the non-zero entries 
in the Cholesky factor in a compressed storage format}
    \item{\code{xlindx}:}{Object of class \code{"numeric"} vector of integer of pointers for lindx}
    \item{\code{nnzl}:}{Object of class \code{"numeric"} number of non-zero
entries, including the diagonal entries, of the Cholesky factor stored in lnz}
    \item{\code{lnz}:}{Object of class \code{"numeric"} contains the entries of the Cholesky factor }
    \item{\code{xlnz}:}{Object of class \code{"numeric"} column pointer for the Cholesky factor stored in lnz}
    \item{\code{invp}:}{Object of class \code{"numeric"} vector of integer of inverse permutation vector}
    \item{\code{perm}:}{Object of class \code{"numeric"} vector of integer of permutation vector }
    \item{\code{xsuper}:}{Object of class \code{"numeric"} array containing the supernode partioning}
    \item{\code{ierr}:}{Object of class \code{"numeric"} error flag }
    \item{\code{time}:}{Object of class \code{"numeric"} execution time }
  }
}

\section{Methods}{
  \describe{
    \item{backsolve}{\code{signature(r = "matrix.csr.chol")}: ... }
  }
}
\seealso{\code{\link{chol}}, \code{\link{backsolve}}}
\keyword{classes}

