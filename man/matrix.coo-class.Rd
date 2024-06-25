\name{matrix.coo-class}
\title{Class "matrix.coo" -- Sparse Matrices in [Coo]rdinate Format}
\docType{class}
\alias{matrix.coo-class}
\description{
  Class for sparse matrices stored in coordinate aka
  \dQuote{triplet} format, storing for each non-zero entry \code{x[i,j]}
  the triplet \code{(i,j, x[i,j])}, in slots \code{(ia, ja, ra)}.
}
\section{Objects from the Class}{
  Objects can be created by calls of the form \code{new("matrix.coo", ...)},
  but typically rather by \code{\link{as.matrix.coo}()}.
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
    \item{as.matrix.coo}{\code{signature(x = "matrix.coo")}: ... }
    \item{as.matrix.csr}{\code{signature(x = "matrix.coo")}: ... }
    \item{as.matrix}{\code{signature(x = "matrix.coo")}: ... }
    \item{dim}{\code{signature(x = "matrix.coo")}: ... }
  }
}

\seealso{ \code{\link{matrix.csr-class}}}

\examples{
 try(  new("matrix.coo") ) # fails currently {FIXME!}     # the 1x1 matrix [0]

 ## correponds to base  matrix()
 mcoo <- new("matrix.coo", ra=NA_real_, ia = 1L, ja = 1L, dimension = c(1L, 1L))
 mcoo # currently *does* print but wrongly:  as.matrix.csr(<matrix.coo>) fails to keep NA !!
 co2 <- new("matrix.coo", ra = c(-Inf, -2, 3, Inf),
            ia = c(1L,1:3), ja = 2L*(1:4), dimension = c(7L, 12L))
 co2 # works fine (as has no NA)

 ## Sparse Diagonal (from "numeric"):
 as(1:5, "matrix.diag.csr") # a sparse version of  diag(1:5)
}

\keyword{classes}

