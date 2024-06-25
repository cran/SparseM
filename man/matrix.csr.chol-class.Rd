\name{matrix.csr.chol-class}
\docType{class}
\alias{matrix.csr.chol-class}
\title{Class "matrix.csr.chol" (Block Sparse Cholesky Decomposition)}
\description{
  A class of objects returned from Ng and Peyton's (1993) block
  sparse Cholesky algorithm.
}
\section{Objects from the Class}{
  Objects may be created by calls of the form \code{new("matrix.csr.chol",
    ...)}, but typically result from
  \code{\link[=SparseM.solve]{chol}(<matrix.csr>)}.
}
\section{Slots}{
  \describe{
    \item{\code{nrow}:}{an \code{integer}, the number of rows of the
      original matrix, or in the linear system of equations.}
    \item{\code{nnzlindx}:}{Object of class \code{numeric}, number of non-zero elements in \code{lindx}}
    \item{\code{nsuper}:}{an \code{integer}, the number of supernodes of the decomposition}
    \item{\code{lindx}:}{Object of class \code{integer}, vector of integer
      containing, in column major order, the row subscripts of the non-zero entries
      in the Cholesky factor in a compressed storage format}
    \item{\code{xlindx}:}{Object of class \code{integer}, vector of integer of pointers for \code{lindx}}
    \item{\code{nnzl}:}{of class \code{"numeric"}, an integer, the number of non-zero
      entries, including the diagonal entries, of the Cholesky factor stored in \code{lnz}}
    \item{\code{lnz}:}{a \code{numeric} vector of the entries of the Cholesky factor}
    \item{\code{xlnz}:}{an \code{integer} vector, the column pointers for the Cholesky factor stored in \code{lnz}}
    \item{\code{invp}:}{inverse permutation vector, \code{integer}}
    \item{\code{perm}:}{permutation vector, \code{integer}}
    \item{\code{xsuper}:}{Object of class \code{integer}, array containing the supernode partioning}
    \item{\code{det}:}{\code{\link{numeric}}, the determinant of the Cholesky factor}
    \item{\code{log.det}:}{\code{\link{numeric}}, the log determinant of the Cholesky factor}
    \item{\code{ierr}:}{an \code{integer}, the error flag (from Fortran's \file{src/chol.f})}
    \item{\code{time}:}{\code{\link{numeric}}, unused (always \code{0.}) currently.}
  }
}
\details{
  Note that the \code{perm} and notably \code{invp} maybe important to back
  permute rows and columns of the decompositions, see the Examples, and our
  \code{\link[=SparseM.solve]{chol}} help page.
}

\section{Methods}{
  \describe{
    \item{as.matrix.csr}{\code{signature(x = "matrix.csr.chol",
	upper.tri=TRUE)}: to get the sparse (\code{"matrix.csr"}) upper
      triangular matrix corresponding to the Cholesky decomposition.}
    \item{backsolve}{\code{signature(r = "matrix.csr.chol")}: for computing
      \eqn{R^{-1} b} when the Cholesky decomposition is \eqn{A = R'R}.}
  }
}
\seealso{Base \R's \code{\link{chol}} and \pkg{SparseM}'s
  \code{\link[=SparseM.solve]{chol}}, notably for examples;
  \code{\link{backsolve}}
}
\examples{
x5g <- new("matrix.csr",
          ra = c(300, 130, 5, 130, 330,
                 125, 10, 5, 125, 200, 70,
                 10, 70, 121.5, 1e30),
          ja = c(1:3, 1:4, 1:4, 2:5),
          ia = c(1L, 4L, 8L, 12L, 15L, 16L),
          dimension = c(5L, 5L))
(m5g <- as.matrix(x5g)) # yes, is symmetric, and positive definite:
eigen(m5g, only.values=TRUE)$values  # all positive (but close to singular)
ch5g <- chol(x5g)
str(ch5g) # --> the slots of the "matrix.csr.chol" class
mch5g <- as.matrix.csr(ch5g)
print.table(as.matrix(mch5g), zero.print=".") # indeed upper triagonal w/ positive diagonal

## x5 has even more extreme entry at [5,5]:
x5 <- x5g; x5[5,5] <- 2.9e32
m5 <- as.matrix(x5)
(c5 <- chol(m5))# still fine, w/ [5,5] entry = 1.7e16 and other diag.entries in (9.56, 17.32)
ch5 <- chol(x5) # --> warning  "Replaced 3 tiny diagonal entries by 'Large'"
                # gave error for a while
(mmc5 <- as.matrix(as.matrix.csr(ch5)))
        # yes, these replacements were extreme, and the result is "strange'
## Solve the problem (here) specifying non-default  singularity-tuning par 'tiny':
ch5. <- chol(x5, tiny = 1e-33)
(mmc5. <- as.matrix(as.matrix.csr(ch5.))) # looks much better.
## Indeed: R'R  back-permuted *is* the original matrix x5, here m5:
(RtR <- crossprod(mmc5.)[ch5.@invp, ch5.@invp])
          all.equal(m5, RtR, tolerance = 2^-52)
stopifnot(all.equal(m5, RtR, tolerance = 1e-14)) # on F38 Linux, only need tol = 1.25e-16
}
\keyword{classes}

