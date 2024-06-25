\name{SparseM.ops}
\title{Basic Linear Algebra for Sparse Matrices}
\alias{Ops.matrix.csr}
\alias{Ops.matrix.diag.csr}
\alias{\%*\%-methods}
\alias{\%*\%,ANY,ANY-method}
\alias{\%*\%,matrix.csr,matrix.csr-method}
\alias{\%*\%,matrix.csr,matrix-method}
\alias{\%*\%,matrix.csr,numeric-method}
\alias{\%*\%,matrix,matrix.csr-method}
\alias{\%*\%,numeric,matrix.csr-method}
\alias{kronecker-methods}
%\alias{kronecker,ANY,ANY-method}
\alias{kronecker,matrix.csr,matrix.csr-method}
\alias{kronecker,matrix.csr,matrix-method}
\alias{kronecker,matrix.csr,numeric-method}
\alias{kronecker,matrix,matrix.csr-method}
\alias{kronecker,numeric,matrix.csr-method}
\alias{+,matrix.csr-method}
\alias{-,matrix.csr-method}
\alias{*,matrix.csr-method}
\alias{/,matrix.csr-method}
\alias{^,matrix.csr-method}
\alias{\%\%,matrix.csr-method}
\alias{\%/\%,matrix.csr-method}
\alias{>,matrix.csr-method}
\alias{>=,matrix.csr-method}
\alias{<,matrix.csr-method}
\alias{<=,matrix.csr-method}
\alias{==,matrix.csr-method}
\alias{!=,matrix.csr-method}
\alias{&,matrix.csr-method}
\alias{|,matrix.csr-method}
\alias{norm}
\alias{norm,ANY-method}
\alias{norm,matrix.csr,character-method}
\alias{norm,matrix.csr,missing-method}
\alias{det}% our own, calling _our_ determinant()
\alias{determinant}
\alias{determinant,matrix.csr.chol,missing-method}
\alias{determinant,matrix.csr.chol,logical-method}
\alias{determinant,matrix.csr,missing-method}
\alias{determinant,matrix.csr,logical-method}
%
\alias{t,ANY-method}
\alias{t,matrix.csr-method}
\alias{t,matrix.csc-method}
\alias{t,matrix.coo-method}
\alias{diag,ANY-method}
\alias{diag,matrix.csr-method}
\alias{diag<-,ANY-method}
\alias{diag<-,matrix.csr-method}
\alias{diag<-,matrix.diag.csr-method}
\alias{diff,matrix.csr-method}
\alias{diff<-,ANY-method}
\alias{diff<-,matrix.csr-method}
\alias{diag.assign,matrix.csr-method}
\alias{ncol,matrix.csr-method}
\alias{nrow,matrix.csr-method}
\alias{dim,ANY-method}
\alias{dim,matrix.csr-method}
\alias{dim,matrix.csc-method}
\alias{dim,matrix.ssr-method}
\alias{dim,matrix.ssc-method}
\alias{dim,matrix.coo-method}
\alias{rbind.matrix.csr}
\alias{cbind.matrix.csr}
\alias{[.matrix.csr}
\alias{[.matrix.diag.csr}
\alias{[<-.matrix.csr}
\alias{[<-.matrix.diag.csr}
\alias{[.matrix.coo}
\alias{[<-.matrix.coo}
\description{
  Basic linear algebra operations for sparse matrices,
  mostly of class \code{matrix.csr}.
}
\arguments{
\item{x}{matrix of class \code{matrix.csr}.}
\item{y}{matrix of class \code{matrix.csr} or a dense matrix or vector.}
\item{value}{replacement values.}
\item{i,j}{vectors of elements to extract or replace.}
\item{nrow}{optional number of rows for the result.}
\item{lag}{an integer indicating which lag to use.}
\item{differences}{an integer indicating the order of the difference.}
}
\details{Linear algebra operations for matrices of class
\code{matrix.csr} are designed to behave exactly as for
regular matrices.  In particular, matrix multiplication, kronecker
product, addition,
subtraction and various logical operations should work as with the conventional
dense form of matrix storage, as does indexing, rbind, cbind, and diagonal
assignment and extraction.  The method diag may be used to extract the
diagonal of a \code{matrix.csr} object, to create a sparse diagonal see
\code{SparseM.ontology}.

The function \code{\link[base]{determinant}} computes the (log) determinant,
of the argument, returning a \code{"det"} object as the base function.
This is typically \emph{preferred} over using the function \code{det()}
which is a simple wrapper for \code{determinant()}, in a way it will work
for our sparse matrices, as well.
\code{determinant()} computes the determinant of the argument
matrix.  If the matrix is of class \code{matrix.csr} then it must
be symmetric, or an error will be returned.  If the matrix is of
class \code{matrix.csr.chol} then the (pre-computed) determinant of the Cholesky
factor is returned, i.e., the product of the diagonal elements.

The function \code{norm}, i.e. \code{norm(x, type)}, by default computes
the \dQuote{sup} (or \code{"M"}aximum norm, i.e., the maximum of the matrix
elements.  Optionally, this \code{type = "sup"} (equivalently, \code{type = "M"}) norm can
be replaced by the Hilbert-Schmidt, \code{type = "HS"} or equivalently, \code{type = "F"}
norm, or the \code{type = "l1"}, norm.
Note that for historical reasons, the default \code{type} differs from \R's
own \code{\link[base]{norm}()}, see the examples using \code{B}, below.
The \code{"sup" === "M"} and \code{"HS" === "F"} equivalences have been
introduced in \pkg{SparseM} version 1.84.
}
\references{
  Koenker, R and Ng, P. (2002).
  SparseM: A Sparse Matrix Package for \R,\cr
  \url{http://www.econ.uiuc.edu/~roger/research/home.html}
}
\seealso{
  \code{\link{slm}} for sparse linear model fitting.
  \code{\link{SparseM.ontology}} for coercion and other class relations
  involving our sparse matrix classes.
}
\examples{
n1 <- 10
n2 <- 10
p <- 6
y <- rnorm(n1)
a <- round(rnorm(n1*p), 2)
a[abs(a) < 0.5] <- 0
A <- matrix(a, n1,p)
A.csr <- as.matrix.csr(A)

B <- matrix(c(1.5, 0,  0,  0, -1.4,    0,    0,   0,   0, -1.4,
              2,   0, -1,  0,  0,      2.1, -1.9, 1.4, 0,  0,
              0,-2.3,  0,  0, -1.9,    0,    0,   0,   0, -1.4,
              0,   0,  0,  0,  0,     -3,    0,   1.3, 0,  1.1,
              0,   0,  0,  0,  2,      0,    0,   0,  -1,  0,
              0,   0, -1.6,0,  0,      0,    0,   0,  -1.7,0),
            10L, 6L)
rcond(B) # 0.21 .. i.e., quite well conditioned
B.csr <- as.matrix.csr(B)
B.csr

## norm() : different 'type' for base R and  SparseM:
(nR <- vapply(c("1", "I", "F", "M", "2"), norm, 0, x = B))
##        1        I        F        M        2
## 8.400000 5.300000 7.372923 3.000000 4.464801
(nSpM <- vapply(c("sup","M",  "HS","F",  "l1"), norm, 0, x = B.csr))
##      sup         M        HS         F        l1
## 3.000000  3.000000  7.372923  7.372923 30.000000
stopifnot(all.equal(unname(nSpM[c("M", "F")]),
                    unname(nR  [c("M", "F" )]), tolerance = 1e-14))

# matrix transposition and multiplication
BtB <- crossprod(B) # == t(B) \%*\% B  {more efficiently}

A.csr \%*\% t(B.csr)
BtBs <- t(B.csr) \%*\% B.csr
BtBs
stopifnot(all.equal(    BtB,  as.matrix(BtBs),  tolerance = 1e-14),
          all.equal(det(BtB), print(det(BtBs)), tolerance = 1e-14))

# matrix  o  vector
stopifnot(all.equal(y \%*\%  A , y   \%*\% A.csr) ,
          all.equal(A \%*\% 1:6, A.csr \%*\% 1:6)
)

# kronecker product - via kronecker() methods:
A.csr \%x\% matrix(1:4,2,2)
}

\keyword{algebra}
