\name{SparseM.ops}
\alias{%*%}
\alias{%*%.default}
\alias{Ops.matrix.csr}
\alias{SparseM.ops}
\alias{\%*\%.matrix.csr}
\alias{+.matrix.csr}
\alias{-.matrix.csr}
\alias{*.matrix.csr}
\alias{/.matrix.csr}
\alias{^.matrix.csr}
\alias{\%\%.matrix.csr}
\alias{\%/\%.matrix.csr}
\alias{>.matrix.csr}
\alias{>=.matrix.csr}
\alias{<.matrix.csr}
\alias{<=.matrix.csr}
\alias{==.matrix.csr}
\alias{!=.matrix.csr}
\alias{&.matrix.csr}
\alias{|.matrix.csr}
\alias{t.matrix.csr} 
\alias{t.matrix.csc} 
\alias{diag.matrix.csr} 
\alias{diag.assign.matrix.csr} 
\alias{ncol.matrix.csr} 
\alias{nrow.matrix.csr} 
\alias{dim.matrix.csr} 
\alias{dim.matrix.csc} 
\alias{dim.matrix.ssr} 
\alias{dim.matrix.ssc} 
\alias{rbind.matrix.csr} 
\alias{cbind.matrix.csr} 
\alias{[.matrix.csr}
\alias{[<-.matrix.csr}
\title{Basic Linear Algebra for Sparse Matrices}
\description{Basic linear algebra operations for sparse matrices
of class \code{matrix.csr}.  
}
\usage{
t(x); diag(x, nrow); diag(x) <- value; ncol(x); nrow(x); dim(x);rbind(...);
cbind(...);x[i,j]; x \%*\% y; x \%\% y; x \%/\% y; x + y; x - y; x * y;
x / y; x ^ y; x > y; x >= y; x < y; x <= y; x == y; x != y; x & y; x | y
}
\arguments{
\item{x}{matrix of class \code{matrix.csr}.}
\item{y}{matrix of class \code{matrix.csr} or a dense vector.}
\item{value}{replacement values.}
\item{i,j}{vectors of elements to extract or replace.}
\item{nrow}{optional number of rows for the result.}
}
\details{Linear algebra operations for matrices of class 
\code{matrix.csr} are designed to behave exactly as for 
regular matrices.  In particular, matrix multiplication, addition, 
subtraction and various logical operations work as with the conventional
dense form of matrix storage, as does indexing, rbind, cbind, and diagonal
assignment and extraction.
 
}
\references{ 
Koenker, R and Ng, P. (2002).  SparseM:  A Sparse Matrix Package for \R,\cr
\url{http://www.econ.uiuc.edu/~roger/research}
}
\seealso{
\code{slm} for sparse linear model fitting.  
\code{SparseM.ontology} for coercion and other class relations involving the
sparse matrix classes.
}
\examples{
n1 <- 10
n2 <- 10
p <- 6
y <- rnorm(n1)
a <- rnorm(n1*p)
a[abs(a)<0.5] <- 0
A <- matrix(a,n1,p)
A.csr <- as.matrix.csr(A)
b <- rnorm(n2*p)
b[abs(b)<1.0] <- 0
B <- matrix(b,n2,p)
B.csr <- as.matrix.csr(B)

# matrix transposition and multiplication
A.csr\%*\%t(B.csr)

}

\keyword{algebra}
