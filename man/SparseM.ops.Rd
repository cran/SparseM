\name{SparseM.ops}
\alias{Ops.matrix.csr}
\alias{Ops.matrix.diag.csr}
\alias{%*%-methods}
\alias{%*%,ANY,ANY-method}
\alias{%*%,matrix.csr,matrix.csr-method}
\alias{%*%,matrix.csr,matrix-method}
\alias{%*%,matrix.csr,numeric-method}
\alias{%*%,matrix,matrix.csr-method}
\alias{%*%,numeric,matrix.csr-method}
\alias{%x%-methods}
\alias{%x%,ANY,ANY-method}
\alias{%x%,matrix.csr,matrix.csr-method}
\alias{%x%,matrix.csr,matrix-method}
\alias{%x%,matrix.csr,numeric-method}
\alias{%x%,matrix,matrix.csr-method}
\alias{%x%,numeric,matrix.csr-method}
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
\alias{[.matrix.coo}
\alias{[<-.matrix.coo}
\title{Basic Linear Algebra for Sparse Matrices}
\description{Basic linear algebra operations for sparse matrices
of class \code{matrix.csr}.  
}
\usage{ }
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

# kronecker product
A.csr \%x\% matrix(1:4,2,2)

}

\keyword{algebra}
