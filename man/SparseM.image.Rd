\name{SparseM.image}
\alias{image,matrix.csr-method}
\alias{image.matrix.csc-method}
\alias{image,ANY-method}
\alias{SparseM.image}
\alias{image}
\title{Image Plot for Sparse Matrices}
\description{Display the pattern of non-zero entries of 
a matrix of class \code{matrix.csr} or \code{matrix.csc}}
\usage{
image(x, \dots)
}
\arguments{
\item{x}{a matrix of class \code{matrix.csr} or \code{matrix.csc}.}
\item{\dots}{additional arguments.}
}
\details{The pattern of the non-zero entries of a sparse matrix 
is displayed.  By default nonzero entries of the matrix appear
as gray blocks and zero entries as white background.}
\references{
Koenker, R and Ng, P. (2002).  SparseM:  A Sparse Matrix Package for \R,\cr
\url{http://www.econ.uiuc.edu/~roger/research}}
\seealso{
\code{SparseM.ops}, 
\code{SparseM.solve}, 
\code{SparseM.ontology} 
}
\examples{
a <- rnorm(20*5)
A <- matrix(a,20,5)
A[row(A)>col(A)+4|row(A)<col(A)+3] <- 0
b <- rnorm(20*5)
B <- matrix(b,20,5)
B[row(A)>col(A)+2|row(A)<col(A)+2] <- 0
image(as.matrix.csr(A)\%*\%as.matrix.csr(t(B)))
}
\keyword{hplot}
\keyword{algebra}
