\name{SparseM.solve}
\alias{SparseM.solve}
\alias{chol.matrix.csr}              
\alias{chol.matrix.csc}              
\alias{chol}              
\alias{chol.default}              
\alias{backsolve}
\alias{backsolve.default}
\alias{backsolve.matrix.csr.chol}
\alias{backsolve.matrix.csr.chol}
\alias{solve.matrix.csr}
\title{Linear Equation Solving for Sparse Matrices}
\description{
\code{chol.matrix.csr} performs a Cholesky
decomposition of a symmetric positive definite sparse matrix \code{x} 
of class \code{matrix.csr}. \cr
\code{backsolve.matrix.csr} performs a triangular back-fitting to compute
the solutions of a system of linear equations.  \cr
\code{solve.matrix.csr} 
combines \code{chol.matrix.csr} and \code{backsolve.matrix.csr} and will
compute the inverse of a matrix if the right-hand-side is missing.
}
\usage{
\method{chol}{matrix.csr}(x, cachsz=64, nsubmax, nnzlmax, tmpmax, \dots)
\method{backsolve}{matrix.csr.chol}(r, x, cachsz=64, \dots)
\method{solve}{matrix.csr}(a, b, \dots)
}
\arguments{
\item{X,a}{symmetric positive definite matrix of class \code{matrix.csr}.}
\item{r}{object of class \code{matrix.csr.chol} returned by the function
\code{chol.matrix.csr}.}
\item{x,b}{vector(regular matrix) of right-hand-side(s) of a system of linear equations.}
\item{cachsz}{size of the cache memory; default to 64.}
\item{nsubmax,nnzlmax,tmpmax}{dimensions of work arrays, in normal
operation these are determined inside the algorithm.}
\item{\dots}{additional arguments to be passed to \code{chol.matrix.csr}.}
}
\details{
\code{chol.matrix.csr} performs a Cholesky decomposition of
a symmetric positive definite sparse matrix \code{x} of class 
\code{matrix.csr} using the block sparse Cholesky algorithm of Ng and 
Peyton (1993).
\code{backsolve.matrix.csr} does triangular back-fitting to compute
the solutions of a system of linear equations.  For systems of linear equations
that only vary on the right-hand-side, the result from \code{chol}
can be reused.
\code{solve.matrix.csr} 
combines \code{chol.matrix.csr} and \code{backsolve.matrix.csr}, and will
compute the inverse of a matrix if the right-hand-side is missing.
}
\references{ 
Koenker, R and Ng, P. (2002).  SparseM:  A Sparse Matrix Package for \R,\cr
\url{http://www.econ.uiuc.edu/~roger/research}

Ng, E. G. and B. W. Peyton (1993), "Block sparse Cholesky algorithms
on advanced uniprocessor computers", \emph{SIAM J. Sci. Comput.}, \bold{14},
pp. 1034-1056.

}
\seealso{
\code{slm} for sparse version of \code{lm} \cr
}
\examples{
# lsq.rra is real rectangular stored in csc (compressed sparse column) format
read.matrix.hb(system.file("HBdata","lsq.rra",package = "SparseM"))-> hb.o
class(hb.o) # -> [1] "matrix.csc.hb"
model.matrix.hb(hb.o)->design.o
class(design.o) # -> "matrix.csr"
dim(design.o) # -> [1] 1850  712
y <- model.response.hb(hb.o) # extract the rhs
length(y) # [1] 1850
t(design.o)\%*\%design.o -> XpX
t(design.o)\%*\%y -> Xpy
chol(XpX)->chol.o
backsolve(chol.o,Xpy)-> b1 # least squares solutions in two steps
solve(XpX,Xpy) -> b2 # least squares estimates in one step
}

\keyword{algebra}
