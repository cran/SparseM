\name{SparseM.solve}
\alias{SparseM.solve}
\alias{chol,ANY-method}
\alias{chol,matrix.csr-method}
\alias{chol,matrix.csc-method}              
\alias{chol,matrix-method}              
\alias{chol}              
\alias{chol.default}              
\alias{backsolve-methods}
\alias{backsolve,ANY-method}
\alias{backsolve}
\alias{solve}
\alias{backsolve,matrix.csr.chol-method}
\alias{solve,ANY-method}
\alias{solve,matrix.csr-method}
\title{Linear Equation Solving for Sparse Matrices}
\description{
\code{chol} performs a Cholesky
decomposition of a symmetric positive definite sparse matrix \code{x} 
of class \code{matrix.csr}. \cr
\code{backsolve} performs a triangular back-fitting to compute
the solutions of a system of linear equations.  \cr
\code{solve} 
combines \code{chol} and \code{backsolve} and will
compute the inverse of a matrix if the right-hand-side is missing.
}
\usage{
chol(x, \dots)
backsolve(r, x, k, upper.tri, transpose)
solve(a, b, \dots)
}
\arguments{
\item{a}{symmetric positive definite matrix of class \code{matrix.csr}.}
\item{r}{object of class \code{matrix.csr.chol} returned by the function
\code{chol}.}
\item{x,b}{vector(regular matrix) of right-hand-side(s) of a system of linear equations.}
\item{k}{inherited from the generic; not used here.}
#\item{pivot}{inherited from the generic; not used here.}
#\item{nsubmax,nnzlmax,tmpmax}{dimensions of work arrays, in normal operation these are determined inside the algorithm.}
\item{upper.tri}{inherited from the generic; not used here.}
\item{transpose}{inherited from the generic; not used here.}
\item{\dots}{further arguments passed to or from other methods.}
}
\details{
\code{chol} performs a Cholesky decomposition of
a symmetric positive definite sparse matrix \code{x} of class 
\code{matrix.csr} using the block sparse Cholesky algorithm of Ng and 
Peyton (1993).
\code{backsolve} does triangular back-fitting to compute
the solutions of a system of linear equations.  For systems of linear equations
that only vary on the right-hand-side, the result from \code{chol}
can be reused.
\code{solve} 
combines \code{chol} and \code{backsolve}, and will
compute the inverse of a matrix if the right-hand-side is missing.
The determinant of the Cholesky factor is returned providing a
means to efficiently compute the determinant of sparse positive
definite symmetric matrices.
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
data(lsq)
class(lsq) # -> [1] "matrix.csc.hb"
model.matrix(lsq)->design.o
class(design.o) # -> "matrix.csr"
dim(design.o) # -> [1] 1850  712
y <- model.response(lsq) # extract the rhs
length(y) # [1] 1850
t(design.o)\%*\%design.o -> XpX
t(design.o)\%*\%y -> Xpy
chol(XpX)->chol.o
backsolve(chol.o,Xpy)-> b1 # least squares solutions in two steps
solve(XpX,Xpy) -> b2 # least squares estimates in one step
}

\keyword{algebra}
