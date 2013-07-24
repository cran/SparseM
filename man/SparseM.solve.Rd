\name{SparseM.solve}
\alias{SparseM.solve}
\alias{chol,ANY-method}
\alias{chol,matrix.csr-method}
\alias{chol,matrix.csc-method}
\alias{chol,matrix-method}
\alias{chol}
\alias{backsolve-methods}
\alias{backsolve,ANY-method}
\alias{backsolve}
\alias{forwardsolve}
\alias{solve}
\alias{backsolve,matrix.csr.chol-method}
\alias{forwardsolve,matrix.csr.chol-method}
\alias{solve,ANY-method}
\alias{solve,matrix.csr-method}
\title{Linear Equation Solving for Sparse Matrices}
\description{
\code{chol} performs a Cholesky
decomposition of a symmetric positive definite sparse matrix \code{x}
of class \code{matrix.csr}. \cr
\code{backsolve} performs a triangular back-fitting to compute
the solutions of a system of linear equations in one step. \cr
\code{backsolve} and \code{forwardsolve} can also split the functionality of
\code{backsolve} into two steps. \cr
\code{solve} combines \code{chol} and \code{backsolve} and will
compute the inverse of a matrix if the right-hand-side is missing.
}
\usage{
chol(x, \dots)%%MM chol(x, pivot = FALSE, \dots)
\S4method{backsolve}{matrix.csr.chol}(r, x, k = NULL, upper.tri = NULL,
          transpose = NULL, twice = TRUE, \dots)
forwardsolve(l, x, k = ncol(l), upper.tri = FALSE, transpose = FALSE)
solve(a, b, \dots)
}
\arguments{
\item{a}{symmetric positive definite matrix of class \code{matrix.csr}.}
\item{r}{object of class \code{matrix.csr.chol} returned by the function \code{chol}.}
\item{l}{object of class \code{matrix.csr.chol} returned by the function \code{chol}.}
\item{x,b}{vector(regular matrix) of right-hand-side(s) of a system of linear equations.}
\item{k}{inherited from the generic; not used here.}
%%MM \item{pivot}{inherited from the generic; not used here.}
\item{upper.tri}{inherited from the generic; not used here.}
\item{transpose}{inherited from the generic; not used here.}
\item{twice}{Logical flag:  If true backsolve solves twice, see below.}
\item{\dots}{further arguments passed to or from other methods.}
}
\details{
\code{chol} performs a Cholesky decomposition of
a symmetric positive definite sparse matrix \code{a} of class
\code{matrix.csr} using the block sparse Cholesky algorithm of Ng and
Peyton (1993).  The structure of the resulting \code{matrix.csr.chol}
object is relatively complicated.  If necessary it can be coerced back
to a \code{matrix.csr} object as usual with \code{as.matrix.csr}.
\code{backsolve} does triangular back-fitting to compute
the solutions of a system of linear equations.  For systems of linear equations
that only vary on the right-hand-side, the result from \code{chol}
can be reused.  Contrary to the behavior of \code{backsolve} in base R,
the default behavior of  \code{backsolve(C,b)} when C is a \code{matrix.csr.chol} object
is to produce a solution to the system \eqn{Ax = b} where \code{C <- chol(A)}, see
the example section.  When the flag \code{twice} is \code{FALSE} then backsolve
solves the system \eqn{Cx = b}, up to a permutation  -- see the comments below.
The command \code{solve} combines \code{chol} and \code{backsolve}, and will
compute the inverse of a matrix if the right-hand-side is missing.
The determinant of the Cholesky factor is returned providing a
means to efficiently compute the determinant of sparse positive
definite symmetric matrices.

There are several integer storage parameters that are set by default in the call
to the Cholesky factorization, these can be overridden in any of the above
functions and will be passed by the usual "dots" mechanism.  The necessity
to do this is usually apparent from error messages like:  Error  
in local(X...) increase tmpmax.   For example, one can use,
\code{solve(A,b, tmpmax = 100*nrow(A))}.  The current default for tmpmax
is \code{50*nrow(A)}. Some experimentation may be needed to
select appropriate values, since they are highly problem dependent.  See
the code of chol() for further details on the current defaults.
}

\note{
Because the sparse Cholesky algorithm re-orders the positive
definite sparse matrix \code{A}, the value of
\code{x <- backsolve(C, b)} does \emph{not} equal the solution to the
triangular system \eqn{Cx = b}, but is instead the solution to the
system \eqn{CPx = Pb} for some permutation matrix \eqn{P}
(and analogously for \code{x <- forwardsolve(C, b)}).  However, a
little algebra easily shows that
\code{backsolve(C, forwardsolve(C, b), twice = FALSE)} \emph{is} the solution
to the equation \eqn{Ax=b}.  Finally, if \code{C <- chol(A)}  for some
sparse covariance matrix \code{A}, and z is a conformable standard normal vector,
then the product  \code{y <- as.matrix.csr(C) \%*\% z} is normal with covariance
matrix \code{A} irrespective of the permutation of the Cholesky factor.
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

t(design.o) \%*\% design.o -> XpX
t(design.o) \%*\% y -> Xpy
chol(XpX) -> chol.o

b1 <- backsolve(chol.o,Xpy) # least squares solutions in two steps
b2 <- solve(XpX,Xpy)        # least squares estimates in one step
b3 <- backsolve(chol.o, forwardsolve(chol.o, Xpy),
                twice = FALSE) # in three steps
## checking that these three are indeed equal :
stopifnot(all.equal(b1, b2), all.equal(b2, b3))
}
\keyword{algebra}
