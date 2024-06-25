\name{SparseM.solve}
\title{Linear Equation Solving via Cholesky Decomposition for Sparse Matrices}
\alias{SparseM.solve}
\alias{chol,ANY-method}
\alias{chol,matrix.csr-method}
\alias{chol,matrix.csc-method}
\alias{chol,matrix-method}
\alias{chol}
% \alias{backsolve}
\alias{backsolve-methods}
\alias{backsolve,ANY-method}
\alias{backsolve,matrix.csr.chol-method}
\alias{forwardsolve}
\alias{forwardsolve,matrix.csr.chol-method}
\alias{solve}
\alias{solve,ANY-method}
\alias{solve,matrix.csr-method}
\description{ \describe{
    \item{\code{chol()}}{performs a Cholesky decomposition of a symmetric
      positive definite sparse matrix \code{x} of class \code{matrix.csr}.}
    \item{\code{backsolve()}}{performs a triangular back-fitting to compute
      the solutions of a system of linear equations in one step.}
    \item{\code{backsolve()} and \code{forwardsolve()}}{can also split the functionality of
      \code{backsolve} into two steps.}
    \item{\code{solve()}}{combines \code{chol()} and \code{backsolve()} to
      compute the inverse of a matrix if the right-hand-side is missing.}
  }
}
\usage{
chol(x, \dots)
\S4method{chol}{matrix.csr}(x, pivot = FALSE,
    nsubmax, nnzlmax, tmpmax,
    eps = .Machine$double.eps, tiny = 1e-30, Large = 1e128, warnOnly = FALSE,
    cacheKb = 1024L, level = 8L, \dots)

\S4method{backsolve}{matrix.csr.chol}(r, x, k, upper.tri, transpose,
          twice = TRUE, drop = TRUE, \dots)
\S4method{forwardsolve}{matrix.csr.chol}(l, x, k, upper.tri, transpose)
\S4method{solve}{matrix.csr}(a, b, \dots)
}
\arguments{
\item{a}{symmetric positive definite matrix of class \code{"matrix.csr"}.}
\item{r, l}{object of class \code{"matrix.csr.chol"} as returned by the
  \code{chol()} method.}
\item{x}{\describe{
    \item{For \code{chol()}:}{One of the sparse matrix classes,
      \code{"matrix.csr"} or \code{"matrix.csc"};}
    \item{For \code{{back,forward,}solve()}:}{vector or regular matrix of
      right-hand-side(s) of a system of linear equations.}
}}
\item{b}{vector or matrix right-hand-side(s) to solve for.}
\item{k}{inherited from the generic; not used here.}
\item{pivot}{inherited from the generic; not used here.}
\item{nsubmax, nnzlmax, tmpmax}{positive integer numbers with smart
  defaults; do \emph{not} set unless you know what you are doing!}
\item{eps}{positive tolerance for checking symmetry; change with caution.}
\item{tiny}{positive tolerance for checking diagonal entries to be
  \dQuote{essentially zero} and hence to be replaced by \code{Large}, during
  Cholesky decomposition.  Chaning this value may help in close to
  singular cases, see \sQuote{Examples}.}
\item{Large}{large positive number, \dQuote{essentially infinite}, to
  replace tiny diagonal entries during Cholesky.}
\item{warnOnly}{\code{\link{logical}}; when set to true, a result is
  returned with a \code{\link{warning}} instead of an error (via
  \code{\link{stop}()}); notably in close to singular cases.}
\item{cacheKb}{a positive integer, specifying an approximate size of the machine's
  cache memory in kilo (1024) bytes (\sQuote{Kb}); used to be hard wired to 64.}
\item{level}{level of loop unrolling while performing numerical
  factorization; an integer in \code{c(1, 2, 4, 8)}; used to be hard wired to 8.}
\item{upper.tri, transpose}{inherited from the generic; not used here.}
\item{twice}{\code{\link{logical}} flag:  If true, \code{backsolve()}
  solves twice, see below.}
\item{drop}{\code{\link{logical}} flag:  If true, \code{backsolve()}
  returns \code{\link{drop}(.)}, i.e., a vector instead of a column-1 matrix.}
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
There is no explicit checking for positive definiteness of the matrix
so users are advised to ensure that this condition is satisfied.
Messages such as "insufficient space" may indicate that one is trying
to factor a singular matrix.
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
  Koenker, R and Ng, P. (2002)
  SparseM: A Sparse Matrix Package for \R. \url{http://www.econ.uiuc.edu/~roger/research/home.html}

  Ng, E. G. and B. W. Peyton (1993)
  Block sparse Cholesky algorithms on advanced uniprocessor computers.
  \emph{SIAM J. Scientific Computing} \bold{14}, 1034--1056.
}

\seealso{
  \code{\link{slm}()} for a sparse version of \pkg{stats} package's \code{\link[stats]{lm}()}.
}
\examples{
data(lsq)
class(lsq) # -> [1] "matrix.csc.hb"
model.matrix(lsq)->design.o
class(design.o) # -> "matrix.csr"
dim(design.o) # -> [1] 1850  712
y <- model.response(lsq) # extract the rhs
length(y) # [1] 1850

X <- as.matrix(design.o)
c(object.size(X) / object.size(design.o)) ## X is 92.7 times larger

t(design.o) \%*\% design.o -> XpX
t(design.o) \%*\% y -> Xpy
chol(XpX) -> chol.o

determinant(chol.o)

b1 <- backsolve(chol.o,Xpy) # least squares solutions in two steps
b2 <- solve(XpX,Xpy)        # least squares estimates in one step
b3 <- backsolve(chol.o, forwardsolve(chol.o, Xpy),
                twice = FALSE) # in three steps
## checking that these three are indeed equal :
stopifnot(all.equal(b1, b2), all.equal(b2, b3))
}
\keyword{algebra}
