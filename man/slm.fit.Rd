\name{slm.fit}
\alias{slm.fit}
\alias{slm.fit}
\alias{slm.wfit}
\alias{slm.fit.csr}
\title{Internal slm fitting functions}
\description{
 Fitting functions for sparse linear model fitting.
}
\usage{
slm.fit(x,y,method, \dots)
slm.wfit(x,y,weights,\dots)
slm.fit.csr(x, y, \dots)
}
\arguments{
\item{x}{ design matrix.  }
  \item{y}{ vector of response observations.  }
  \item{method}{ only \code{csr} is supported currently }
  \item{weights}{an optional vector of weights to be used in the fitting
          process. If specified, weighted least squares is used with
          weights `weights' (that is, minimizing 
	\deqn{\sum w_i*e_i^2}{sum w e^2}}
  \item{\dots}{additional arguments.}
}
\details{
   \code{slm.fit} and \code{slm.wfit} call \code{slm.fit.csr} 
   to do Cholesky decomposition
   and then backsolve to obtain the least squares estimated coefficients.
   These functions can be called directly if the user is willing to
   specify the design matrix in \code{matrix.csr} form.  This is often
   advantageous in large problems to reduce memory requirements.
}
\item{weights}{
    vector of observation weights to be used by `slm.wfit'; 
    if supplied, the algorithm fits
    to minimize the sum of the weights multiplied by the
    squared residuals. The length of weights must be the same as
    the number of observations.  The weights must be nonnegative
    and it is strongly recommended that they be strictly positive,
    since zero weights are ambiguous.
  }
\item{\dots}{
   additional arugments.
  }
}
}
\value{
A list of class \code{slm} consisting of:
  \item{coef}{estimated coefficients}
  \item{chol}{cholesky object from fitting}
  \item{residuals}{residuals}
  \item{fitted}{fitted values}
  \item{terms}{terms}
  \item{call}{call}
  ...
}
\references{
Koenker, R and Ng, P. (2002).  SparseM:  A Sparse Matrix Package for \R,\cr
\url{http://www.econ.uiuc.edu/~roger/research}}
}
\author{ Roger Koenker }

\seealso{  \code{\link{slm}} }
\example{
# lsq.rra is real rectangular matrix stored in compressed sparse column format
read.matrix.hb(system.file("HBdata","lsq.rra",package = "SparseM"))-> hb.o
X <- model.matrix.hb(hb.o) #extract the design matrix
y <- model.response.hb(hb.o) # extract the rhs
class(X) # -> "matrix.csr"
class(y) # -> NULL
slm.fit(X,y)->slm.fit.o # this is much more efficient in memory usage than slm()
slm(y~as.matrix(X)-1) -> slm.o # this requires X to be transformed into dense mode
cat("Difference between `slm.fit' and `slm' estimated coefficients =",sum(abs(slm.fit.o$coef-slm.o$coef)),"\n")
}

\keyword{ regression }
