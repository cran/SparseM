\name{slm.methods}
\alias{slm.methods}
\alias{summary.slm}
\alias{summary.mslm}
\alias{print.summary.slm}
\alias{print.slm}
\alias{fitted.slm}
\alias{residuals.slm}
\alias{coef.slm}
\title{Methods for slm objects}
\description{
  Summarize, print, and extract objects from  \code{\link{slm}} objects.
}
\usage{
summary.slm(object, correlation, \dots)
summary.mslm(object, \dots)
print.slm(x, digits, \dots)
print.summary.slm(x, digits, symbolic.cor, signif.stars, \dots)
fitted.slm(object, \dots)
residuals.slm(object, \dots)
coef.slm(object, \dots)
}
\arguments{
\item{object,x}{ object of class \code{slm}.  }
\item{digits}{minimum number of significant digits to be used for most numbers.}
\item{symbolic.cor}{logical; if \code{TRUE}, the correlation of coefficients
      will be printed. The default is \code{FALSE}}
\item{signif.stars}{logical; if \code{TRUE}, P-values are additionally encoded
          visually as ``significance stars'' in order to help scanning
          of long coefficient tables.  It defaults to the
          `show.signif.stars' slot of `options'. }
\item{correlation}{logical; if \code{TRUE}, the correlation matrix of the
          estimated parameters is returned and printed.}
\item{\dots}{additional arguments passed to methods.}
}
\value{
  \code{print.slm} and \code{print.summary.slm} return invisibly.
  \code{fitted.slm}, \code{residuals.slm}, and \code{coef.slm} 
  return the corresponding components of the \code{slm} object.
}
\references{
Koenker, R and Ng, P. (2002).  SparseM:  A Sparse Matrix Package for \R,\cr
\url{http://www.econ.uiuc.edu/~roger/research}
}

\author{ Roger Koenker }

\seealso{  
\code{slm} 
}
\examples{
data(lsq)
X <- model.matrix(lsq) #extract the design matrix
y <- model.response(lsq) # extract the rhs
X1 <- as.matrix(X)
slm.time <- unix.time(slm(y~X1-1) -> slm.o) # pretty fast
cat("slm time =",slm.time,"\n")
cat("slm Results: Reported Coefficients Truncated to 5  ","\n")
sum.slm <- summary(slm.o)
sum.slm$coef <- sum.slm$coef[1:5,]
sum.slm
fitted(slm.o)[1:10]
residuals(slm.o)[1:10]
coef(slm.o)[1:10]
}

\keyword{ regression }
