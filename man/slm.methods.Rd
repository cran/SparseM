\name{slm.methods}
\alias{slm.methods}
\alias{summary.slm}
\alias{print.summary.slm}
\alias{fitted.slm}
\alias{residuals.slm}
\alias{coef.slm}
\title{Methods for slm objects}
\description{
  Summarize, print, and extract objects from  \code{\link{slm}} objects.
}
\usage{
\method{summary}{slm}(object, correlation = FALSE, \dots)
\method{print}{summary.slm}(x, digits = max(3, getOption("digits") - 3), 
	symbolic.cor = p > 4, signif.stars = getOption("show.signif.stars"), 
	\dots)
\method{fitted}{slm}(object, \dots)
\method{residuals}{slm}(object, \dots)
\method{coef}{slm}(object, \dots)
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

\keyword{ regression }
