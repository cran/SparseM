\name{slm-class}
\docType{class}
\alias{slm-class}
\title{Class "slm"}
\description{A sparse extension of \code{lm}}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("slm", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{coefficients}:}{Object of class \code{"numeric"} estimated coefficients}
    \item{\code{chol}:}{Object of class \code{"matrix.csr.chol"} Cholesky object from fitting}
    \item{\code{residuals}:}{Object of class \code{"numeric"} residuals}
    \item{\code{fitted}:}{Object of class \code{"numeric"} fitted values}
  }
}
\section{Extends}{
Class \code{"lm"}, directly.
Class \code{"oldClass"}, by class "lm".
}
\section{Methods}{
  \describe{
    \item{coef}{\code{signature(object = "slm")}: ... }
    \item{fitted}{\code{signature(object = "slm")}: ... }
    \item{residuals}{\code{signature(object = "slm")}: ... }
    \item{summary}{\code{signature(object = "slm")}: ... }
  }
}
\seealso{ \code{\link{slm}}}
\keyword{classes}

