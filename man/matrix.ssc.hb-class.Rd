\name{matrix.ssc.hb-class}
\docType{class}
\alias{matrix.ssc.hb-class}
\title{Class "matrix.ssc.hb"}
\description{A new class consists of the coefficient matrix and
the right-hand-side of a linear system of equations,
initial guess of the solution and the exact solutions if they exist stored
in external files using the Harwell-Boeing format.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("matrix.ssc.hb", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{ra}:}{Object of class \code{numeric}, ra component of the csc or ssc format of the coefficient matrix, X.}
    \item{\code{ja}:}{Object of class \code{integer}, ja component of the csc or s
sc format of the coefficient matrix, X.}
    \item{\code{ia}:}{Object of class \code{integer}, ia component of the csc or ssc format of the coefficient matrix, X.}
    \item{\code{rhs.ra}:}{Object of class \code{numeric}, ra component of the right-hand-side, y, if stored in csc or
        ssc format; right-hand-side stored in dense vector or matrix otherwise.}
    \item{\code{guess}:}{Object of class \code{numeric or NULL} vector of the initial guess of the solutions if they exist;
        a null vector otherwise.}
    \item{\code{xexact}:}{Object of class \code{numeric or NULL} vector of the exact solutions, b, if they exist; a null vector otherwise.}
    \item{\code{dimension}:}{Object of class \code{integer}, dimenson of the coefficient matrix, X.}
    \item{\code{rhs.dim}:}{Object of class \code{integer}, dimenson of the right-hand-side, y.}
    \item{\code{rhs.mode}:}{Object of class \code{character or NULL} storage mode of the right-hand-side; can be full storage or
        same format as the coefficient matrix.}
  }
}
\section{Extends}{
Class \code{"matrix.csc.hb"}, directly.
}
\section{Methods}{
  \describe{
    \item{model.matrix}{\code{signature(object = "matrix.ssc.hb")}: ... }
  }
}
\seealso{\code{\link{model.matrix}}, \code{\link{model.response}}, 
\code{\link{read.matrix.hb}}, \code{\link{matrix.csc.hb-class}}}
\keyword{classes}

