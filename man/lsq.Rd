\name{lsq}
\alias{lsq}
\title{Least Squares Problems in Surveying}
\description{One of the four matrices from the least-squares solution of 
problems in surveying that were used by Michael Saunders and Chris Paige 
in the testing of LSQR}
\usage{data(lsq)}
\format{
A list of class \code{matrix.csc.hb} or \code{matrix.ssc.hb} depending
on how the coefficient matrix is stored with the following components:
\itemize{
  \item{ra }{ra component of the csc or ssc format of the coefficient matrix, X.
}   \item{ja }{ja component of the csc or ssc format of the coefficient matrix, X.
}   \item{ia }{ia component of the csc or ssc format of the coefficient matrix, X.
}   \item{rhs.ra }{ra component of the right-hand-side, y, if stored in csc or
        ssc format; right-hand-side stored in dense vector or matrix otherwise.}
  \item{rhs.ja }{ja component of the right-hand-side, y, if stored in csc or 
        ssc format; a null vector otherwise.}
  \item{rhs.ia }{ia component of the right-hand-side, y, if stored in csc or
        ssc format; a null vector otherwise.}
  \item{xexact}{vector of the exact solutions, b, if they exist; a null vector o
therwise.}  \item{guess}{vector of the initial guess of the solutions if they exist;
        a null vector otherwise.}
  \item{dim}{dimenson of the coefficient matrix, X.}
  \item{rhs.dim}{dimenson of the right-hand-side, y.}
  \item{rhs.mode}{storage mode of the right-hand-side; can be full storage or
        same format as the coefficient matrix.}
}
}
\seealso{\code{read.matrix.hb}}
\examples{
data(lsq)
class(lsq) # -> [1] "matrix.csc.hb"
model.matrix(lsq)->X
class(X) # -> "matrix.csr"
dim(X) # -> [1] 1850  712
y <- model.response(lsq) # extract the rhs
length(y) # [1] 1850 
}
\references{
Koenker, R and Ng, P. (2002).  SparseM:  A Sparse Matrix Package for \R,\cr
\url{http://www.econ.uiuc.edu/~roger/research}

Matrix Market, \url{http://math.nist.gov/MatrixMarket/data/Harwell-Boeing/lsq/lsq.html}
}
\keyword{datasets}
