\name{SparseM.ontology}
\title{Sparse Matrix Class Ontology}
\alias{SparseM.ontology}
\alias{initialize,ANY-method}
\alias{initialize,matrix.csr-method}
\alias{initialize,matrix.coo-method}
\alias{coerce,vector,matrix.diag.csr-method}
\alias{coerce,matrix,matrix.csr-method}
\alias{coerce,matrix.csr,matrix.diag.csr-method}
\alias{coerce,vector,matrix.csr-method}
\alias{coerce,numeric,matrix.diag.csr-method}
\alias{as.matrix}
\alias{as.matrix,ANY-method}
\alias{as.matrix,matrix.csr-method}
\alias{as.matrix,matrix.csc-method}
\alias{as.matrix,matrix.ssc-method}
\alias{as.matrix,matrix.ssr-method}
\alias{as.matrix,matrix.coo-method}
%
\alias{as.matrix.coo}
\alias{as.matrix.coo,ANY-method}
\alias{as.matrix.coo,matrix.csr-method}
% \alias{as.matrix.coo,matrix.csc-method}
% \alias{as.matrix.coo,matrix.ssr-method}
% \alias{as.matrix.coo,matrix.ssc-method}
% \alias{as.matrix.coo,matrix.coo-method}
\alias{as.matrix.csr}
\alias{as.matrix.csr,ANY-method}
\alias{as.matrix.csr,matrix.csc-method}
\alias{as.matrix.csr,matrix.ssr-method}
\alias{as.matrix.csr,matrix.ssc-method}
\alias{as.matrix.csr,matrix.coo-method}
\alias{as.matrix.csr,matrix.csr.chol-method}
\alias{as.matrix.csc}
\alias{as.matrix.csc,ANY-method}
\alias{as.matrix.csc,matrix.csr-method}
\alias{as.matrix.csc,matrix.csc-method}
\alias{as.matrix.csc,matrix.ssr-method}
\alias{as.matrix.csc,matrix.ssc-method}
\alias{as.matrix.csc,matrix.coo-method}
\alias{as.matrix.ssr}
\alias{as.matrix.ssr,ANY-method}
\alias{as.matrix.ssr,matrix.coo-method}
\alias{as.matrix.ssr,matrix.csc-method}
\alias{as.matrix.ssr,matrix.csr-method}
\alias{as.matrix.ssr,matrix.ssc-method}
\alias{as.matrix.ssr,matrix.ssr-method}
\alias{as.matrix.ssc}
\alias{as.matrix.ssc,ANY-method}
\alias{as.matrix.ssc,matrix.coo-method}
\alias{as.matrix.ssc,matrix.csc-method}
\alias{as.matrix.ssc,matrix.csr-method}
\alias{as.matrix.ssc,matrix.ssc-method}
\alias{as.matrix.ssc,matrix.ssr-method}
% regular functions i.e. *not* methods:
\alias{is.matrix.coo}
\alias{is.matrix.csc}
\alias{is.matrix.csr}
\alias{is.matrix.ssc}
\alias{is.matrix.ssr}
% "non-sense" (have {matrix.csr-class} in ./matrix.csr-class.Rd), but \link{}ed from pkg-space:
\alias{matrix.csr}
%
\description{
  This group of functions evaluates and coerces changes in class structure.
}
\usage{
as.matrix.csr(x, nrow, ncol, eps = .Machine$double.eps, \dots)
\S4method{as.matrix.csr}{matrix.csr.chol}(x, nrow, ncol, eps, upper.tri=TRUE, \dots)
\S4method{as.matrix.csc}{matrix.csr}(x, nrow = 1, ncol = 1, eps = .Machine$double.eps)
\S4method{as.matrix.ssr}{matrix.coo}(x, nrow = 1, ncol = 1, eps = .Machine$double.eps)
\S4method{as.matrix.ssc}{matrix.csc}(x, nrow = 1, ncol = 1, eps = .Machine$double.eps)
\S4method{as.matrix.coo}{matrix.csr}(x, nrow = 1, ncol = 1, eps = .Machine$double.eps)

is.matrix.csr(x)
is.matrix.csc(x)
is.matrix.ssr(x)
is.matrix.ssc(x)
is.matrix.coo(x)
}
\arguments{
\item{x}{is a matrix, or vector object,  of either dense or sparse form}
\item{nrow}{number of rows of matrix }
\item{ncol}{number of columns of matrix }
\item{eps}{A tolerance parameter:  elements of x such that abs(x) < eps set to zero.
  This argument is only relevant when coercing matrices from dense to sparse form. Defaults to
  \code{eps = .Machine$double.eps} }
\item{upper.tri}{\code{\link{logical}}, to choose upper or lower triangular matrix result.}
\item{\dots}{other arguments}
}
\details{
The function \code{matrix.csc} acts like \code{matrix} to coerce a vector object to
a sparse matrix object of class \code{matrix.csr}.
This aspect of the code is in the process of conversion from S3 to S4 classes.
For the most part the S3 syntax prevails.  An exception is the code to
coerce vectors to diagonal matrix form which uses \code{as(v,"matrix.diag.csr"}.
The generic functions \code{as.matrix.xxx} coerce a matrix \code{x} into
a matrix of storage class \code{matrix.xxx}. The argument matrix \code{x}
may be of conventional dense form, or of any of the four supported
classes:  \code{matrix.csr, matrix.csc, matrix.ssr, matrix.ssc}.
The generic functions \code{is.matrix.xxx} evaluate whether the
argument is of class \code{matrix.xxx}.  The function
\code{as.matrix} transforms a matrix of any sparse class into conventional
dense form.  The primary storage class for sparse matrices is the
compressed sparse row \code{matrix.csr} class.
An \emph{n} by \emph{m} matrix \emph{A} with real elements \eqn{a_{ij}}{a_{ij}},
stored in \code{matrix.csr} format consists of three arrays:

\itemize{
\item \code{ra}: a real array of \emph{nnz} elements containing the non-zero
elements of \emph{A}, stored in row order. Thus, if \emph{i<j}, all elements of row \emph{i}
precede elements from row \emph{j}. The order of elements within the rows is immaterial.

\item \code{ja}: an integer array of \emph{nnz} elements containing the column
indices of the elements stored in \code{ra}.

\item \code{ia}: an integer array of \emph{n+1} elements containing pointers to
the beginning of each row in the arrays \code{ra} and \code{ja}. Thus
\code{ia[i]} indicates the position in the arrays \code{ra} and
\code{ja} where the \emph{i}th row begins. The last, \emph{(n+1)}st, element of
\code{ia} indicates where the \emph{n+1} row would start, if it existed.
}

The compressed sparse column class  \code{matrix.csc} is defined in
an analogous way, as are  the \code{matrix.ssr}, symmetric sparse row, and
\code{matrix.ssc}, symmetric sparse column classes.
}
\note{
\code{as.matrix.ssr} and \code{as.matrix.ssc} should ONLY be used with
symmetric matrices.

\code{as.matrix.csr(x)}, when \code{x} is an object of class \code{matrix.csr.chol}
(that is, an object returned by a call to \code{chol(a)} when \code{a}
is an object of class \code{matrix.csr} or \code{matric.csc}),
by default returns an upper triangular matrix, which
is \emph{not} consistent with the result of \code{chol} in the \pkg{base}
package.  To get an lower triangular \code{matric.csr} matrix, use either
\code{as.matrix.csr(x, upper.tri = FALSE)} or
\code{t(as.matrix.csr(x))}.
}
\references{
  Koenker, R and Ng, P. (2002)
  \emph{SparseM: A Sparse Matrix Package for \R}.
  \url{http://www.econ.uiuc.edu/~roger/research/home.html}
}
\seealso{
\code{SparseM.hb} for handling Harwell-Boeing sparse matrices.
}
\examples{
t(m5 <- as.matrix.csr(c(-1:1,0,0)))
t(M4 <- as.matrix.csc(c(0:2,0), 4))
 (S3 <- as.matrix.ssr(diag(x = 0:2))) # *symmetric*
stopifnot(identical(dim(m5), c(5L, 1L)),
          identical(dim(M4), c(4L, 1L)),
          identical(dim(S3), c(3L, 3L)))

n1 <- 10
p <- 5
a <- round(rnorm(n1*p), 2)
a[abs(a) < 0.7] <- 0
A <- matrix(a,n1,p)
B <- t(A) \%*\% A      % 
A.csr <- as.matrix.csr(A)
A.csc <- as.matrix.csc(A)
B.ssr <- as.matrix.ssr(B)
B.ssc <- as.matrix.ssc(B)
stopifnot(exprs = {
  is.matrix.csr(A.csr) # -> TRUE
  is.matrix.csc(A.csc) # -> TRUE
  is.matrix.ssr(B.ssr) # -> TRUE
  is.matrix.ssc(B.ssc) # -> TRUE
})
as.matrix(A.csr)
as.matrix(A.csc)
as.matrix(B.ssr)
as.matrix(B.ssc)
as.matrix.csr(0, 2,3)    # sparse matrix of all zeros
## Diagonal (sparse) : 
as(4,   "matrix.diag.csr") # identity matrix of dimension 4
as(2:0, "matrix.diag.csr") # diagonal 3x3 matrix
}
\keyword{algebra}
