\name{SparseM.ontology}
\alias{SparseM.ontology}
\alias{matrix.csr}
\alias{as.matrix.matrix.csr} 
\alias{as.matrix.csr}
\alias{is.matrix.csr} 
\alias{as.matrix.csc}
\alias{is.matrix.csc} 
\alias{as.matrix.ssr}
\alias{is.matrix.ssr} 
\alias{as.matrix.ssc}
\alias{is.matrix.ssc} 
\alias{as.matrix.coo}
\alias{is.matrix.coo} 
\alias{matrix.csc}
\alias{as.matrix.csc.matrix.csr} 
\alias{as.matrix.csc.matrix.ssc}  
\alias{as.matrix.csc.matrix.ssr} 
\alias{as.matrix.csr.matrix.csc} 
\alias{as.matrix.csr.matrix.ssc}  
\alias{as.matrix.csr.matrix.ssr} 
\alias{as.matrix.csr.matrix.coo} 
\alias{as.matrix.coo.matrix.csr} 
\alias{as.matrix.matrix.csc}      
\alias{as.matrix.matrix.ssc}     
\alias{as.matrix.matrix.ssr}      
\alias{as.matrix.matrix.coo}      
\alias{matrix.ssc}            
\alias{as.matrix.ssc.matrix.csc} 
\alias{as.matrix.ssc.matrix.csr}  
\alias{as.matrix.ssc.matrix.ssr} 
\alias{matrix.ssr}             
\alias{as.matrix.ssr.matrix.csc}  
\alias{as.matrix.ssr.matrix.csr} 
\alias{as.matrix.ssr.matrix.ssc}  
\title{Sparse Matrix Class}
\description{
This group of functions evaluates and coerces changes in class structure.
}
\usage{
as.matrix.csr(x, nrow = 1, ncol = 1, eps = .Machine$double.eps)
as.matrix.csc(x, nrow = 1, ncol = 1, eps = .Machine$double.eps)
as.matrix.ssr(x, nrow = 1, ncol = 1, eps = .Machine$double.eps)
as.matrix.ssc(x, nrow = 1, ncol = 1, eps = .Machine$double.eps)
is.matrix.csr(x, ...)
is.matrix.csc(x, ...)
is.matrix.ssr(x, ...)
is.matrix.ssc(x, ...)
}
\arguments{
\item{x}{is a matrix, or vector object,  of either dense or sparse form}
\item{nrow}{number of rows of matrix }
\item{ncol}{number of columns of matrix }
\item{eps}{A tolerance parameter:  elements of x such that abs(x) < eps set to zero.
This argument is only relevant when coercing matrices from dense to sparse form. Defaults to
\code{eps = .Machine$double.eps} }
\item{...}{other arguments}

}
\details{
The function \code{matrix.csc} acts like \code{matrix} to coerce a vector object to
a sparse matrix object of class \code{matrix.csr}.  
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

The compressed sparse column class  \code{matrix.csc} is defined in
an analogous way, as are  the \code{matrix.ssr}, symmetric sparse row, and
\code{matrix.ssc}, symmetric sparse column classes.
}
\note{
\code{as.matrix.ssr} and \code{as.matrix.ssc} should ONLY be used with 
symmetric matrices.
}
\references{
Koenker, R and Ng, P. (2002).  SparseM:  A Sparse Matrix Package for \R,\cr
\url{http://www.econ.uiuc.edu/~roger/research}
}
\seealso{
\code{SparseM.hb} for handling Harwell-Boeing sparse matrices.
}
\examples{
n1 <- 10
p <- 5
a <- rnorm(n1*p)
a[abs(a)<0.5] <- 0
A <- matrix(a,n1,p)
B <- t(A)\%*\%A
A.csr <- as.matrix.csr(A)
A.csc <- as.matrix.csc(A)
B.ssr <- as.matrix.ssr(B)
B.ssc <- as.matrix.ssc(B)
is.matrix.csr(A.csr) # -> TRUE
is.matrix.csc(A.csc) # -> TRUE
is.matrix.ssr(B.ssr) # -> TRUE
is.matrix.ssc(B.ssc) # -> TRUE
as.matrix(A.csr)
as.matrix(A.csc)
as.matrix(B.ssr)
as.matrix(B.ssc)
as.matrix.csr(rep(0,9),3,3) #sparse matrix of all zeros
}
\keyword{algebra}
