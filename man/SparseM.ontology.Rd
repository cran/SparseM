\name{SparseM.ontology}
\alias{SparseM.ontology}
\alias{matrix.csr}
\alias{matrix.csc}
\alias{initialize,ANY-method}
\alias{initialize,matrix.csr-method}
\alias{initialize,matrix.coo-method}
\alias{coerce,vector,matrix.diag.csr-method}
\alias{coerce,matrix,matrix.csr-method}
\alias{coerce,matrix.csr,matrix.diag.csr-method}
\alias{coerce,vector,matrix.csr-method}
\alias{coerce,numeric,matrix.diag.csr-method}
\alias{as.matrix,ANY-method}
\alias{as.matrix,matrix.csr-method}
\alias{as.matrix,csr-method}
\alias{is.matrix,csr-methods} 
\alias{as.matrix,csc-methods}
\alias{is.matrix,csc-methods} 
\alias{as.matrix,ssr-methods}
\alias{is.matrix,ssr-methods} 
\alias{as.matrix,ssc-methods}
\alias{is.matrix,ssc-methods} 
\alias{as.matrix,coo-methods}
\alias{is.matrix,coo-methods} 
\alias{as.matrix.csr,ANY-method}
\alias{as.matrix.csr,matrix.csc-method}
\alias{as.matrix.csr,matrix.ssr-method}
\alias{as.matrix.csr,matrix.ssc-method}
\alias{as.matrix.csr,matrix.coo-method}
\alias{as.matrix.csr,matrix.csr.chol-method}
\alias{as.matrix.csc,ANY-method}
\alias{as.matrix.csc,matrix.csr-method}
\alias{as.matrix.csc,matrix.csc-method}
\alias{as.matrix.csc,matrix.ssr-method}
\alias{as.matrix.csc,matrix.ssc-method}
\alias{as.matrix.csc,matrix.coo-method}
\alias{as.matrix.ssr,ANY-method}
\alias{as.matrix.ssr,matrix.csc-method}
\alias{as.matrix.ssr,matrix.ssr-method}
\alias{as.matrix.ssr,matrix.ssc-method}
\alias{as.matrix.ssr,matrix.coo-method}
\alias{as.matrix.ssc,ANY-method}
\alias{as.matrix.ssc,matrix.csr-method}
\alias{as.matrix.ssc,matrix.csc-method}
\alias{as.matrix.ssc,matrix.ssr-method}
\alias{as.matrix.ssc,matrix.ssc-method}
\alias{as.matrix.ssc,matrix.coo-method}
\alias{as.matrix.coo,ANY-method}
\alias{as.matrix.coo,matrix.csr-method}
\alias{as.matrix.coo,matrix.csc-method}
\alias{as.matrix.coo,matrix.ssr-method}
\alias{as.matrix.coo,matrix.ssc-method}
\alias{as.matrix.coo,matrix.coo-method}
\alias{as.matrix,matrix.csc-method}      
\alias{as.matrix,matrix.ssc-method}     
\alias{as.matrix,matrix.ssr-method}      
\alias{as.matrix,matrix.coo-method}      
\alias{matrix.ssc}            
\alias{as.matrix.ssc,matrix.csc-method} 
\alias{as.matrix.ssc.matrix.csr-method}  
\alias{as.matrix.ssc,matrix.ssr-method} 
\alias{matrix.ssr}             
\alias{as.matrix.ssr,matrix.csc-method}  
\alias{as.matrix.ssr,matrix.csr-method} 
\alias{as.matrix.ssr,matrix.ssc-method}  
\alias{is.matrix.csr}
\alias{is.matrix.csc}
\alias{is.matrix.ssr}
\alias{is.matrix.ssc}
\alias{is.matrix.coo}
\alias{is.matrix.csr}
\alias{as.matrix.csr}
\alias{as.matrix.csc}
\alias{as.matrix.ssr}
\alias{as.matrix.ssc}
\alias{as.matrix.coo}
\title{Sparse Matrix Class}
\description{
This group of functions evaluates and coerces changes in class structure.
}
\usage{
\method{as}{matrix.csr}(x, nrow = 1, ncol = 1, eps = .Machine$double.eps, ...)
\method{as}{matrix.csc}(x, nrow = 1, ncol = 1, eps = .Machine$double.eps, ...)
\method{as}{matrix.ssr}(x, nrow = 1, ncol = 1, eps = .Machine$double.eps, ...)
\method{as}{matrix.ssc}(x, nrow = 1, ncol = 1, eps = .Machine$double.eps, ...)
\method{is}{matrix.csr}(x, ...)
\method{is}{matrix.csc}(x, ...)
\method{is}{matrix.ssr}(x, ...)
\method{is}{matrix.ssc}(x, ...)
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
as(4,"matrix.diag.csr") #identity matrix of dimension 4
}
\keyword{algebra}
