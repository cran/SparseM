\name{SparseM.hb}
\alias{SparseM.hb}
\alias{read.matrix.hb}
\alias{model.matrix,ANY-method}
\alias{model.matrix,matrix.csc.hb-method}
\alias{model.matrix,matrix.ssc.hb-method}
\alias{model.matrix.matrix.ssc.hb}
\alias{model.response,ANY-method}
\alias{model.response,matrix.csc.hb-method}
\alias{model.response,matrix.ssc.hb-method}
\alias{model.guess,matrix.csc.hb-method}
\alias{model.guess,matrix.ssc.hb-method}
\alias{model.xexact,matrix.csc.hb-method}
\alias{model.xexact,matrix.ssc.hb-method}
\alias{model.response}
\alias{model.guess}
\alias{model.xexact}
\alias{model.matrix}
\title{Harwell-Boeing Format Sparse Matrices}
\description{
Read, and extract components of data in Harwell-Boeing sparse matrix format.
}
\usage{
read.matrix.hb(file)
model.matrix(object, \dots)
model.response(data,type)
}
\arguments{
  \item{file}{file name to read from or }
  \item{data, object}{an object of either 'matrix.csc.hb' or 'matrix.ssc.hb' class}
  \item{type}{One of `"any"', `"numeric"', `"double"'. Using the either of
          latter two coerces the result to have storage mode `"double"'}
  \item{\dots}{additional arguments to model.matrix}
}
\details{
Sparse coefficient matrices in the Harwell-Boeing format are stored in 
80-column records. Each file begins with a multiple line header block 
followed by two, three or four data blocks. The header block contains 
summary information on the storage formats and storage requirements. 
The data blocks contain information of the sparse coefficient matrix and 
data for the right-hand-side of the linear system of equations, 
initial guess of the solution and the exact solutions if they exist.
The function \code{model.matrix} extracts the X matrix component.
The function \code{model.response} extracts the y vector (or matrix).
The function \code{model.guess} extracts the guess vector. 
The function \code{model.xexact} extracts the xexact vector. 
This function is written in R replacing a prior implementation based
on iohb.c which had memory fault difficulties.  The function write.matrix.hb
has been purged; users wishing to write matrices in Harwell-Boeing format
are advised to convert SparseM matrices to Matrix classes and use writeHB
from the Matrix package.  Contributions of code to facilitate this conversion
would be appreciated!
}
\value{
The function \code{read.matrix.hb} returns a list of class 
\code{matrix.csc.hb} or \code{matrix.ssc.hb} depending
on how the coefficient matrix is stored in the \code{file}.
  \item{ra }{ra component of the csc or ssc format of the coefficient matrix, X.}
  \item{ja }{ja component of the csc or ssc format of the coefficient matrix, X.}
  \item{ia }{ia component of the csc or ssc format of the coefficient matrix, X.}
  \item{rhs.ra }{ra component of the right-hand-side, y, if stored in csc or 
	ssc format; right-hand-side stored in dense vector or matrix otherwise.}
  \item{rhs.ja }{ja component of the right-hand-side, y, if stored in csc or 
	ssc format; a null vector otherwise.}
  \item{rhs.ia }{ia component of the right-hand-side, y, if stored in csc or 
	ssc format; a null vector otherwise.}
  \item{xexact}{vector of the exact solutions, b, if they exist; a null vector otherwise.}
  \item{guess}{vector of the initial guess of the solutions if they exist; 
	a null vector otherwise.}
  \item{dimension}{dimenson of the coefficient matrix, X.}
  \item{rhs.dim}{dimenson of the right-hand-side, y.}
  \item{rhs.mode}{storage mode of the right-hand-side; can be full storage or
	same format as the coefficient matrix, for the moment the only allowed
	mode is "F" for full, or dense mode.}

The function \code{model.matrix} returns the X matrix of class \code{matrix.csr}.
The function \code{model.response} returns the y vector (or matrix).
The function \code{model.guess} returns the guess  vector (or matrix).
The function \code{model.xexact} returns the xexact vector (or matrix).
}
\references{ 
Duff, I.S., Grimes, R.G. and Lewis, J.G. (1992) 
User's Guide for Harwell-Boeing Sparse Matrix Collection at 
\url{http://math.nist.gov/MatrixMarket/collections/hb.html}}
\author{Pin Ng}
\seealso{
\code{slm} for sparse version of \code{lm} \cr
\code{SparseM.ops} for operators on class \code{matrix.csr} \cr
\code{SparseM.solve} for linear equation solving for class \code{matrix.csr} \cr
\code{SparseM.image} for image plotting of class \code{matrix.csr} \cr
\code{SparseM.ontology} for coercion of class \code{matrix.csr} \cr

}

\examples{ 
Xy <- read.matrix.hb(system.file("extdata","lsq.rra",package = "SparseM"))
class(Xy) # -> [1] "matrix.csc.hb"
X <- model.matrix(Xy)->X
class(X) # -> "matrix.csr"
dim(X) # -> [1] 1850  712
y <- model.response(Xy) # extract the rhs
length(y) # [1] 1850
Xy <- read.matrix.hb(system.file("extdata","rua_32_ax.rua",package = "SparseM"))
X <- model.matrix(Xy)
y <- model.response(Xy) # extract the rhs
g <- model.guess(Xy) # extract the guess 
a <- model.xexact(Xy) # extract the xexact 
fit <- solve(t(X) \%*\% X, t(X) \%*\% y) # compare solution with xexact solution 
}
\keyword{IO}
