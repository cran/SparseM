\name{SparseM.hb}
\alias{SparseM.hb}
\alias{read.matrix.hb}
\alias{write.matrix.hb}
\alias{model.matrix,ANY-method}
\alias{model.matrix,matrix.csc.hb-method}
\alias{model.matrix,matrix.ssc.hb-method}
\alias{model.matrix.matrix.ssc.hb}
\alias{model.response,ANY-method}
\alias{model.response,matrix.csc.hb-method}
\alias{model.response,matrix.ssc.hb-method}
\alias{model.response}
\alias{model.matrix}
\title{Harwell-Boeing Format Sparse Matrices}
\description{
Read, write and extract components of data in Harwell-Boeing  
sparse matrix format.
}
\usage{
read.matrix.hb(filename)
write.matrix.hb(filename = "hb.out", X, title, key, mxtype, rhs = NULL, 
guess = FALSE,  xsol = FALSE, ptrfmt = "(16I5)", indfmt = "(16I5)",
valfmt = "(5D16.9)", rhsfmt = "(5D16.9)")
model.matrix(object, \dots)
model.response(data,type)
}
\arguments{
  \item{filename}{file name to read from or write to }
  \item{data, object}{an object of either \code{matrix.csc.hb} or \code{matrix.ssc.hb} 
class}
  \item{type}{One of `"any"', `"numeric"', `"double"'. Using the either of
          latter two coerces the result to have storage mode
          `"double"'}
  \item{X}{coefficient matrix stored in csc (for unsymmetric and rectangular
        matrix) or ssc (for symmetric matrix) format}
  \item{title}{72-character title for the matrix}
  \item{key}{8-character identifier for the matrix}
  \item{mxtype}{3-character identifier for type of the coefficient matrix;\cr
        First character: currently only "R" for real matrix is supported;
        second character: "S" for symmetric, "U" for unsymmetric" 
		and "R" for rectangular matrix;
        third character: currently only "A" for assembled matrix is supported }
  \item{rhs}{vector or matrix of right-hand-side(s) including starting guesses
        and solution vectors if present stored in full storage mode}
  \item{guess}{logical flag for the presence of initial guess of the solutions;
        if TRUE, the values of initial guess are appended to the end of
        \code{rhs}}
  \item{xsol}{logical flag for the presence of exact solutions; if TRUE, the
        values of the exact solutions are appended to the end of \code{rhs}}
  \item{ptrfmt}{printing format for the column pointers}
  \item{indfmt}{printing format for the row indices}
  \item{valfmt}{printing format for the values}
  \item{rhsfmt}{printing format for the right-hand-sides}
  \item{\dots}{additional arguments}
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
}
\value{
The function \code{read.matrix.hb} returns a list of class \code{matrix.csc.hb} or \code{matrix.ssc.hb} depending
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
	same format as the coefficient matrix.}

The function \code{model.matrix} returns the X matrix of class \code{matrix.csr}.
The function \code{model.response} returns the y vector (or matrix).
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
\dontrun{read.matrix.hb(system.file("HBdata","lsq.rra",package = "SparseM"))-> hb.o
class(hb.o) # -> [1] "matrix.csc.hb"
model.matrix(hb.o)->X
class(X) # -> "matrix.csr"
dim(X) # -> [1] 1850  712
y <- model.response(hb.o) # extract the rhs
length(y) # [1] 1850
}
}
\keyword{IO}
