\name{slm}
\alias{slm}
\title{Fit a linear regression model using sparse matrix algebra}
\description{
This is a function to illustrate the use of sparse linear algebra
to solve a linear least squares problem using Cholesky decomposition.
The syntax and output attempt to emulate \code{lm()} but may
fail to do so fully satisfactorily.  Ideally, this would eventually
become a method for \code{lm}.  The main obstacle to this step is
that it would be necessary to have a model.matrix function that
returned an object in sparse csr form.  For the present, the objects
represented in the formula must be in dense form.  If the user wishes
to specify fitting with a design matrix that is already in sparse form,
then the lower level function \code{slm.fit()} should be used. 
}
\usage{
slm(formula, data, weights, na.action, method = "csr", contrasts = NULL, ...)
}
\arguments{
\item{formula}{
    a formula object, with the response on the left of a \code{~} operator,
    and the terms, separated by \code{+} operators, on the right.  As in
    \code{lm()}, the response variable in the formula can be matrix valued.
  }
  \item{data}{
    a data.frame in which to interpret the variables
    named in the formula, or in the subset and the weights argument.
    If this is missing, then the variables in the formula should be on the
    search list.  This may also be a single number to handle some special
    cases -- see below for details.
  }
\item{weights}{
    vector of observation weights; if supplied, the algorithm fits
    to minimize the sum of the weights multiplied into the
    absolute residuals. The length of weights must be the same as
    the number of observations.  The weights must be nonnegative
    and it is strongly recommended that they be strictly positive,
    since zero weights are ambiguous.
  }
  \item{na.action}{
    a function to filter missing data.
    This is applied to the model.frame after any subset argument has been used.
    The default (with \code{na.fail}) is to create an error if any missing 
    values are found.  A possible alternative is \code{na.omit}, which
    deletes observations that contain one or more missing values.
  }
\item{method}{there is only one method based on Cholesky factorization}
\item{contrasts}{
    a list giving contrasts for some or all of the factors
    default = \code{NULL} appearing in the model formula.
    The elements of the list should have the same name as the variable
    and should be either a contrast matrix (specifically, any full-rank
    matrix with as many rows as there are levels in the factor),
    or else a function to compute such a matrix given the number of levels.
  }
\item{...}{
    additional arguments for the fitting routines
  }
}
\value{
A list of class \code{slm} consisting of:
  \item{coefficients}{estimated coefficients}
  \item{chol}{cholesky object from fitting}
  \item{residuals}{residuals}
  \item{fitted}{fitted values}
  \item{terms}{terms}
  \item{call}{call}
  ...
}
\references{
Koenker, R and Ng, P. (2002).  SparseM:  A Sparse Matrix Package for \R, \cr
\url{http://www.econ.uiuc.edu/~roger/research}
}
\author{ Roger Koenker }

\seealso{
\code{slm.methods} for methods \code{summary}, \code{print}, \code{fitted},
\code{residuals} and \code{coef} associated with class \code{slm},
and \code{slm.fit} for lower level fitting functions.  The latter functions
are of special interest if you would like to pass a sparse form of the
design matrix directly to the fitting process.}
\examples{
data(lsq)
X <- model.matrix(lsq) #extract the design matrix
y <- model.response(lsq) # extract the rhs
X1 <- as.matrix(X)
slm.time <- unix.time(slm(y~X1-1) -> slm.o) # pretty fast
lm.time <- unix.time(lm(y~X1-1) -> lm.o) # very slow
cat("slm time =",slm.time,"\n")
cat("slm Results: Reported Coefficients Truncated to 5  ","\n")
sum.slm <- summary(slm.o)
sum.slm$coef <- sum.slm$coef[1:5,]
sum.slm
cat("lm time =",lm.time,"\n")
cat("lm Results: Reported Coefficients Truncated to 5  ","\n")
sum.lm <- summary(lm.o)
sum.lm$coef <- sum.lm$coef[1:5,]
sum.lm
}

\keyword{regression}
