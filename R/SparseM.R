#--------------------------------------------------------------------
".First.lib" <-
function(lib, pkg) {
   library.dynam("SparseM", pkg, lib)
   print("SparseM library loaded")}
#--------------------------------------------------------------------
"is.matrix.csr" <- function(x, ...) inherits(x,"matrix.csr")
#--------------------------------------------------------------------
"is.matrix.csc" <- function(x, ...) inherits(x,"matrix.csc")
#--------------------------------------------------------------------
"is.matrix.ssr" <- function(x, ...) inherits(x,"matrix.ssr")
#--------------------------------------------------------------------
"is.matrix.ssc" <- function(x, ...) inherits(x,"matrix.ssc")
#--------------------------------------------------------------------
"as.matrix.csr" <- function(x, ...) UseMethod("as.matrix.csr")
#--------------------------------------------------------------------
"as.matrix.csc" <- function(x, ...) UseMethod("as.matrix.csc")
#--------------------------------------------------------------------
"as.matrix.ssr" <- function(x, ...) UseMethod("as.matrix.ssr")
#--------------------------------------------------------------------
"as.matrix.ssc" <- function(x, ...) UseMethod("as.matrix.ssc")
#--------------------------------------------------------------------
"as.matrix.csr.default" <- function (x, ...) 
{
    if (is.matrix.csr(x)) 
        x
    else matrix.csr(x, ...)
}
#--------------------------------------------------------------------
"as.matrix.csr.matrix.csc" <- function(x, ...)
{
	x <- t(x)
	x$dim <- rev(dim(x))
	class(x) <- "matrix.csr"
	x
}
#--------------------------------------------------------------------
"as.matrix.csr.matrix.ssr" <- function(x, ...)
{
	.ssr.csr(x)
}
#--------------------------------------------------------------------
"as.matrix.csr.matrix.ssc" <- function(x, ...)
{
	.ssr.csr(x)
}
#--------------------------------------------------------------------
"as.matrix.matrix.csr" <-
function(x){
	nrow <- x$dim[1]
	ncol <- x$dim[2]
	if(length(x$ra)==1 && is.na(x$ra)){ #trap zero matrix
		dns <- matrix(0,nrow=nrow,ncol=ncol)
		return(dns)
		}
	nan <- is.nan(x$ra)
        infty <- is.infinite(x$ra) & x$ra >0
        ninfty <- is.infinite(x$ra) & x$ra <0
        uniq <- rnorm(3)
        while(any(uniq %in% x$ra[!(nan|infty|ninfty)]))
                uniq <- rnorm(3)
	x$ra[nan] <- uniq[1]
        x$ra[infty] <- uniq[2]
        x$ra[ninfty] <- uniq[3]
	z <- .Fortran("csrdns",
		as.integer(nrow),
		as.integer(ncol),
		as.double(x$ra),
		as.integer(x$ja),
		as.integer(x$ia),
		dns = double(nrow*ncol),
		ndns = as.integer(nrow),
		ierr = integer(1),
		PACKAGE = "SparseM")
	if(z$ierr != 0) stop("insufficient space for dns")
	dns <- matrix(z$dns,nrow=nrow,ncol=ncol)
	dns[dns==uniq[1]] <- NaN
	dns[dns==uniq[2]] <- Inf
	dns[dns==uniq[3]] <- -Inf
	return(dns)
}
#--------------------------------------------------------------------
"as.matrix.csc.matrix.csr" <- function(x, ...)
{
	x <- t(x)
	x$dim <- rev(dim(x))
	class(x) <- "matrix.csc"
	x
}
#--------------------------------------------------------------------
"as.matrix.ssr.matrix.csr" <- function(x, ...)
{
	.csr.ssr(x)
}
#--------------------------------------------------------------------
"as.matrix.ssc.matrix.csr" <- function(x, ...)
{
	x <- as.matrix.csc(x)
	x <- .csc.ssc(x)
	class(x) <- "matrix.ssc"
	x
}
#--------------------------------------------------------------------
"as.matrix.csc.default" <- function(x, ...)
{
    if (is.matrix.csc(x)) x
    else as.matrix.csc(as.matrix.csr(x, ...))
}
#--------------------------------------------------------------------
"as.matrix.csc.matrix.ssr" <- function(x, ...)
{
	as.matrix.csc(as.matrix.csr(x))
}
#--------------------------------------------------------------------
"as.matrix.csc.matrix.ssc" <- function(x, ...)
{
	as.matrix.csc(as.matrix.csr(x))
}
#--------------------------------------------------------------------
"as.matrix.ssr.default" <- function(x, ...)
{
	if (is.matrix.ssr(x)) x
	else as.matrix.ssr(as.matrix.csr(x, ...))
}
#--------------------------------------------------------------------
"as.matrix.ssr.matrix.csc" <- function(x, ...)
{
	as.matrix.ssr(as.matrix.csr(x))
}
#--------------------------------------------------------------------
"as.matrix.ssr.matrix.ssc" <- function(x, ...)
{
	as.matrix.ssr(as.matrix.csr(x))
}
#--------------------------------------------------------------------
"as.matrix.ssc.default" <- function(x, ...)
{
	if (is.matrix.ssc(x)) x
	else as.matrix.ssc(as.matrix.csc(x, ...))
}
#--------------------------------------------------------------------
"as.matrix.ssc.matrix.csc" <- function(x, ...)
{
	.csc.ssc(x)
}
#--------------------------------------------------------------------
"as.matrix.ssc.matrix.ssr" <- function(x, ...)
{
	as.matrix.ssc(as.matrix.csr(x))
}
#--------------------------------------------------------------------
"matrix.csr" <-
function(x, nrow = 1, ncol = 1, eps = .Machine$double.eps, ...){
         if (!is.matrix(x)) {
		x <-return(as.matrix.csr(matrix(x,nrow,ncol)))
		}
	dimx <- dim(x)
	nnz <- sum(abs(x)>eps)
	if(nnz==0){
	        z<-list(ra=0,ja=1,ia=c(1,rep(2,dimx[1])),dim=dimx)
		class(z) <- "matrix.csr"
		return(z)
        }
	z <- .Fortran("csr",
		as.double(x),
		ra=double(nnz),
		ja=integer(nnz),
		ia=integer(dimx[1]+1),
		as.integer(dimx[1]),
		as.integer(dimx[2]),
		nnz=as.integer(nnz),
		as.double(eps),
		PACKAGE = "SparseM")
	if(nnz!=z$nnz)warning("nnz values inconsistent")
	nnz <- z$nnz
	z <- list(ra = z$ra[1:nnz], ja = z$ja[1:nnz], ia = z$ia, dim = dimx)
	class(z) <- "matrix.csr"
	return(z)
}
#--------------------------------------------------------------------
"dim.matrix.csr" <-
function(x){x$dim}
#--------------------------------------------------------------------
"ncol.matrix.csr" <-
function(x){dim(x)[2]}
#--------------------------------------------------------------------
"nrow.matrix.csr" <-
function(x){dim(x)[1]}
#--------------------------------------------------------------------
"t.matrix.csr" <-
function(x){
# transpose a sparse matrix stored in csr format; a.k.a. transform the matrix
# stored in csr format to csc format
	nrow <- x$dim[1]
	ncol <- x$dim[2]
	nnz <- x$ia[nrow+1]-1
	z <- .Fortran("csrcsc2",
		as.integer(nrow),
		as.integer(ncol),
		as.integer(1),
		as.integer(1),
		as.double(x$ra),
		as.integer(x$ja),
		as.integer(x$ia),
		ao=double(nnz),
		jao=integer(nnz),
		iao=integer(ncol+1),
		PACKAGE = "SparseM")
	z <- list(ra = z$ao, ja = z$jao, ia = z$iao, dim = rev(x$dim))
	class(z) <- "matrix.csr"
	return(z)
}
#--------------------------------------------------------------------
"as.matrix.matrix.csc" <-function (x) 
{
	as.matrix(as.matrix.csr(x))
}
#--------------------------------------------------------------------
"as.matrix.matrix.ssr" <-function (x) 
{
	as.matrix(as.matrix.csr(x))
}
#--------------------------------------------------------------------
"as.matrix.matrix.ssc" <-function (x) 
{
	as.matrix(as.matrix.csr(x))
}
#--------------------------------------------------------------------
"dim.matrix.csc" <-
function (x) 
{
    x$dim
}
#--------------------------------------------------------------------
"dim.matrix.ssr" <-
function (x) 
{
    x$dim
}
#--------------------------------------------------------------------
"dim.matrix.ssc" <-
function (x) 
{
    x$dim
}
#--------------------------------------------------------------------
"t.matrix.csc" <- function (x) 
{
    nrow <- x$dim[1]
    ncol <- x$dim[2]
    nnz <- x$ia[ncol + 1] - 1
    z <- .Fortran("csrcsc2", as.integer(ncol), as.integer(nrow), 
        as.integer(1), as.integer(1), as.double(x$ra), as.integer(x$ja), 
        as.integer(x$ia), ao = double(nnz), jao = integer(nnz), 
        iao = integer(nrow + 1), PACKAGE = "SparseM")
    z <- list(ra = z$ao, ja = z$jao, ia = z$iao, dim = rev(x$dim))
    class(z) <- "matrix.csc"
    return(z)
}
#--------------------------------------------------------------------
"model.matrix.hb" <- function(object, ...) UseMethod("model.matrix.hb")
#--------------------------------------------------------------------
"model.matrix.hb.default" <- function(object, ...){
	object
}
#--------------------------------------------------------------------
"model.matrix.hb.matrix.csc.hb" <- function(object, ...){
	object <- list(ra=object$ra,ja=object$ja,ia=object$ia,dim=object$dim)
	class(object) <- "matrix.csc"
	object <- as.matrix.csr(object)
	object
}
#--------------------------------------------------------------------
"model.matrix.hb.matrix.ssc.hb" <- function(object, ...){
	object <- .ssc.csc(object)
	class(object) <- "matrix.csc"
	object <- as.matrix.csr(object)
	object
}
#--------------------------------------------------------------------
"model.response.hb" <- function(x){
	if(is.null(x$rhs.mode)) stop("Right-hand side doesn't exist")
	if (x$rhs.mode == "F")
		z <- x$rhs.ra
	else{
		z <- list(ra=x$rhs.ra,ja=x$rhs.ja,ia=x$rhs.ia,dim=x$rhs.dim)
		class(z) <- "matrix.csc"
		z <- as.matrix.csr(z)
		}
	z
}
#--------------------------------------------------------------------
".ssr.csr" <- function(x){
	nrow <- x$dim[1]
	nnza <- x$ia[nrow+1]-1
	nnzao <- 2*nnza #can be set smaller
	z <- .Fortran("ssrcsr",
		job = as.integer(0),
		value2 = as.integer(1),
		nrow = as.integer(nrow),
		a = as.double(x$ra),
		ja = as.integer(x$ja),
		ia = as.integer(x$ia),
		nzmax = as.integer(nnzao),
		ao = double(nnzao),
		jao = integer(nnzao),
		iao = integer(nrow+1),
		indu = integer(nrow),
		iwk = integer(nrow+1),
		ierr = integer(1),
		PACKAGE = "SparseM")
	if(z$ierr != 0) stop("Not enough space")		
	nnz <- z$iao[nrow+1]-1
	z <- list(ra=z$ao[1:nnz],ja=z$jao[1:nnz],ia=z$iao,dim=x$dim)
	class(z) <- "matrix.csr"
	return(z)
}
#--------------------------------------------------------------------
".ssc.csc" <- function(x){
	z <- .ssr.csr(x)
	class(z) <- "matrix.csc"
	return(z)
}
#--------------------------------------------------------------------
".csr.ssr" <- function(x){
	nrow <- x$dim[1]
	nnza <- ceiling((x$ia[nrow+1]-1)/2)+nrow
	if(nrow!=x$dim[2]) 
                stop("Cannot convert an asymmetric matrix into `matrix.ssr' class")

        if(sum(abs((t(as.matrix.csr(x))-as.matrix.csr(x))$ra))!=0)
                stop("Cannot convert an asymmetric matrix into `matrix.ssr' class")
	z <- .Fortran("csrssr",
		as.integer(nrow),
		as.double(x$ra),
		as.integer(x$ja),
		as.integer(x$ia),
		as.integer(nnza),
		ao = as.double(x$ra),
		jao = as.integer(x$ja),
		iao = as.integer(x$ia),
		ierr = integer(1),
		PACKAGE = "SparseM")
	if(z$ierr != 0) stop("Not enough space. This is usually caused by trying to convert an asymmetric matrix into ssr format")
	nnza <- z$iao[nrow+1]-1
	z <- list(ra=z$ao[1:nnza],ja=z$jao[1:nnza],ia=z$iao,dim=x$dim)
	class(z) <- "matrix.ssr"
	z
}
#--------------------------------------------------------------------
".csc.ssc" <- function(x){
	nrow <- x$dim[2]
	nnza <- ceiling((x$ia[nrow+1]-1)/2)+nrow
	if(nrow!=x$dim[1])
                stop("Cannot convert an asymmetric matrix into `matrix.ssc' class")
        if(sum(abs((t(as.matrix.csr(x))-as.matrix.csr(x))$ra))!=0)
                stop("Cannot convert an asymmetric matrix into `matrix.ssc' class")
	z <- .Fortran("cscssc",
		as.integer(nrow),
		as.double(x$ra),
		as.integer(x$ja),
		as.integer(x$ia),
		as.integer(nnza),
		ao = as.double(x$ra),
		jao = as.integer(x$ja),
		iao = as.integer(x$ia),
		ierr = integer(1),
		PACKAGE = "SparseM")
	if(z$ierr != 0) stop("Not enough space. This is usually caused by trying to convert an asymmetric matrix into ssc format")
	nnza <- z$iao[nrow+1]-1
	z <- list(ra=z$ao[1:nnza],ja=z$jao[1:nnza],ia=z$iao,dim=x$dim)
	class(z) <- "matrix.ssc"
	return(z)
}
#--------------------------------------------------------------------
image.matrix.csr <- function(x,col=c("white","gray"),xlab="column",ylab="row", ...){
#plot non-zero elements of a sparse matrix in csr format
n <- x$dim[1]
p <- x$dim[2]
z <- matrix(0,n,p)
column <- x$ja
row <- rep(n:1,diff(x$ia))
z[cbind(row,column)] <- 1
image(x=1:p,y=-(n:1),t(z),axes=FALSE, col=col,xlab=xlab,ylab=ylab)
axis(1,pretty(1:p))
axis(2,pretty(-(n:1)),labels=rev(pretty(1:n)))
box()
}
#--------------------------------------------------------------------
"image.matrix.csc" <- 
function (x, ...) 
{
	x <- as.matrix.csr(x)
	image(x)
}
#--------------------------------------------------------------------
"rbind.matrix.csr" <- function(...) {
# Very preliminary function to rbind matrix.csr objects no name handling
    allargs <- list(...)
    allargs <- allargs[sapply(allargs, length) > 0]
    n <- length(allargs)
    if (n == 0)
       stop("nothing to rbind")
    nms <- names(allargs)
    Ncol <- ncol(allargs[[1]])
    if(n>1){
      if(!all(sapply(allargs,is.matrix.csr)))stop("some args not csr format")
      if(!all(sapply(allargs,ncol) == Ncol))
        stop("args have differing numbers of columns")
        }
    if (is.null(nms))
        nms <- character(length(allargs))
    cl <- NULL
    perm <- rows <- rlabs <- vector("list", n)
    Nrow <- nia <- 0
    value <- clabs <- NULL
    all.levs <- list()
    ra <- ja <- ia <- dim <- NULL
    for(i in 1:n){
        xi <- allargs[[i]]
        ra <- c(ra,xi$ra)
        ja <- c(ja,xi$ja)
        ia <- c(ia[-(Nrow+1)],nia + xi$ia)
        nia <- ia[length(ia)]-1
        Nrow <- Nrow + length(xi$ia)-1
        }
z <- list(ra=ra, ja=ja, ia = ia, dim = c(Nrow,Ncol))
class(z) <- "matrix.csr"
return(z)
}
#--------------------------------------------------------------------
"cbind.matrix.csr" <- function(...)
{
# Very preliminary function to cbind matrix.csr objects no name handling
    allargs <- list(...)
    allargs <- allargs[sapply(allargs, length) > 0]
    n <- length(allargs)
    if (n == 0)
        return(structure(list(), class = "data.frame", row.names = character()))
    nms <- names(allargs)
    allargs <- sapply(allargs,t,simplify=FALSE)
    Ncol <- ncol(allargs[[1]])
    if(n>1){
      if(!all(sapply(allargs,is.matrix.csr)))stop("some args not csr format")
      if(!all(sapply(allargs,ncol) == Ncol))
        stop("args have differing numbers of rows")
        }
    if (is.null(nms))
        nms <- character(length(allargs))
    cl <- NULL
    perm <- rows <- rlabs <- vector("list", n)
    Nrow <- nia <- 0
    value <- clabs <- NULL
    all.levs <- list()
    ra <- ja <- ia <- dim <- NULL
for(i in 1:n){
        xi <- allargs[[i]]
        ra <- c(ra,xi$ra)
        ja <- c(ja,xi$ja)
        ia <- c(ia[-(Nrow+1)],nia + xi$ia)
        nia <- ia[length(ia)]-1
        Nrow <- Nrow + length(xi$ia)-1
        }
z <- list(ra=ra, ja=ja, ia = ia, dim = c(Nrow,Ncol))
class(z) <- "matrix.csr"
z <- t(z)
return(z)
}
#--------------------------------------------------------------------
"read.matrix.hb" <-
function (filename) 
{
	hb1.o <- .C("read_HB1", 
		infile = as.character(filename),
		M = integer(1),
		N = integer(1),
		nnz = integer(1),
		Nrhs = integer(1),
		mxtype = character(1),
		Rhstype = character(1),
		PACKAGE = "SparseM"
		)
	mxtype = hb1.o$mxtype
        if(substr(mxtype,1,1)!="R") stop("Doesn't handle non-real matrices")
        if(substr(mxtype,2,2) == "S")
                format <- "ssc"
        else if(substr(mxtype,2,2) == "R" | substr(mxtype,2,2) == "U")
                format <- "csc"
        else
                stop("Doesn't handle matrices other than symmetric or rectangular!")
        if(substr(mxtype,3,3)!="A") stop("Doesn't handle elemental matrices!")
	hb2.o <- .C("read_HB2",
		infile = as.character(filename),
		M = integer(1),
		N = integer(1),
		nnz = integer(1),
		integer(1),
		integer(1),
		double(1),
		colptr = integer(hb1.o$N+1),
		rowind = integer(hb1.o$nnz),
		val = double(hb1.o$nnz),
		PACKAGE = "SparseM"
		)
	rhs.mode <- NULL
	if(hb1.o$Nrhs > 0){
		rhs.mode <- "F"
		Gflag <- "No"
		Xflag <- "No"
		Rhstype <- hb1.o$Rhstype
		if(substr(Rhstype,1,1) != "F") 
			stop("Right-hand side has to be in full storage mode.")
		hb3.o <- .C("read_HB3",
			infile = as.character(filename),
			M = as.integer(hb1.o$M),
			Nrhs = as.integer(hb1.o$Nrhs),
			rhs = double(hb1.o$Nrhs*hb1.o$M),
			rhsflag = as.character("F"),
			PACKAGE = "SparseM"
			)
		if(substr(Rhstype,2,2) == "G"){
			hb4.o <- .C("read_HB3",
				infile = as.character(filename),
				M = as.integer(hb1.o$M),
				Nrhs = as.integer(hb1.o$Nrhs),
				rhs = double(hb1.o$Nrhs*hb1.o$M),
				rhsflag = as.character("G"),
				PACKAGE = "SparseM"
				)
			Gflag <- "Yes"
			}
		if(substr(Rhstype,3,3) == "X"){
			hb5.o <- .C("read_HB3",
				infile = as.character(filename),
				M = as.integer(hb1.o$M),
				Nrhs = as.integer(hb1.o$Nrhs),
				rhs = double(hb1.o$Nrhs*hb1.o$M),
				rhsflag = as.character("X"),
				PACKAGE = "SparseM"
				)
			Xflag <- "Yes"
			}
		rd.o <- list(ra = hb2.o$val, ja = hb2.o$rowind, ia = hb2.o$colptr, 
			rhs.ra = hb3.o$rhs, guess = switch(Gflag, Yes = hb4.o$rhs, No = NULL), 
			xexact = switch(Xflag, Yes = hb5.o$rhs, No = NULL), dim = c(hb1.o$M, hb1.o$N),
			rhs.dim = c(hb1.o$M, hb1.o$Nrhs), rhs.mode=rhs.mode)
		}
	else
		rd.o <- list(ra = hb2.o$val, ja = hb2.o$rowind, ia = hb2.o$colptr,
			dim=c(hb1.o$M,hb1.o$N),rhs.mode=rhs.mode)
        if(format=="csc")
                class(rd.o) <- c("matrix.csc.hb")
        else
                class(rd.o) <- c("matrix.ssc.hb")

   return(rd.o)
}
#--------------------------------------------------------------------
"write.matrix.hb" <- function(filename="hb.out",X,title,key,mxtype,rhs=NULL,
		guess=FALSE,xsol=FALSE,ptrfmt="(16I5)",indfmt="(16I5)",
		valfmt="(1P,5D16.9)", rhsfmt="(1P,5D16.9)"){
        if(!substr(mxtype,1,1)%in%c("r","R")) stop("The first character of `mxtype' can only be 'R'")
        if(!substr(mxtype,2,2)%in%c("s","S","u","U","r","R")) stop("The second character of `mxtype' can only be `S',`U' or 'R'")
        if(!substr(mxtype,3,3)%in%c("a","A")) stop("The third character of `mxtype' can only be `A'")
        if(substr(mxtype,2,2)%in%c("s","S") && !is.matrix.ssc(X)) stop("Matrix X has to be in in ssc format")
        if(substr(mxtype,2,2)%in%c("u","U","r","R") && !is.matrix.csc(X)) stop("Matrix X has to be in in csc format")
	M <- X$dim[1]
	N <- X$dim[2]
	nnz <- length(X$ra)
	nrhs <- 0
	guesol <- NULL
	Rhs <- Guess <- Exact <- rep(0,M)
        if(!missing(rhs)){
		guesol <- paste(guesol,"F",sep="")
                idiv <- 1
		Rhs <- rhs[1:(M*idiv)]
                if(guess){
                        idiv <- idiv + 1
                        guesol <- paste(guesol,"G",sep="")
			Guess <- rhs[(M*(idiv-1)+1):(M*idiv)]
                        }
                if(xsol){
                        idiv <- idiv+1
                        guesol <- paste(guesol,"X",sep="")
			Exact <- rhs[(M*(idiv-1)+1):(M*idiv)]
                        }
                nrhs <- length(rhs)/M/idiv
                if(length(rhs)%%(M*idiv) != 0) stop("The length of `rhs' is not a multiple of the number of equations")
                }
	.C("write_HB1",
		as.character(filename),
		as.integer(M),
		as.integer(N),
		as.integer(nnz),
		as.integer(X$ia),
		as.integer(X$ja),
		as.double(X$ra),
		as.integer(nrhs),
		as.double(Rhs),
		as.double(Guess),
		as.double(Exact),
		as.character(title),
		as.character(key),
		as.character(mxtype),
		as.character(guesol),
		as.character(ptrfmt),
		as.character(indfmt),
		as.character(valfmt),
		as.character(rhsfmt),
		PACKAGE = "SparseM"
		)
}
#--------------------------------------------------------------------
"Ops.matrix.csr" <- function(e1,e2){
	if(missing(e2)){
		e1.op <- switch(.Generic,
			"+" = e1,
			"-" = structure(list(ra=-e1$ra,ja=e1$ja,ia=e1$ia,dim=e1$dim),class="matrix.csr"),
			"!" = .matrix.csr.compl(e1),
			stop(paste("Unary operator \"",.Generic,"\""," is undefined for class \"matrix.csr\"",sep=""))
			)
	return(e1.op)
		}
	e1.op.e2 <- {
		switch(.Generic,
		"+" = .matrix.csr.addsub(e1,e2,1),
		"-" = .matrix.csr.addsub(e1,e2,-1),
		"*" = .matrix.csr.elemul(e1,e2),
		"/" = .matrix.csr.elediv(e1,e2),
		"^" = .matrix.csr.expo(e1,e2),
		"%/%" = .matrix.csr.intdiv(e1,e2),
		"%%" = .matrix.csr.mod(e1,e2),
		">" = .matrix.csr.relation(e1,e2,"gt"),
		">=" = .matrix.csr.relation(e1,e2,"ge"),
		"<" = .matrix.csr.relation(e1,e2,"lt"),
		"<=" = .matrix.csr.relation(e1,e2,"le"),
		"==" = .matrix.csr.relation(e1,e2,"eq"),
		"!=" = .matrix.csr.relation(e1,e2,"ne"),
		"&" = {z <- .matrix.csr.elemul(e1,e2);z$ra <- rep(1,length(z$ja));z},
		"|" = {z <- .matrix.csr.addsub(e1,e2,1);z$ra <- rep(1,length(z$ja));z},
		stop(paste("Binary operator \"",.Generic,"\""," is undefined for class \"matrix.csr\"",sep=""))
		)
		}
	return(e1.op.e2)
}
#--------------------------------------------------------------------
".matrix.csr.compl" <- function(e1){
	nrow <- e1$dim[1]
	ncol <- e1$dim[2]
	nnz <- e1$ia[nrow+1]-1
	nz <- nrow*ncol - nnz
	if(length(e1$ra) == 1 && is.na(e1$ra)){ #trap zero matrix
		z <- list(ra=rep(1,nz),ja=rep(1:ncol,nrow),ia=seq(1,nz+1,by=ncol),dim=e1$dim)
		}
	else{
		z <- .Fortran("nzero",
			as.double(e1$ra),
			as.integer(e1$ja),
			as.integer(e1$ia),
			as.integer(nrow),
			as.integer(ncol),
			as.integer(nnz),
			as.integer(nz),
			ra = double(nz),
			ja = integer(nz),
			ia = integer(nrow+1),
			logical(ncol),
			PACKAGE = "SparseM")
		z <- list(ra=z$ra,ja=z$ja,ia=z$ia,dim=e1$dim)
		}
	class(z) <- "matrix.csr"
	z
}
#--------------------------------------------------------------------
".matrix.csr.addsub" <-
function(A,B,s){
#matrix addition/subtraction of two sparse csr matrices 
nrow <- A$dim[1]
ncol <- A$dim[2]
Bcol <- B$dim[2]
Brow <- B$dim[1]
if(ncol != Bcol | nrow != Brow)stop("matrices not conformable for addition")
nnza <- A$ia[nrow+1]-1
nnzb <- B$ia[nrow+1]-1
z <- .Fortran("csort",
	as.integer(nrow),
	ra = as.double(A$ra),
	ja = as.integer(A$ja),
	ia = as.integer(A$ia),
	integer(max(nrow+1,2*nnza)),
	as.logical(TRUE),
	PACKAGE = "SparseM")
A$ra <- z$ra
A$ja <- z$ja
A$ia <- z$ia
z <- .Fortran("csort",
	as.integer(nrow),
	ra = as.double(B$ra),
	ja = as.integer(B$ja),
	ia = as.integer(B$ia),
	integer(max(nrow+1,2*nnzb)),
	as.logical(TRUE),
	PACKAGE = "SparseM")
B$ra <- z$ra
B$ja <- z$ja
B$ia <- z$ia
nnzmax <- length(union(A$ja+A$dim[2]*(rep(1:A$dim[1],diff(A$ia))-1),
	B$ja+B$dim[2]*(rep(1:B$dim[1],diff(B$ia))-1)))+1
z <- .Fortran("aplsb",
	as.integer(nrow),
	as.integer(ncol),
	as.double(A$ra),
	as.integer(A$ja),
	as.integer(A$ia),
	as.double(s),
	as.double(B$ra),
	as.integer(B$ja),
	as.integer(B$ia),
	ra = double(nnzmax),
	ja = integer(nnzmax),
	ia = integer(nrow+1),
	as.integer(nnzmax),
	ierr = integer(1),
	PACKAGE = "SparseM")
if(z$ierr != 0) stop("insufficient space for sparse matrix addition")
nnz <- z$ia[nrow+1]-1
z <- list(ra=z$ra[1:nnz],ja=z$ja[1:nnz],ia=z$ia,dim=c(nrow,ncol))
class(z) <- "matrix.csr"
return(z)
}
#--------------------------------------------------------------------
".matrix.csr.elemul" <-
function(A,B){
#element-wise matrix multiplication of two sparse csr matrices 
if(is.numeric(A) && length(A) == 1)
	z <- list(ra=A*B$ra,ja=B$ja,ia=B$ia,dim=B$dim)
else if(is.numeric(B) && length(B) == 1)
	z <- list(ra=B*A$ra,ja=A$ja,ia=A$ia,dim=A$dim)
else if(is.matrix.csr(A) || is.matrix.csr(B) || is.matrix(A) || is.matrix(B)){
	if(is.matrix(A)) A <- as.matrix.csr(A)
	if(is.matrix(B)) B <- as.matrix.csr(B)
	nrow <- A$dim[1]
	ncol <- A$dim[2]
	Bcol <- B$dim[2]
	Brow <- B$dim[1]
	if(ncol != Bcol | nrow != Brow)stop("matrices not conformable for element-by-element multiplication")
	nnza <- A$ia[nrow+1]-1
	nnzb <- B$ia[nrow+1]-1
	z <- .Fortran("csort",
		as.integer(nrow),
		ra = as.double(A$ra),
		ja = as.integer(A$ja),
		ia = as.integer(A$ia),
		integer(max(nrow+1,2*nnza)),
		as.logical(TRUE),
		PACKAGE = "SparseM")
	A$ra <- z$ra
	A$ja <- z$ja
	A$ia <- z$ia
	z <- .Fortran("csort",
		as.integer(nrow),
		ra = as.double(B$ra),
		ja = as.integer(B$ja),
		ia = as.integer(B$ia),
		integer(max(nrow+1,2*nnzb)),
		as.logical(TRUE),
		PACKAGE = "SparseM")
	B$ra <- z$ra
	B$ja <- z$ja
	B$ia <- z$ia
	nnzmax <- length(intersect(A$ja+A$dim[2]*(rep(1:A$dim[1],diff(A$ia))-1),
		B$ja+B$dim[2]*(rep(1:B$dim[1],diff(B$ia))-1)))+1
	z <- .Fortran("aemub",
		as.integer(nrow),
		as.integer(ncol),
		as.double(A$ra),
		as.integer(A$ja),
		as.integer(A$ia),
		as.double(B$ra),
		as.integer(B$ja),
		as.integer(B$ia),
		ra = double(nnzmax),
		ja = integer(nnzmax),
		ia = integer(nrow+1),
		as.integer(nnzmax),
		ierr = integer(1),
		PACKAGE = "SparseM")
	if(z$ierr != 0) stop("insufficient space for element-wise sparse matrix multiplication")
	nnz <- z$ia[nrow+1]-1
	z <- list(ra=z$ra[1:nnz],ja=z$ja[1:nnz],ia=z$ia,dim=c(nrow,ncol))
	}
else stop("Arguments have to be class \"matrix.csr\" or numeric")
class(z) <- "matrix.csr"
return(z)
}
#--------------------------------------------------------------------
".matrix.csr.elediv" <- function(A,B){
# Element-wise matrix division of two sparse csr matrices 
# This operation is not efficient storage-wise for sparse matrices
if(is.numeric(A) && length(A) == 1)
	z <- list(ra=A/B$ra,ja=B$ja,ia=B$ia,dim=B$dim)
else if(is.numeric(B) && length(B) == 1)
	z <- list(ra=B/A$ra,ja=A$ja,ia=A$ia,dim=A$dim)
else if(is.matrix.csr(A) || is.matrix.csr(B) || is.matrix(A) || is.matrix(B)){
	if(is.matrix(A)) A <- as.matrix.csr(A)
	if(is.matrix(B)) B <- as.matrix.csr(B)
	nrow <- A$dim[1]
	ncol <- A$dim[2]
	Bcol <- B$dim[2]
	Brow <- B$dim[1]
	if(ncol != Bcol | nrow != Brow)stop("matrices not conformable for element-by-element division")
	nnza <- A$ia[nrow+1]-1
	nnzb <- B$ia[nrow+1]-1
	z <- .Fortran("csort",
		as.integer(nrow),
		ra = as.double(A$ra),
		ja = as.integer(A$ja),
		ia = as.integer(A$ia),
		integer(max(nrow+1,2*nnza)),
		as.logical(TRUE),
		PACKAGE = "SparseM")
	A$ra <- z$ra
	A$ja <- z$ja
	A$ia <- z$ia
	z <- .Fortran("csort",
		as.integer(nrow),
		ra = as.double(B$ra),
		ja = as.integer(B$ja),
		ia = as.integer(B$ia),
		integer(max(nrow+1,2*nnzb)),
		as.logical(TRUE),
		PACKAGE = "SparseM")
	B$ra <- z$ra
	B$ja <- z$ja
	B$ia <- z$ia
	nnzmax <- length(union(A$ja+A$dim[2]*(rep(1:A$dim[1],diff(A$ia))-1),
		B$ja+B$dim[2]*(rep(1:B$dim[1],diff(B$ia))-1)))+1
	z <- .Fortran("aedib",
		as.integer(nrow),
		as.integer(ncol),
		as.double(A$ra),
		as.integer(A$ja),
		as.integer(A$ia),
		as.double(B$ra),
		as.integer(B$ja),
		as.integer(B$ia),
		ra = double(nnzmax),
		ja = integer(nnzmax),
		ia = integer(nrow+1),
		as.integer(nnzmax),
		ierr = integer(1),
		PACKAGE = "SparseM")
	if(z$ierr != 0) stop("insufficient space for element-wise sparse matrix multiplication")
	nnz <- z$ia[nrow+1]-1
	z1 <- vector("numeric",nrow*ncol)
        idx1 <- z$ja[1:nnz]+ncol*(rep(1:nrow,diff(z$ia))-1)
        idx2 <- union(A$ja+A$dim[2]*(rep(1:A$dim[1],diff(A$ia))-1),
                B$ja+B$dim[2]*(rep(1:B$dim[1],diff(B$ia))-1))
        idx3 <- setdiff(1:(nrow*ncol),idx2)
	z1[idx1] <- z$ra[1:nnz]
	z1[idx3] <- NaN
	z <- list(ra=z1,ja=rep(1:ncol,nrow),ia=seq(1,nrow*ncol+1,by=ncol),dim=c(nrow,ncol))
	}
else stop("Arguments have to be class \"matrix.csr\" or numeric")
class(z) <- "matrix.csr"
return(z)
}
#--------------------------------------------------------------------
".matrix.csr.expo" <- function(A,B){
# This operation is not efficient storage-wise for sparse matrices
	if(is.matrix.csr(A))
		A <- as.matrix(A)
	if(is.matrix.csr(B))
		B <- as.matrix(B)
	AB <- A^B
	nan <- is.nan(AB)
        infty <- is.infinite(AB) & AB >0
        ninfty <- is.infinite(AB) & AB <0
        uniq <- rnorm(3)
        while(any(uniq %in% AB[!(nan|infty|ninfty)]))
                uniq <- rnorm(3)
	AB[nan] <- uniq[1]
        AB[infty] <- uniq[2]
        AB[ninfty] <- uniq[3]
        AB <- as.matrix.csr(AB)
        AB$ra[AB$ra==uniq[1]] <- NaN
        AB$ra[AB$ra==uniq[2]] <- Inf
        AB$ra[AB$ra==uniq[3]] <- -Inf
        class(AB) <- "matrix.csr"
	AB
}
#--------------------------------------------------------------------
".matrix.csr.intdiv" <- function(A,B){
# This operation is not efficient storage-wise for sparse matrices
	if(is.matrix.csr(A))
		A <- as.matrix(A)
	if(is.matrix.csr(B))
		B <- as.matrix(B)
	AB <- A%/%B
	nan <- is.nan(AB)
        infty <- is.infinite(AB) & AB >0
        ninfty <- is.infinite(AB) & AB <0
        uniq <- rnorm(3)
        while(any(uniq %in% AB[!(nan|infty|ninfty)]))
                uniq <- rnorm(3)
	AB[nan] <- uniq[1]
        AB[infty] <- uniq[2]
        AB[ninfty] <- uniq[3]
        AB <- as.matrix.csr(AB)
        AB$ra[AB$ra==uniq[1]] <- NaN
        AB$ra[AB$ra==uniq[2]] <- Inf
        AB$ra[AB$ra==uniq[3]] <- -Inf
        class(AB) <- "matrix.csr"
	AB
}
#--------------------------------------------------------------------
".matrix.csr.mod" <- function(A,B){
# This operation is not efficient storage-wise for sparse matrices
	if(is.matrix.csr(A))
		A <- as.matrix(A)
	if(is.matrix.csr(B))
		B <- as.matrix(B)
	AB <- A%%B
	nan <- is.nan(AB)
        infty <- is.infinite(AB) & AB >0
        ninfty <- is.infinite(AB) & AB <0
        uniq <- rnorm(3)
        while(any(uniq %in% AB[!(nan|infty|ninfty)]))
                uniq <- rnorm(3)
	AB[nan] <- uniq[1]
        AB[infty] <- uniq[2]
        AB[ninfty] <- uniq[3]
        AB <- as.matrix.csr(AB)
        AB$ra[AB$ra==uniq[1]] <- NaN
        AB$ra[AB$ra==uniq[2]] <- Inf
        AB$ra[AB$ra==uniq[3]] <- -Inf
        class(AB) <- "matrix.csr"
	AB
}
#--------------------------------------------------------------------
".matrix.csr.relation" <- function(e1,e2,rel){
	if(is.numeric(e2) && length(e2) == 1){
		z <- .csr.relation(e1,e2,rel)
		}
	else if(is.numeric(e1) && length(e1) == 1){
		z <- .csr.relation(e2,e1,rel)
		}
	else { #inefficient implementation
		if (is.matrix(e2)) e1 <- as.matrix(e1)
		else if (is.matrix(e1)) e2 <- as.matrix(e2)
		else {e1 <- as.matrix(e1); e2 <- as.matrix(e2)}
		z <- switch(rel,
			"gt" = e1 > e2,
			"lt" = e1 < e2,
			"ge" = e1 >= e2,
			"le" = e1 <= e2,
			"eq" = e1 == e2,
			"ne" = e1 != e2)
		z <- as.matrix.csr(z)
		}
	z$ra <- rep(1,length(z$ra))
	z
}
#--------------------------------------------------------------------
".csr.relation" <- function(A,drptol,rel){
	nrow <- A$dim[1]
	nnza <- A$ia[nrow+1]-1
	flag <- FALSE
	if(rel=="gt" && drptol >=0){
		relidx <- 1
		flag <- TRUE
		}
	if(rel=="ge" && drptol >=0){
		relidx <- 2
		flag <- TRUE
		}
	if(rel=="lt" && drptol <=0){
		relidx <- 1
		drptol <- -drptol
		A$ra <- -A$ra
		flag <- TRUE
		}
	if(rel=="le" && drptol <=0){
		relidx <- 2
		drptol <- -drptol
		A$ra <- -A$ra
		flag <- TRUE
		}
	if(rel=="eq" && drptol !=0){
		relidx <- 3
		flag <- TRUE
		}
	if(rel=="ne" && drptol ==0){
		relidx <- 4
		flag <- TRUE
		}
	if(flag){
		z <- .Fortran("filter1",
			as.integer(nrow),
			as.integer(relidx),
			as.double(drptol),
			as.double(A$ra),
			as.integer(A$ja),
			as.integer(A$ia),
			ra = double(nnza),
			ja = integer(nnza),
			ia = integer(nrow+1),
			as.integer(nnza),
			ierr = integer(1),
		PACKAGE = "SparseM")
		if(z$ierr !=0) stop("Not enough space")
		nnza <- z$ia[nrow+1]-1
		if(rel == "lt" || rel == "le")
			z <- list(ra=z$ra[1:nnza],ja=z$ja[1:nnza],ia=z$ia,dim=A$dim)
		else
			z <- list(ra=-z$ra[1:nnza],ja=z$ja[1:nnza],ia=z$ia,dim=A$dim)
		}
	else{ #This operation is inefficient storage-wise
		z <- as.matrix.csr(as.matrix(A) > drptol)
		}
	class(z) <- "matrix.csr"
	z
}
#--------------------------------------------------------------------
"[.matrix.csr" <- function(x,rw=1:x$dim[1],cl=1:x$dim[2]){
	sorted <- FALSE
	nrow <- nrow1 <- x$dim[1]
	ncol <- ncol1 <- x$dim[2]
	nnza <- x$ia[nrow+1]-1
	z <- .Fortran("csrcoo",
		as.integer(nrow),
		as.integer(1),
		as.integer(nnza),
		as.double(x$ra),
		as.integer(x$ja),
		as.integer(x$ia),
		nnz = integer(1),
		ao = as.double(x$ra),
		ir = integer(nnza),
		jc = as.integer(x$ja),
		ierr = integer(1),
		PACKAGE = "SparseM")
	xir <- z$ir
	xic <- z$jc
	if(z$ierr !=0) stop("Not enough space")
	if(is.matrix.csr(rw)){
		if(nrow!=rw$dim[1] || ncol!=rw$dim[2]) 
			stop("Indexing matrix has a different dimension than the matrix to be indexed")
		}
	else{
		if(!all(abs(rw)%in%1:nrow)||!all(abs(cl)%in%1:ncol)) 
			stop("Subscripts out of bound")
		}
	if(missing(cl))
		if(is.matrix(rw)){
			case <- "matidx"
			rwidx <- rw[,1]
			clidx <- rw[,2]
			if(any(rw < 0)){ #negative indexing
				if(!all(rw<=0)) stop("Only 0's may mix with negative subscripts")
				else{
					rw <- setdiff(1:(nrow*ncol),abs(rw[,1]+(rw[,2]-1)*nrow))
					clidx <- ceiling(rw/nrow)
					rwidx <- rw - floor((rw-1)/nrow)*nrow
					}
				}
			nsub <- length(rwidx)
			}
		else if(is.vector(rw)){
			case <- "rwclidx"
			if(any(rw < 0)){ #negative indexing
				if(!all(rw<=0)) stop("Only 0's may mix with negative subscripts")
				else
					rw <- setdiff(1:nrow,abs(rw))
				}
			lcl <- length(cl)
       			lrw <- nrow1 <- length(rw)
	                rwidx <- rep(rw,lcl)
       		        clidx <- rep(cl,rep(lrw,lcl))
			m.o <- match(paste(rwidx,clidx),paste(xir,xic))
			midx <- m.o[!is.na(m.o)]
			rwidx <- xir[sort(midx)]
			clidx <- clidx1 <- xic[sort(midx)]
			rwidx1 <- (1:lrw)[match(rwidx,rw)]
			nsub <- length(midx)
			}
		else if(is.matrix.csr(rw)){
			case <- "csridx"
			rw <- t(rw)
			nrw <- rw$dim[1]
			if(ncol != nrw) stop("dimension of the indexing matrix is not the same as the matrix being indexed")
			nnzb <- rw$ia[nrw+1]-1
			z <- .Fortran("csrcoo",
				as.integer(nrw),
				as.integer(1),
				as.integer(nnzb),
				as.double(rw$ra),
				as.integer(rw$ja),
				as.integer(rw$ia),
				nnz = integer(1),
				ao = double(nnzb),
				ir = integer(nnzb),
				jc = integer(nnzb),
				ierr = integer(1),
				PACKAGE = "SparseM")
			if(z$ierr !=0) stop("Not enough space")
			clidx <- z$ir
			rwidx <- rw$ja
			nsub <- length(rwidx)
			}
		else stop("Invalid indexing")
	else if(missing(rw))
		if(is.vector(cl)){
			case <- "rwclidx"
			if(any(cl < 0)){ #negative indexing
				if(!all(cl<=0)) stop("Only 0's may mix with negative subscripts")
				else
					cl <- setdiff(1:ncol,abs(cl))
				}
			lcl <- ncol1 <- length(cl)
       			lrw <- length(rw)
	                rwidx <- rep(rw,lcl)
       		        clidx <- rep(cl,rep(lrw,lcl))
                        m.o <- match(paste(rwidx,clidx),paste(xir,xic))
                        midx <- m.o[!is.na(m.o)]
                        clidx <- xic[midx]
                        rwidx <- rwidx1 <- xir[midx]
                        clidx1 <- (1:lcl)[match(clidx,cl)]
			nsub <- length(midx)
			}
		else stop("The second index has to be a vector")
	else if(missing(rw)&&missing(cl))
		case <- wholeidx
	else {
		case <- "rwclidx"
		if(any(rw < 0)){ #negative indexing
			if(!all(rw<=0)) stop("Only 0's may mix with negative subscripts")
			else
				rw <- setdiff(1:nrow,abs(rw))
			}
		if(any(cl < 0)){ #negative indexing
			if(!all(cl<=0)) stop("Only 0's may mix with negative subscripts")
			else
				cl <- setdiff(1:ncol,abs(cl))
			}
		lcl <- ncol1 <- length(cl)
		lrw <- nrow1 <- length(rw)
		rwidx <- rep(rw,lcl)
		clidx <- rep(cl,rep(lrw,lcl))
		m.o <- match(paste(rwidx,clidx),paste(xir,xic))
		midx <- m.o[!is.na(m.o)]
		rwidx <- xir[sort(midx)]
		clidx <- xic[sort(midx)]
		rwidx1 <- (1:lrw)[match(rwidx,rw)]
		clidx1 <- (1:lcl)[match(clidx,cl)]
		nsub <- length(midx)
		}
	z <- .Fortran("subext",
		as.integer(nsub),
		as.integer(rwidx),
		as.integer(clidx),
		as.double(x$ra),
		as.integer(x$ja),
		as.integer(x$ia),
		as.logical(sorted),
		values = double(nsub),
		iadd = integer(nsub),
		PACKAGE = "SparseM")
	if(case=="rwclidx"){
		z1 <- .Fortran("coocsr",
			as.integer(nrow1),
			as.integer(nsub),
			as.double(z$values),
			as.integer(rwidx1),
			as.integer(clidx1),
			ao = double(nsub),
			jao = integer(nsub),
			iao = integer(nrow1+1),
			PACKAGE = "SparseM")
		x.sub <- list(ra=z1$ao,ja=z1$jao,ia=z1$iao,dim=c(nrow1,ncol1))
		class(x.sub) <- "matrix.csr"
		}
	values <- switch(case,
		"rwclidx" = x.sub,
		"matidx" =, "vecidx" = z$values, 
		"csridx" = z$values,
		"wholeidx" = as.matrix(x)
		)
	return(values)
}
#--------------------------------------------------------------------
"[<-.matrix.csr" <- function(x,rw=1:x$dim[1],cl=1:x$dim[2],value){
	sorted <- FALSE
	nrow <- x$dim[1]
	ncol <- x$dim[2]
	nnza <- x$ia[nrow+1]-1
	if(is.matrix.csr(rw)){
                if(nrow!=rw$dim[1] || ncol!=rw$dim[2])
                        stop("Indexing matrix has a different dimension than the matrix to be indexed")
                }
        else{
                if(!all(abs(rw)%in%1:nrow)||!all(abs(cl)%in%1:ncol))
                        stop("Subscripts out of bound")
                }
	if(missing(cl))
		if(is.matrix(rw)){
			case <- "matidx"
			rwidx <- rw[,1]
			clidx <- rw[,2]
			if(any(rw < 0)){ #negative indexing
				if(!all(rw<=0)) stop("Only 0's may mix with negative subscripts")
				else{
					rw <- setdiff(1:(nrow*ncol),abs(rw[,1]+(rw[,2]-1)*nrow))
					clidx <- ceiling(rw/nrow)
					rwidx <- rw - floor((rw-1)/nrow)*nrow
					}
				}
			nsub <- length(rwidx)
			}
		else if(is.vector(rw)){
			case <- "rwclidx"
			if(any(rw < 0)){ #negative indexing
				if(!all(rw<=0)) stop("Only 0's may mix with negative subscripts")
				else
					rw <- setdiff(1:nrow,abs(rw))
				}
                        lcl <- length(cl)
                        lrw <- length(rw)
                        rwidx <- rep(rw,lcl)
                        clidx <- rep(cl,rep(lrw,lcl))
			nsub <- length(rwidx)
			}
		else if(is.matrix.csr(rw)){
			case <- "csridx"
			rw <- t(rw)
			nrw <- rw$dim[1]
			if(ncol != nrw) stop("dimension of the indexing matrix is not the same as the matrix being indexed")
			nnzb <- rw$ia[nrw+1]-1
			z <- .Fortran("csrcoo",
				as.integer(nrw),
				as.integer(1),
				as.integer(nnzb),
				as.double(rw$ra),
				as.integer(rw$ja),
				as.integer(rw$ia),
				nnz = integer(1),
				ao = double(nnzb),
				ir = integer(nnzb),
				jc = integer(nnzb),
				ierr = integer(1),
				PACKAGE = "SparseM")
			if(z$ierr !=0) stop("Not enough space")
			clidx <- z$ir
			rwidx <- rw$ja
			nsub <- length(rwidx)
			}
		else stop("Invalid indexing")
	else if(missing(rw))
                if(is.vector(cl)){
                        case <- "rwclidx"
                        if(any(cl < 0)){ #negative indexing
                                if(!all(cl<=0)) stop("Only 0's may mix with negative subscripts")
                                else
                                        cl <- setdiff(1:ncol,abs(cl))
                                }
                        lcl <- length(cl)
                        lrw <- length(rw)
                        rwidx <- rep(rw,lcl)
                        clidx <- rep(cl,rep(lrw,lcl))
                        nsub <- length(rwidx)
                        }
                else stop("The second index has to be a vector")
        else if(missing(rw)&&missing(cl))
                case <- wholeidx
	else {
		case <- "rwclidx"
		if(any(rw < 0)){ #negative indexing
			if(!all(rw<=0)) stop("Only 0's may mix with negative subscripts")
			else
				rw <- setdiff(1:nrow,abs(rw))
			}
		if(any(cl < 0)){ #negative indexing
			if(!all(cl<=0)) stop("Only 0's may mix with negative subscripts")
			else
				cl <- setdiff(1:ncol,abs(cl))
			}
		lcl <- length(cl)
		lrw <- length(rw)
		rwidx <- rep(rw,lcl)
		clidx <- rep(cl,rep(lrw,lcl))
		nsub <- length(rwidx)
		}
	z <- .Fortran("subext", # to identify the zero's to be replaced
                as.integer(nsub),
                as.integer(rwidx),
                as.integer(clidx),
                as.double(x$ra),
                as.integer(x$ja),
                as.integer(x$ia),
                as.logical(sorted),
                value = double(nsub),
                iadd = integer(nsub),
		PACKAGE = "SparseM")
	nadd <- sum(z$iadd==0)
	nnzb <- nnza + nadd
	z <- .Fortran("subasg",
		as.integer(nrow),
		as.integer(ncol),
		as.integer(nsub),
		as.integer(nnza),
		as.integer(nnzb),
		as.integer(rwidx),
		as.integer(clidx),
		as.double(x$ra),
		as.integer(x$ja),
		as.integer(x$ia),
		ra = double(nnzb),
		ja = integer(nnzb),
		ia = integer(nrow+1),
		value = as.double(value),
		logical(ncol),
		ierr = integer(1),
		PACKAGE = "SparseM")
        if(z$ierr != 0) stop("not enough space")
	z <- list(ra=z$ra,ja=z$ja,ia=z$ia,dim=x$dim)
	class(z) <- "matrix.csr"
	return(z)
}
#--------------------------------------------------------------------
"%*%" <-
function(e1, e2) UseMethod("%*%")
#--------------------------------------------------------------------
"%*%.default" <-
.Primitive("%*%")
#--------------------------------------------------------------------
"%*%.matrix.csr" <-
function(A,B,nnzmax){
if(is.matrix.csr(B)){
   #matrix multiply two sparse csr matrices 
   nrow <- A$dim[1]
   ncol <- B$dim[2]
   Acol <- A$dim[2]
   Brow <- B$dim[1]
   if(Acol != Brow)
	stop("matrices not conformable for multiplication")
   if(missing(nnzmax)){
	z <- .Fortran("amubdg",
		as.integer(nrow),
		as.integer(Acol),
		as.integer(ncol),
		as.integer(A$ja),
		as.integer(A$ia),
		as.integer(B$ja),
		as.integer(B$ia),
		integer(nrow),
		nnz = integer(1),
		integer(ncol),
		PACKAGE = "SparseM")
	nnzmax <- z$nnz
	}
   z <- .Fortran("amub",
   	as.integer(nrow),
   	as.integer(ncol),
   	as.integer(1),
   	as.double(A$ra),
   	as.integer(A$ja),
   	as.integer(A$ia),
   	as.double(B$ra),
   	as.integer(B$ja),
   	as.integer(B$ia),
   	ra = double(nnzmax),
   	ja = integer(nnzmax),
   	ia = integer(nrow+1),
   	as.integer(nnzmax),
   	integer(ncol),
   	ierr = integer(1),
	PACKAGE = "SparseM")
   nnz <- z$ia[nrow+1]-1
   if(z$ierr != 0) stop("insufficient space for sparse matrix multiplication")
   z <- list(ra=z$ra[1:nnz],ja=z$ja[1:nnz],ia=z$ia,dim=c(nrow,ncol))
   class(z) <- "matrix.csr"
   }
else{
   if(is.matrix(B))stop("Can't multiply sparse times dense matrix (yet)")
   #matrix-vector multiplication: multiply a sparse csr matrix by a vector
   #	A -- csr structure returned from call to function "as.matrix.csr"
   #	B -- vector
   nrow <- A$dim[1]
   ncol <- A$dim[2]
   if(length(B) != ncol)stop("not conformable for multiplication")
   z <- .Fortran("amux",
   	as.integer(nrow),
   	as.double(B),
   	y=double(nrow),
   	as.double(A$ra),
   	as.integer(A$ja),
   	as.integer(A$ia),
	PACKAGE = "SparseM")
   z <- z$y
   }
return(z)
}
#--------------------------------------------------------------------
"chol" <- function(x, ...) UseMethod("chol")
#--------------------------------------------------------------------
"chol.default" <-get("chol", pos=NULL, mode= "function")
#--------------------------------------------------------------------
formals(chol.default) <- c(formals(chol.default), alist(... =))
#--------------------------------------------------------------------
"chol.matrix.csr" <-
function(x,cachsz=64,nsubmax,nnzlmax,tmpmax, ...){
# Interface for a sparse least squares solver via Ng-Peyton's Cholesky
# factorization
#	x -- csr structure returned from call to function "as.matrix.csr"
#       cachsz -- size of the cache memory; it's machine dependent
	nrow <- x$dim[1]
	ncol <- x$dim[2]
	if(nrow!=ncol) stop("Can't perform Cholesky Factorization for Non-square matrix\n")
	nnzdmax <- x$ia[nrow+1]-1
	nnzdsm <- nnzdmax + nrow + 1
	iwmax <- 7*nrow+3
	if(missing(nsubmax)) nsubmax <- nnzdmax
	if(missing(nnzlmax)) nnzlmax <- max(4*nnzdmax,floor(.2*nnzdmax^1.3))
	if(missing(tmpmax)) tmpmax <- 10*nrow
	level <- 8
	z <- .Fortran("chol",
		nrow = as.integer(nrow),
		nnzdmax = as.integer(nnzdmax),
		d = as.double(x$ra),
		jd = as.integer(x$ja),
		id = as.integer(x$ia),
		nnzdsm = as.integer(nnzdsm),
		dsub = double(nnzdsm),
		jdsub = integer(nnzdsm),
		nsubmax = as.integer(nsubmax),
		lindx = integer(nsubmax),	
		xlindx = integer(nrow+1),
		nsuper = integer(1),
		nnzlmax = as.integer(nnzlmax),
		lnz = double(nnzlmax),
		xlnz = integer(nrow+1),
		invp = integer(nrow),
		perm = integer(nrow),
		iwmax = as.integer(iwmax),
		iwork = integer(iwmax),
		colcnt = integer(nrow),
		snode = integer(nrow),
		xsuper = integer(nrow+1),
		split = integer(nrow),
		tmpmax = as.integer(tmpmax),
		tmpvec = double(tmpmax),
		cachsz = as.integer(cachsz),
		level = as.integer(level),
		ierr = integer(1),
		time = double(1),
		PACKAGE = "SparseM")
	if (z$ierr != 0){
        	if(z$ierr == 9) mess <- "singularity problem"
		else if(z$ierr == 4) mess <- "Increase nnzlmax"
		else if(z$ierr == 5) mess <- "Increase nsubmax"
		else if(z$ierr %in% c(8,10)) mess <- "Increase tmpmax"
        	else mess <- "insufficient space"
		if(z$ierr == 9) warning(mess) else stop(mess)
	        }
	nnzl <- z$xlnz[length(z$xlnz)]-1
        nnzlindx <- z$xlindx[nrow+1]-1
	z <- list(nrow=z$nrow,nnzlindx=nnzlindx,nsuper=z$nsuper,
                lindx=z$lindx[1:nnzlindx],xlindx=z$xlindx,nnzl=nnzl,
                lnz=z$lnz[1:nnzl],xlnz=z$xlnz,invp=z$invp,perm=z$perm,
                xsuper=z$xsuper,ierr=z$ierr,time=z$time)
	class(z) <- "matrix.csr.chol"
	return(z)
}
#--------------------------------------------------------------------
"chol.matrix.csc" <- function(x,cachsz=64, ...){
	x <- as.matrix.csr(x)
	x <- chol(x,cachsz=cachsz)
	x
}
#--------------------------------------------------------------------
"backsolve" <- function(r, x, ...) UseMethod("backsolve")
#--------------------------------------------------------------------
"backsolve.default" <-get("backsolve", pos=NULL, mode = "function")
#--------------------------------------------------------------------
formals(backsolve.default) <- c(formals(backsolve.default), alist(... =))
#--------------------------------------------------------------------
"backsolve.matrix.csr.chol" <-
function(r,x,cachsz=64, ...){
# backsolve for Ng-Peyton's Cholesky factorization
#	Solves linear system A b = x where r is chol(A)
# Input:
#	r -- structure returned by chol.matrix.csr
#	x --  rhs  may be a matrix in dense form
	m <- r$nrow
	if(!is.matrix(x)) x <- as.matrix(x) 
	if(nrow(x)!=m)stop("chol not conformable with x")
	p <- ncol(x)
	z <- .Fortran("bckslv",
		m = as.integer(m),
		nnzlindx = as.integer(r$nnzlindx),
		as.integer(r$nsuper),
		as.integer(p),
		as.integer(r$lindx),	
		as.integer(r$xlindx),
		as.integer(r$nnzl),
		as.double(r$lnz),
		as.integer(r$xlnz),
		as.integer(r$invp),
		as.integer(r$perm),
		as.integer(r$xsuper),
		double(m),
		sol = double(m*p),
		as.double(x),
		time = double(1),
		PACKAGE = "SparseM")
	z <- matrix(z$sol,nrow=nrow(x),ncol=ncol(x))
	z <- drop(z)
return(z)
}
#--------------------------------------------------------------------
"diag" <- function(x, ...) UseMethod("diag")
#--------------------------------------------------------------------
"diag.default" <- get("diag", pos= NULL, mode="function")
#--------------------------------------------------------------------
formals(diag.default) <- c(formals(diag.default), alist(... =))
#--------------------------------------------------------------------
"diag.matrix.csr" <-
function (x = 1, nrow, ...) 
{
    if (is.matrix.csr(x) && nargs() == 1) {
        if ((m <- min(dim(x))) == 0) 
            return(numeric(0))
        #y <- c(x)[1 + 0:(m - 1) * (dim(x)[1] + 1)]
	y <- rep(0,m)
        ia <- rep(1:dim(x)[1],diff(x$ia))
        y[x$ja[ia == x$ja]] <- x$ra[ia == x$ja]
	n <- sum(ia == x$ja)
        nms <- dimnames(x)
        if (is.list(nms) && !any(sapply(nms, is.null)) && all((nm <- nms[[1]][1:m]) == 
            nms[[2]][1:m])) 
            names(y) <- nm
        return(y)
    }
    if (is.array(x) && length(dim(x)) != 1) 
        stop("first argument is array, but not matrix.")
    if (missing(x)) 
        n <- nrow
    else if (length(x) == 1 && missing(nrow) ) {
        n <- as.integer(x)
        x <- 1
    }
    else n <- length(x)
    ja <- 1:n
    ra <- ja
    ra[1:n] <- x
    ia <- 1:(n+1)
    y <- list(ra = ra, ja = ja, ia = ia, dim=c(n,n))
    class(y) <- "matrix.csr"
    return(y)
}
#--------------------------------------------------------------------
"diag<-" <- function(x, value) UseMethod("diag<-")
#--------------------------------------------------------------------
"diag<-.default" <-get("diag<-", pos=NULL, mode="function")
#--------------------------------------------------------------------
"diag<-.matrix.csr" <- function(x,value)
{
	dx <- x$dim
	if (length(dx) != 2 || prod(dx) == 1) 
		stop("only matrix diagonals can be replaced")
	nrow <- x$dim[1]
	ncol <- x$dim[2]
	nnza <- x$ia[nrow+1]-1
	ia <- rep(1:nrow,diff(x$ia))
	idx <- x$ja[ia == x$ja]
	m <- nsub <- min(dx)
	ir <- jc <- 1:m
	if (length(value) != 1 && length(value) != m) 
		stop("replacement diagonal has wrong length")
        if (length(value) == 1) value <- rep(value,min(dx))
        sorted <- FALSE
        z <- .Fortran("subext",
                as.integer(nsub),
                as.integer(ir),
                as.integer(jc),
                as.double(x$ra),
                as.integer(x$ja),
                as.integer(x$ia),
                as.logical(sorted),
                value = double(nsub),
                iadd = integer(nsub),
		PACKAGE = "SparseM")
        nadd <- sum(z$iadd==0)
        nnzb <- nnza + nadd
	z <- .Fortran("subasg",
		as.integer(nrow),
		as.integer(ncol),
		as.integer(nsub),
		as.integer(nnza),
		as.integer(nnzb),
		as.integer(ir),
		as.integer(jc),
		as.double(x$ra),
		as.integer(x$ja),
		as.integer(x$ia),
		ra = double(nnzb),
		ja = integer(nnzb),
		ia = integer(nrow+1),
		value = as.double(value),
		logical(ncol),
		ierr = integer(1),
		PACKAGE = "SparseM")
        if(z$ierr != 0) stop("not enough space")
	x <- list(ra = z$ra, ja = z$ja, ia = z$ia, dim=x$dim)
	class(x) <- "matrix.csr"
	x
}
#--------------------------------------------------------------------
solve.matrix.csr <-
function (a, b, ...) 
{
    if(!is.matrix.csr(a))stop("a not in csr format")
    nr <- nrow(a)
    nc <- ncol(a)
    if (nc != nr) 
         stop("only square systems can be solved")
    a <- chol(a,...)
    if (missing(b)) {
        b <- diag(1, nc)
    }
    else
        if(!is.matrix.csr(b))b <- as.matrix(b) 
    backsolve(a,b)
}
#--------------------------------------------------------------------
"slm" <-
function (formula,  data, weights, na.action, method = "csr", 
    contrasts = NULL, ...) 
{
    call <- match.call()
    m <- match.call(expand.dots = FALSE)
    m$method <- m$model <- m$x <- m$y <- m$contrasts <-  m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.frame(sys.parent()))
    if (method == "model.frame") 
        return(m)
    Terms <- attr(m, "terms")
    weights <- model.extract(m, weights)
    Y <- model.extract(m, response)
    X <- as.matrix.csr(model.matrix(Terms, m, contrasts))
    fit <- {
        if (length(weights)) 
            slm.wfit(X, Y,  weights, method, ...)
        else slm.fit(X, Y,  method, ...)
    }
    fit$terms <- Terms
    fit$call <- call
    attr(fit, "na.message") <- attr(m, "na.message")
    class(fit) <-  "slm"
    fit
}
#--------------------------------------------------------------------
"slm.fit" <-
function (x, y,  method = "csr", ...) 
{
    fit <-  slm.fit.csr(x, y,  ...)
    fit$contrasts <- attr(x, "contrasts")
    fit
}
#--------------------------------------------------------------------
"slm.wfit" <-
function (x, y, weights,  ...) 
{
    if (!is.matrix.csr(x)) 
        stop("model matrix must be in sparse csr mode")
    if (!is.numeric(y)) 
        stop("response must be numeric")
    if (any(weights < 0)) 
        stop("negative weights not allowed")
    contr <- attr(x, "contrasts")
    w <- sqrt(weights)
    x <- t(t(x) %*% diag.matrix.csr(w))
    y <- y * w
    fit <- slm.fit.csr(x, y,  ...)
    fit$contrasts <- attr(x, "contrasts")
    fit
}
#--------------------------------------------------------------------
"slm.fit.csr" <-
function (x, y, ...) 
{
    n <- length(y)
    p <- x$dim[2]
    if (n != x$dim[1]) 
        stop("x and y don't match n")
    chol <- chol(t(x)%*%x)
    xy <- t(x) %*% y
    coefficients <- backsolve(chol,xy)
    #if (z$info != 0) 
    #    stop(paste("Error info = ", z$info, "in stepy: singular design"))
    #names(coefficients) <- dimnames(x)[[2]]
    fitted <-  x %*% coefficients
    residuals <- y - fitted
    return(coefficients, chol, residuals, fitted)
}
#--------------------------------------------------------------------
"coef.slm" <-
function (object, ...) 
object$coefficients
#--------------------------------------------------------------------
"fitted.slm" <-
function (object, ...) 
object$fitted
#--------------------------------------------------------------------
"residuals.slm" <- function (object, ...) 
r <- object$residuals
#--------------------------------------------------------------------
"summary.slm" <-
function (object, correlation = FALSE, ...)
{
    Chol <- object$chol
    if (is.null(object$terms) || is.null(Chol))
        stop("invalid 'lm' object:  no terms or chol component")
    n <- length(object$residuals)
    p <- length(object$coefficients)
    rdf <- n - p
    r <- resid(object)
    f <- fitted(object)
    w <- weights(object)
    if (is.null(w)) {
        mss <- if (attr(object$terms, "intercept"))
            sum((f - mean(f))^2)
        else sum(f^2)
        rss <- sum(r^2)
    }
    else {
        mss <- if (attr(object$terms, "intercept")) {
            m <- sum(w * f/sum(w))
            sum(w * (f - m)^2)
        }
        else sum(w * f^2)
        rss <- sum(w * r^2)
        r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    R <- backsolve(Chol,diag(p))
    se <- sqrt(diag(R) * resvar)
    est <- coefficients(object)
    tval <- est/se
    ans <- object[c("call", "terms")]
    ans$residuals <- r
    ans$coefficients <- cbind(est, se, tval, 2 * (1 - pt(abs(tval), rdf)))
    dimnames(ans$coefficients) <- list(names(object$coefficients),
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans$sigma <- sqrt(resvar)
    ans$df <- c(p, rdf, n)
    if (p != attr(object$terms, "intercept")) {
        df.int <- if (attr(object$terms, "intercept"))
            1
        else 0
        ans$r.squared <- mss/(mss + rss)
        ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n -
            df.int)/rdf)
        ans$fstatistic <- c(value = (mss/(p - df.int))/resvar,
            numdf = p - df.int, dendf = rdf)
    }
    ans$cov.unscaled <- R
    dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1,
        1)]
    if (correlation) {
        ans$correlation <- (R * resvar)/outer(se, se)
        dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
    }
    class(ans) <- "summary.slm"
    ans
}
#--------------------------------------------------------------------
"print.summary.slm" <-
function (x, digits = max(3, getOption("digits") - 3), 
	symbolic.cor = p > 4, signif.stars = getOption("show.signif.stars"), 
	...)
{
    cat("\nCall:\n")
    cat(paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    resid <- x$residuals
    df <- x$df
    rdf <- df[2]
    cat(if (!is.null(x$w) && diff(range(x$w)))
        "Weighted ", "Residuals:\n", sep = "")
    if (rdf > 5) {
        nam <- c("Min", "1Q", "Median", "3Q", "Max")
        rq <- if (length(dim(resid)) == 2)
            structure(apply(t(resid), 1, quantile), dimnames = list(nam,
                dimnames(resid)[[2]]))
        else structure(quantile(resid), names = nam)
        print(rq, digits = digits, ...)
    }
    else if (rdf > 0) {
        print(resid, digits = digits, ...)
    }
    else {
        cat("ALL", df[1], "residuals are 0: no residual degrees of freedom!\n")
    }
#    if (nsingular <- df[3] - df[1])
#        cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n",
#            sep = "")
    #else cat("\nCoefficients:\n")
    cat("\nCoefficients:\n")
    print.coefmat(x$coef, digits = digits, signif.stars = signif.stars,
        ...)
    cat("\nResidual standard error:", format(signif(x$sigma,
        digits)), "on", rdf, "degrees of freedom\n")
    if (!is.null(x$fstatistic)) {
        cat("Multiple R-Squared:", formatC(x$r.squared, digits = digits))
        cat(",\tAdjusted R-squared:", formatC(x$adj.r.squared,
            digits = digits), "\nF-statistic:", formatC(x$fstatistic[1],
            digits = digits), "on", x$fstatistic[2], "and", x$fstatistic[3],
            "DF,\tp-value:", formatC(1 - pf(x$fstatistic[1],
                x$fstatistic[2], x$fstatistic[3]), dig = digits),
            "\n")
    }
    correl <- x$correlation
    if (!is.null(correl)) {
        p <- NCOL(correl)
        if (p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            if (symbolic.cor)
                print(symnum(correl)[-1, -p])
            else {
                correl[!lower.tri(correl)] <- NA
                print(correl[-1, -p, drop = FALSE], digits = digits,
                  na = "")
            }
        }
    }
    cat("\n")
    invisible(x)
}
#--------------------------------------------------------------------
