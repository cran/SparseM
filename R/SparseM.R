#--------------------------------------------------------------------
#".First.lib" <- function(lib, pkg) {
#   require(methods)
#   library.dynam("SparseM", pkg, lib)
#   print("SparseM library loaded")
#}
.onLoad <- function(lib, pkg) {
	require(methods)
	print("SparseM library loaded")
}
# /*RSB*/ .onLoad needed instead of .First.lib, no library.dynam()
#--------------------------------------------------------------------
"is.matrix.csr" <- function(x, ...) is(x,"matrix.csr")
#--------------------------------------------------------------------
"is.matrix.csc" <- function(x, ...) is(x,"matrix.csc")
#--------------------------------------------------------------------
"is.matrix.ssr" <- function(x, ...) is(x,"matrix.ssr")
#--------------------------------------------------------------------
"is.matrix.ssc" <- function(x, ...) is(x,"matrix.ssc")
#--------------------------------------------------------------------
"is.matrix.coo" <- function(x, ...) is(x,"matrix.coo")
#--------------------------------------------------------------------
"as.matrix.csr" <-
function(x, nrow = 1, ncol = 1, eps = .Machine$double.eps){
	 if(is.matrix.csr(x)) {x; return(x)}
         if (!is.matrix(x)) {
		if (missing(nrow))
                         nrow <- ceiling(length(x)/ncol)
		else if (missing(ncol))
                         ncol <- ceiling(length(x)/nrow)
		if (length(x) == nrow * ncol)
                         x <- matrix(x, nrow, ncol)
		else{
			if(length(x)==1 && abs(x)<eps) {
				dimx <- c(nrow,ncol)
	        		z<-new("matrix.csr",ra=0,ja=as.integer(1),
					ia=as.integer(c(1:1,rep(2,dimx[1]))), 
					dimension=as.integer(dimx))
				return(z)
				}
			else if((nrow*ncol)%%length(x)!=0){
				R <- ceiling(nrow*ncol/length(x))
                                x <- matrix(rep(x,R), nrow, ncol)
				warning("ncol*nrow indivisable by length(x)")
				}
			else {
				R <- ceiling(nrow*ncol/length(x))
				x <- matrix(rep(x,R), nrow, ncol)
				}
			}
		}
	dimx <- dim(x)
	nnz <- sum(abs(x)>eps)
	if(nnz==0){
	        z<-new("matrix.csr",ra=0,ja=as.integer(1),
			ia=as.integer(c(1:1,rep(2,dimx[1]))), dimension=dimx)
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
	z <- new("matrix.csr",ra = z$ra[1:nnz], ja = z$ja[1:nnz], ia = z$ia, dimension = dimx)
	return(z)
}
#--------------------------------------------------------------------
"as.matrix.csc" <- function(x, nrow = 1, ncol = 1, eps = .Machine$double.eps)
{
    if (is.matrix.csc(x)) x
    else as.matrix.csc(as.matrix.csr(x))
#    else as.matrix.csc(as.matrix.csr(x, nrow = 1, ncol = 1, eps = .Machine$double.eps))
}
#--------------------------------------------------------------------
"as.matrix.ssr" <- function(x, nrow = 1, ncol = 1, eps = .Machine$double.eps)
{
	if (is.matrix.ssr(x)) x
	else as.matrix.ssr(as.matrix.csr(x))
#	else as.matrix.ssr(as.matrix.csr(x, nrow = 1, ncol = 1, eps = .Machine$double.eps))
}
#--------------------------------------------------------------------
"as.matrix.ssc" <- function(x, nrow = 1, ncol = 1, eps = .Machine$double.eps)
{
	if (is.matrix.ssc(x)) x
	else as.matrix.ssc(as.matrix.csc(x))
#	else as.matrix.ssc(as.matrix.csc(x, nrow = 1, ncol = 1, eps = .Machine$double.eps))
}
#--------------------------------------------------------------------
"as.matrix.coo" <- function(x, nrow = 1, ncol = 1, eps = .Machine$double.eps)
{
	if (is.matrix.coo(x) && missing(nrow) && missing(ncol)) x
	else as.matrix.coo(as.matrix.csr(x))
}
#--------------------------------------------------------------------
#"ncol.matrix.csr" <-
#function(x){dim(x)[2]}
#--------------------------------------------------------------------
#"nrow.matrix.csr" <-
#function(x){dim(x)[1]}
#--------------------------------------------------------------------
".ssr.csr" <- function(x){
	nrow <- x@dimension[1]
	nnza <- x@ia[nrow+1]-1
	nnzao <- 2*nnza #can be set smaller
	z <- .Fortran("ssrcsr",
		job = as.integer(0),
		value2 = as.integer(1),
		nrow = as.integer(nrow),
		a = as.double(x@ra),
		ja = as.integer(x@ja),
		ia = as.integer(x@ia),
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
	z <- new("matrix.csr",ra=z$ao[1:nnz],ja=z$jao[1:nnz],ia=z$iao,dimension=x@dimension)
	return(z)
}
#--------------------------------------------------------------------
".csr.ssr" <- function(x){
	nrow <- x@dimension[1]
	nnza <- ceiling((x@ia[nrow+1]-1)/2)+nrow
	if(nrow!=x@dimension[2]) 
                stop("Cannot convert an asymmetric matrix into `matrix.ssr' class")

        if(sum(abs((t(as.matrix.csr(x))-as.matrix.csr(x))@ra))!=0)
                stop("Cannot convert an asymmetric matrix into `matrix.ssr' class")
	z <- .Fortran("csrssr",
		as.integer(nrow),
		as.double(x@ra),
		as.integer(x@ja),
		as.integer(x@ia),
		as.integer(nnza),
		ao = as.double(x@ra),
		jao = as.integer(x@ja),
		iao = as.integer(x@ia),
		ierr = integer(1),
		PACKAGE = "SparseM")
	if(z$ierr != 0) stop("Not enough space. This is usually caused by trying to convert an asymmetric matrix into ssr format")
	nnza <- z$iao[nrow+1]-1
	z <- new("matrix.ssr",ra=z$ao[1:nnza],ja=z$jao[1:nnza],ia=z$iao,dimension=x@dimension)
	z
}
#--------------------------------------------------------------------
".csc.ssc" <- function(x){
	nrow <- x@dimension[2]
	nnza <- ceiling((x@ia[nrow+1]-1)/2)+nrow
	if(nrow!=x@dimension[1])
                stop("Cannot convert an asymmetric matrix into `matrix.ssc' class")
        if(sum(abs((t(as.matrix.csr(x))-as.matrix.csr(x))@ra))!=0)
                stop("Cannot convert an asymmetric matrix into `matrix.ssc' class")
	z <- .Fortran("cscssc",
		as.integer(nrow),
		as.double(x@ra),
		as.integer(x@ja),
		as.integer(x@ia),
		as.integer(nnza),
		ao = as.double(x@ra),
		jao = as.integer(x@ja),
		iao = as.integer(x@ia),
		ierr = integer(1),
		PACKAGE = "SparseM")
	if(z$ierr != 0) stop("Not enough space. This is usually caused by trying to convert an asymmetric matrix into ssc format")
	nnza <- z$iao[nrow+1]-1
	z <- new("matrix.ssc",ra=z$ao[1:nnza],ja=z$jao[1:nnza],ia=z$iao,dimension=x@dimension)
	return(z)
}
#--------------------------------------------------------------------
".csr.coo" <- function (x) {
    nrow <- x@dimension[1]
    ncol <- x@dimension[2]
    nnza <- length(x@ra)
    z <- .Fortran("csrcoo", as.integer(nrow), as.integer(1),
        as.integer(nnza), as.double(x@ra), as.integer(x@ja),
        as.integer(x@ia), nnz = integer(1), ao = as.double(x@ra),
        ir = integer(nnza), jc = as.integer(x@ja), ierr = integer(1),
        PACKAGE = "SparseM")
    if (z$ierr != 0)
        stop("Not enough space.")
    z <- new("matrix.coo", ra = x@ra, ja = x@ja, ia = z$ir, dimension = x@dimension)
    return(z)
}
#--------------------------------------------------------------------
"rbind.matrix.csr" <- function(...) {
# Very preliminary function to rbind matrix.csr objects no name handling
    allargs <- list(...)
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
        ra <- c(ra,xi@ra)
        ja <- as.integer(c(ja,xi@ja))
        ia <- as.integer(c(ia[-(Nrow+1)],nia + xi@ia))
        nia <- ia[length(ia)]-1
        Nrow <- Nrow + length(xi@ia)-1
        }
z <- new("matrix.csr", ra=ra, ja=ja, ia = ia, dimension = as.integer(c(Nrow,Ncol)))
return(z)
}
#--------------------------------------------------------------------
"cbind.matrix.csr" <- function(...)
{
# Very preliminary function to cbind matrix.csr objects no name handling
    allargs <- list(...)
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
        ra <- c(ra,xi@ra)
        ja <- as.integer(c(ja,xi@ja))
        ia <- as.integer(c(ia[-(Nrow+1)],nia + xi@ia))
        nia <- ia[length(ia)]-1
        Nrow <- Nrow + length(xi@ia)-1
        }
z <- new("matrix.csr",ra=ra, ja=ja, ia = ia, dimension = as.integer(c(Nrow,Ncol)))
z <- t(z)
return(z)
}
#--------------------------------------------------------------------
"read.matrix.hb" <-
function (filename) 
{
	if(TRUE){
		cat("This function is undergoing therapy.\n")
		return()
		}
	hb1.o <- .C("read_HB1", 
		infile = as.character(filename),
		M = integer(1),
		N = integer(1),
		nnz = integer(1),
		Nrhs = integer(1),
		mxtype = character(1),
		Rhstype = character(1),
		errflg = integer(1),
		PACKAGE = "SparseM"
		)
	if(hb1.o$errflg == -1) stop(paste("Can't find",filename,sep = " "))
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
		if(format == "csc")
			rd.o <- new("matrix.csc.hb", ra = hb2.o$val, ja = hb2.o$rowind, ia = hb2.o$colptr, 
			rhs.ra = hb3.o$rhs, guess = switch(Gflag, Yes = hb4.o$rhs, No = NULL), 
			xexact = switch(Xflag, Yes = hb5.o$rhs, No = NULL), dimension = c(hb1.o$M, hb1.o$N),
			rhs.dim = c(hb1.o$M, hb1.o$Nrhs), rhs.mode=rhs.mode)
		else
			rd.o <- new("matrix.ssc.hb", ra = hb2.o$val, ja = hb2.o$rowind, ia = hb2.o$colptr, 
			rhs.ra = hb3.o$rhs, guess = switch(Gflag, Yes = hb4.o$rhs, No = NULL), 
			xexact = switch(Xflag, Yes = hb5.o$rhs, No = NULL), dimension = c(hb1.o$M, hb1.o$N),
			rhs.dim = c(hb1.o$M, hb1.o$Nrhs), rhs.mode=rhs.mode)
		}
	else
        	if(format=="csc")
			rd.o <- new("matrix.csc.hb",ra = hb2.o$val, ja = hb2.o$rowind, ia = hb2.o$colptr,
			dimension=c(hb1.o$M,hb1.o$N),rhs.mode=rhs.mode)
        	else
			rd.o <- new("matrix.ssc.hb",ra = hb2.o$val, ja = hb2.o$rowind, ia = hb2.o$colptr,
			dimension=c(hb1.o$M,hb1.o$N),rhs.mode=rhs.mode)

   return(rd.o)
}
#--------------------------------------------------------------------
"write.matrix.hb" <-
function (filename = "hb.out", X, title, key, mxtype, rhs = NULL, 
    guess = FALSE, xsol = FALSE, ptrfmt, indfmt,
    valfmt = "(1P,5D16.9)", rhsfmt = "(1P,5D16.9)") 
{
     if(TRUE){
           cat("This function is undergoing therapy.\n")
           return()
           }
    if (!substr(mxtype, 1, 1) %in% c("r", "R")) 
        stop("The first character of `mxtype' can only be 'R'")
    if (!substr(mxtype, 2, 2) %in% c("s", "S", "u", "U", "r", 
        "R")) 
        stop("The second character of `mxtype' can only be `S',`U' or 'R'")
    if (!substr(mxtype, 3, 3) %in% c("a", "A")) 
        stop("The third character of `mxtype' can only be `A'")
    if (substr(mxtype, 2, 2) %in% c("s", "S") && !is.matrix.ssc(X)) 
        stop("Matrix X has to be in in ssc format")
    if (substr(mxtype, 2, 2) %in% c("u", "U", "r", "R") && !is.matrix.csc(X)) 
        stop("Matrix X has to be in in csc format")
    nch <- nchar(as.integer(max(X@ia)))+1
    if(missing(ptrfmt)) ptrfmt <- paste("(",floor(80/nch),"I",nch,")",sep="")
    nch.msg <- strsplit(ptrfmt,"I")[[1]][2]
    if(substr(nch.msg,1,nchar(nch.msg)-1) < nch-1) stop("Bad format specification")
    nch <- nchar(as.integer(max(X@ja)))+1
    if(missing(indfmt)) indfmt <- paste("(",floor(80/nch),"I",nch,")",sep="")
    nch.msg <- strsplit(indfmt,"I")[[1]][2]
    if(substr(nch.msg,1,nchar(nch.msg)-1) < nch-1) stop("Bad format specification")
    M <- X@dimension[1]
    N <- X@dimension[2]
    nnz <- length(X@ra)
    nrhs <- 0
    guesol <- ""
    Rhs <- Guess <- Exact <- rep(0, M)
    if (!missing(rhs)) {
        guesol <- paste(guesol, "F", sep = "")
        idiv <- 1
        Rhs <- rhs[1:(M * idiv)]
        if (guess) {
            idiv <- idiv + 1
            guesol <- paste(guesol, "G", sep = "")
            Guess <- rhs[(M * (idiv - 1) + 1):(M * idiv)]
        }
        if (xsol) {
            idiv <- idiv + 1
            guesol <- paste(guesol, "X", sep = "")
            Exact <- rhs[(M * (idiv - 1) + 1):(M * idiv)]
        }
        nrhs <- length(rhs)/M/idiv
        if (length(rhs)%%(M * idiv) != 0) 
            stop("The length of `rhs' is not a multiple of the number of equations")
    }
    .C("write_HB1", as.character(filename), as.integer(M), as.integer(N), 
        as.integer(nnz), as.integer(X@ia), as.integer(X@ja), 
        as.double(X@ra), as.integer(nrhs), as.double(Rhs), as.double(Guess), 
        as.double(Exact), as.character(title), as.character(key), 
        as.character(mxtype), as.character(guesol), as.character(ptrfmt), 
        as.character(indfmt), as.character(valfmt), as.character(rhsfmt), 
        PACKAGE = "SparseM")
    invisible(X)
}
#--------------------------------------------------------------------
"Ops.matrix.csr" <- function(e1,e2){
	if(missing(e2)){
		e1.op <- switch(.Generic,
			"+" = e1,
			"-" = new("matrix.csr",ra=-e1@ra,ja=e1@ja,ia=e1@ia,dimension=e1@dimension),
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
		"&" = {z <- .matrix.csr.elemul(e1,e2);z@ra <- rep(1,length(z@ja));z},
		"|" = {z <- .matrix.csr.addsub(e1,e2,1);z@ra <- rep(1,length(z@ja));z},
		stop(paste("Binary operator \"",.Generic,"\""," is undefined for class \"matrix.csr\"",sep=""))
		)
		}
	return(e1.op.e2)
}
#--------------------------------------------------------------------
"Ops.matrix.diag.csr" <- function(e1,e2){
        if(missing(e2)){
                e1.op <- switch(.Generic,
                        "+" = e1,
                        "-" = new("matrix.csr",ra=-e1@ra,ja=e1@ja,ia=e1@ia,dimension=e1@dimension),
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
                "&" = {z <- .matrix.csr.elemul(e1,e2);z@ra <- rep(1,length(z@ja));z},
                "|" = {z <- .matrix.csr.addsub(e1,e2,1);z@ra <- rep(1,length(z@ja));z},
                stop(paste("Binary operator \"",.Generic,"\""," is undefined for class \"matrix.csr\"",sep=""))
                )
                }
        return(e1.op.e2)
}

#--------------------------------------------------------------------
".matrix.csr.compl" <- function(e1){
	nrow <- e1@dimension[1]
	ncol <- e1@dimension[2]
	nnz <- e1@ia[nrow+1]-1
	nz <- nrow*ncol - nnz
	if(nz == 0){ # full matrix
		z <- as.matrix.csr(0,nrow,ncol)
		return(z)
		}
	if(nnz == 1 && e1@ra == 0){ # zero matrix
		z <- as.matrix.csr(1,nrow,ncol)
		return(z)
		}
	if(length(e1@ra) == 1 && is.na(e1@ra)){ #trap zero matrix
		z <- list(ra=rep(1,nz),ja=rep(1:ncol,nrow),ia=seq(1,nz+1,by=ncol),dim=e1@dimension)
		}
	else{
		z <- .Fortran("nzero",
			as.double(e1@ra),
			as.integer(e1@ja),
			as.integer(e1@ia),
			as.integer(nrow),
			as.integer(ncol),
			as.integer(nnz),
			as.integer(nz),
			ra = double(nz),
			ja = integer(nz),
			ia = integer(nrow+1),
			logical(ncol),
			PACKAGE = "SparseM")
		z <- new("matrix.csr",ra=z$ra,ja=z$ja,ia=z$ia,dimension=e1@dimension)
		}
	z
}
#--------------------------------------------------------------------
".matrix.csr.addsub" <-
function(A,B,s){
#matrix addition/subtraction of two sparse csr matrices 
if(is.matrix(A)) A <- as.matrix.csr(A)
if(is.matrix(B)) B <- as.matrix.csr(B)
nrow <- A@dimension[1]
ncol <- A@dimension[2]
Bcol <- B@dimension[2]
Brow <- B@dimension[1]
if(ncol != Bcol | nrow != Brow)stop("matrices not conformable for addition")
nnza <- A@ia[nrow+1]-1
nnzb <- B@ia[nrow+1]-1
nnzmax <- length(union(A@ja+A@dimension[2]*(rep(1:A@dimension[1],diff(A@ia))-1),
	B@ja+B@dimension[2]*(rep(1:B@dimension[1],diff(B@ia))-1)))+1
z <- .Fortran("aplsb",
	as.integer(nrow),
	as.integer(ncol),
	as.integer(1),
	as.double(A@ra),
	as.integer(A@ja),
	as.integer(A@ia),
	as.double(s),
	as.double(B@ra),
	as.integer(B@ja),
	as.integer(B@ia),
	ra = double(nnzmax),
	ja = integer(nnzmax),
	ia = integer(nrow+1),
	as.integer(nnzmax),
	integer(ncol),
	ierr = integer(1),
	PACKAGE = "SparseM")
if(z$ierr != 0) stop("insufficient space for sparse matrix addition")
nnz <- z$ia[nrow+1]-1
z <- new("matrix.csr",ra=z$ra[1:nnz],ja=z$ja[1:nnz],ia=z$ia,dimension=c(nrow,ncol))
return(z)
}
#--------------------------------------------------------------------
".matrix.csr.elemul" <- function(A,B){
if(is.vector(A)) {
        if(length(A) == 1){
                if(A==0) return(as.matrix.csr(0,nrow(B),ncol(B)))
                else{B@ra <- A*B@ra;return(B)}
                }
        else if(length(A) == nrow(B))
		return(as(A,"matrix.diag.csr") %*% B)
        else if(length(A) == ncol(B))
		return(B %*% as(A,"matrix.diag.csr"))
        else
                stop("A and B not conformable for element-by-element multiplication")
        }
else if(is.vector(B)) {
        if(length(B) == 1){
                if(B==0) return(as.matrix.csr(0,nrow(A),ncol(A)))
                else{A@ra <- B*A@ra;return(A)}
                }
        else if(length(B) == nrow(A))
		return(as(B,"matrix.diag.csr") %*% A)
        else if(length(B) == ncol(A))
		return(A %*% as(B,"matrix.diag.csr"))
        else
                stop("A and B not conformable for element-by-element multiplication")
        }
if(is.matrix(A))
        A <- as.matrix.csr(A)
else if(is.matrix(B))
        B <- as.matrix.csr(B)
if(!(is.matrix.csr(A) && is.matrix.csr(B)))
        stop("Arguments must be of class:  vector, matrix or matrix.csr")
else
	Arow <- nrow(A)
        Acol <- ncol(A)
        Brow <- nrow(B)
        Bcol <- ncol(B)
        if(Acol != Bcol | Arow != Brow)
                stop("A and B not conformable for element-by-element multiplication")
        nnza <- A@ia[Arow+1]-1
        nnzb <- B@ia[Arow+1]-1
        nnzmax <- length(intersect(A@ja+A@dimension[2]*(rep(1:A@dimension[1],diff(A@ia))-1),
                B@ja+B@dimension[2]*(rep(1:B@dimension[1],diff(B@ia))-1)))+1
        z <- .Fortran("aemub",
                as.integer(Arow),
                as.integer(Acol),
                as.double(A@ra),
                as.integer(A@ja),
                as.integer(A@ia),
                as.double(B@ra),
                as.integer(B@ja),
                as.integer(B@ia),
                ra = double(nnzmax),
                ja = integer(nnzmax),
                ia = integer(Arow+1),
                integer(Acol),
                double(Acol),
                as.integer(nnzmax),
                ierr = integer(1),
                PACKAGE = "SparseM")
	if(z$ierr != 0)
                stop("insufficient space for element-wise sparse matrix multiplication")
        nnz <- z$ia[Arow+1]-1
        if(identical(z$ra,0)){#trap zero matrix
                z$ja <- as.integer(1)
                z$ia <- as.integer(c(1,rep(2,nrow)))
                }
z <- new("matrix.csr",ra=z$ra[1:nnz],ja=z$ja[1:nnz],ia=z$ia,dimension=c(Arow,Acol))
return(z)
}
#--------------------------------------------------------------------
".matrix.csr.elediv" <- function(A,B){
# Element-wise matrix division of two sparse csr matrices 
if(is.numeric(A) && length(A) == 1)
        z <- new("matrix.csr",ra=A/B@ra,ja=B@ja,ia=B@ia,dimension=B@dimension)
else if(is.numeric(B) && length(B) == 1)
        z <- new("matrix.csr",ra=A@ra/B,ja=A@ja,ia=A@ia,dimension=A@dimension)
else if(is.matrix.csr(A) || is.matrix.csr(B) || is.matrix(A) || is.matrix(B)){
        if(is.matrix(A)) A <- as.matrix.csr(A)
        if(is.matrix(B)) B <- as.matrix.csr(B)
        nrow <- A@dimension[1]
        ncol <- A@dimension[2]
        Bcol <- B@dimension[2]
        Brow <- B@dimension[1]
        if(ncol != Bcol | nrow != Brow)stop("matrices not conformable for element-by-element division")
        nnza <- A@ia[nrow+1]-1
        nnzb <- B@ia[nrow+1]-1
	nnzmax <- length(union(A@ja+A@dimension[2]*(rep(1:A@dimension[1],diff(A@ia))-1),
                B@ja+B@dimension[2]*(rep(1:B@dimension[1],diff(B@ia))-1)))+1
        z <- .Fortran("aedib",
                as.integer(nrow),
                as.integer(ncol),
		as.integer(1),
                as.double(A@ra),
                as.integer(A@ja),
                as.integer(A@ia),
                as.double(B@ra),
                as.integer(B@ja),
                as.integer(B@ia),
                ra = double(nnzmax),
                ja = integer(nnzmax),
                ia = integer(nrow+1),
                as.integer(nnzmax),
		integer(ncol),
		double(ncol),
                ierr = integer(1),
                PACKAGE = "SparseM")
        if(z$ierr != 0) stop("insufficient space for element-wise sparse matrix division")
        nnz <- z$ia[nrow+1]-1
        z1 <- vector("numeric",nrow*ncol)
        idx1 <- z$ja[1:nnz]+ncol*(rep(1:nrow,diff(z$ia))-1)
        idx2 <- union(A@ja+A@dimension[2]*(rep(1:A@dimension[1],diff(A@ia))-1),
                B@ja+B@dimension[2]*(rep(1:B@dimension[1],diff(B@ia))-1))
        idx3 <- setdiff(1:(nrow*ncol),idx2)
        z1[idx1] <- z$ra[1:nnz]
        z1[idx3] <- NaN
        z <- new("matrix.csr",ra=z1,ja=as.integer(rep(1:ncol,nrow)),
		ia=as.integer(seq(1,nrow*ncol+1,by=ncol)),dimension=as.integer(c(nrow,ncol)))
        }
else stop("Arguments have to be class \"matrix.csr\" or numeric")
return(z)
}
#--------------------------------------------------------------------
".matrix.csr.expo" <- function(A,B){
# Performs element-wise exponentiation on sparse matrices
if(is.numeric(A) && length(A) == 1)
        z <- new("matrix.csr",ra=A/B@ra,ja=B@ja,ia=B@ia,dimension=B@dimension)
else if(is.numeric(B) && length(B) == 1)
        z <- new("matrix.csr",ra=B/A@ra,ja=A@ja,ia=A@ia,dimension=A@dimension)
else if(is.matrix.csr(A) || is.matrix.csr(B) || is.matrix(A) || is.matrix(B)){
        if(is.matrix(A)) A <- as.matrix.csr(A)
        if(is.matrix(B)) B <- as.matrix.csr(B)
        nrow <- A@dimension[1]
        ncol <- A@dimension[2]
        Bcol <- B@dimension[2]
        Brow <- B@dimension[1]
        if(ncol != Bcol | nrow != Brow)stop("matrices not conformable for element-by-element division")
        nnza <- A@ia[nrow+1]-1
        nnzb <- B@ia[nrow+1]-1
	nnzmax <- length(union(A@ja+A@dimension[2]*(rep(1:A@dimension[1],diff(A@ia))-1),
                B@ja+B@dimension[2]*(rep(1:B@dimension[1],diff(B@ia))-1)))+1
        z <- .Fortran("aeexpb",
                as.integer(nrow),
                as.integer(ncol),
		as.integer(1),
                as.double(A@ra),
                as.integer(A@ja),
                as.integer(A@ia),
                as.double(B@ra),
                as.integer(B@ja),
                as.integer(B@ia),
                ra = double(nnzmax),
                ja = integer(nnzmax),
                ia = integer(nrow+1),
                as.integer(nnzmax),
		integer(ncol),
		double(ncol),
                ierr = integer(1),
                PACKAGE = "SparseM")
        if(z$ierr != 0) stop("insufficient space for element-wise sparse matrix exponentiation")
        nnz <- z$ia[nrow+1]-1
        z1 <- vector("numeric",nrow*ncol)
        idx1 <- z$ja[1:nnz]+ncol*(rep(1:nrow,diff(z$ia))-1)
	idxA <- A@ja+A@dimension[2]*(rep(1:A@dimension[1],diff(A@ia))-1)
	idxB <- B@ja+B@dimension[2]*(rep(1:B@dimension[1],diff(B@ia))-1)
        idx2 <- union(idxA,idxB)
        idx3 <- setdiff(1:(nrow*ncol),idx2)
	idx4 <- setdiff(idxB,idxA)
	idxInf <- (idxB[B@ra<0])[!is.element(idxB[B@ra<0],idxA)]
        z1[idx1] <- z$ra[1:nnz]
        z1[idx3] <- 1
	z1[idxInf] <- Inf  #this is needed because Fortran returns NA instead of Inf for 0^-1
        z <- new("matrix.csr",ra=z1,ja=as.integer(rep(1:ncol,nrow)),
		ia=as.integer(seq(1,nrow*ncol+1,by=ncol)),dimension=as.integer(c(nrow,ncol)))
        }
else stop("Arguments have to be class \"matrix.csr\" or numeric")
return(z)
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
        AB@ra[AB@ra==uniq[1]] <- NaN
        AB@ra[AB@ra==uniq[2]] <- Inf
        AB@ra[AB@ra==uniq[3]] <- -Inf
	as(AB,"matrix.csr")
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
        AB@ra[AB@ra==uniq[1]] <- NaN
        AB@ra[AB@ra==uniq[2]] <- Inf
        AB@ra[AB@ra==uniq[3]] <- -Inf
	as(AB,"matrix.csr")
	AB
}
#--------------------------------------------------------------------
".matrix.csr.relation" <- function(e1,e2,rel){
	if(is.numeric(e2) && length(e2) == 1){
		z <- .csr.relation(e1,e2,rel)
		}
	else if(is.numeric(e1) && length(e1) == 1){
		z <- .csr.relation(-e2,-e1,rel)
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
	if(!(length(z@ra)==1 && z@ra==0 && z@ja==1)) #trap returned matrix with all zeros
		z@ra <- rep(1,length(z@ra))
	z
}
#--------------------------------------------------------------------
".csr.relation" <- function(A,drptol,rel){
	nrow <- A@dimension[1]
	nnza <- A@ia[nrow+1]-1
	flag <- FALSE
	if((rel=="gt" || rel=="le") && drptol >=0){
		relidx <- 1
		flag <- TRUE
		}
	if((rel=="lt" || rel=="ge") && drptol <=0){
		relidx <- 1
		drptol <- -drptol
		A@ra <- -A@ra
		flag <- TRUE
		}
	if(rel=="eq" && drptol !=0){
		relidx <- 3
		flag <- TRUE
		}
	if((rel=="ne" || rel=="eq") && drptol ==0){
		relidx <- 4
		flag <- TRUE
		}
	if(flag){
		z <- .Fortran("filter1",
			as.integer(nrow),
			as.integer(relidx),
			as.double(drptol),
			as.double(A@ra),
			as.integer(A@ja),
			as.integer(A@ia),
			ra = double(nnza),
			ja = integer(nnza),
			ia = integer(nrow+1),
			as.integer(nnza),
			ierr = integer(1),
		PACKAGE = "SparseM")
		if(z$ierr !=0) stop("Not enough space")
		nnza <- z$ia[nrow+1]-1
		if(nnza==0){ # trap returned matrix with all zeros
			z$ra <- 0
			z$ja <- as.integer(1)
			z$ia <- as.integer(c(1,rep(2,A@dimension[1])))
			}
		if(rel == "lt")
			z <- new("matrix.csr",ra=z$ra[1:nnza],ja=z$ja[1:nnza],ia=z$ia,dimension=A@dimension)
		else if(rel == "gt")
			z <- new("matrix.csr",ra=-z$ra[1:nnza],ja=z$ja[1:nnza],ia=z$ia,dimension=A@dimension)
		else if(rel == "le"){
			z <- new("matrix.csr",ra=-z$ra[1:nnza],ja=z$ja[1:nnza],ia=z$ia,dimension=A@dimension)
			z <- !z
			}
		else if(rel == "ge"){
			z <- new("matrix.csr",ra=z$ra[1:nnza],ja=z$ja[1:nnza],ia=z$ia,dimension=A@dimension)
			z <- !z
			}
		else if(rel == "ne")
			z <- new("matrix.csr",ra=-z$ra[1:nnza],ja=z$ja[1:nnza],ia=z$ia,dimension=A@dimension)
		else if(drptol == 0){
			z <- new("matrix.csr",ra=-z$ra[1:nnza],ja=z$ja[1:nnza],ia=z$ia,dimension=A@dimension)
			z <- !z
			}
		else
			z <- new("matrix.csr",ra=-z$ra[1:nnza],ja=z$ja[1:nnza],ia=z$ia,dimension=A@dimension)
		}
	else{ #This operation is inefficient storage-wise
		if(rel == "gt")
			z <- as.matrix.csr(as.matrix(A) > drptol)
		else if(rel == "ge")
			z <- as.matrix.csr(as.matrix(A) >= drptol)
		else if(rel == "lt")
			z <- as.matrix.csr(as.matrix(A) < drptol)
		else if(rel == "le")
			z <- as.matrix.csr(as.matrix(A) <= drptol)
		else
			z <- as.matrix.csr(as.matrix(A) != drptol)
		}
	z
}
#--------------------------------------------------------------------
"[.matrix.csr" <- function (x, rw = 1:x@dimension[1], cl = 1:x@dimension[2])
{
        x <- as.matrix.coo(x)
	y <- x[rw,cl]
	if(is(y,"matrix.coo"))
		as.matrix.csr(y)
	else
		y
}
#--------------------------------------------------------------------
"[<-.matrix.csr" <- function (x, rw = 1:x@dimension[1], cl = 1:x@dimension[2], value)
{
        x <- as.matrix.coo(x)
#        value <- as.matrix.coo(value)
        x[rw,cl]<-value
        as.matrix.csr(x)
}
#--------------------------------------------------------------------
"[.matrix.diag.csr" <- function (x, rw = 1:x@dimension[1], cl = 1:x@dimension[2])
{
        x <- as.matrix.coo(as(x,"matrix.csr"))
	y <- x[rw,cl]
	if(is(y,"matrix.coo"))
		as(as.matrix.csr(y),"matrix.diag.csr")
	else
		y
}
#--------------------------------------------------------------------
".matmul.matrix.csr" <- function(x,y){
if(is.matrix.csr(x)){
	if(is.matrix.csr(y)){
		#matrix multiply two sparse csr matrices 
		nrow <- x@dimension[1]
		ncol <- y@dimension[2]
		Acol <- x@dimension[2]
		Brow <- y@dimension[1]
		if(Acol != Brow)
			stop("matrices not conformable for multiplication")
		z <- .Fortran("amubdg",
			as.integer(nrow),
			as.integer(Acol),
			as.integer(ncol),
			as.integer(x@ja),
			as.integer(x@ia),
			as.integer(y@ja),
			as.integer(y@ia),
			integer(nrow),
			nnz = integer(1),
			integer(ncol),
			PACKAGE = "SparseM")
		nnzmax <- z$nnz
		z <- .Fortran("amub",
		   	as.integer(nrow),
		   	as.integer(ncol),
		   	as.integer(1),
		   	as.double(x@ra),
		   	as.integer(x@ja),
		   	as.integer(x@ia),
		   	as.double(y@ra),
		   	as.integer(y@ja),
		   	as.integer(y@ia),
		   	ra = double(nnzmax),
		   	ja = integer(nnzmax),
		   	ia = integer(nrow+1),
		   	as.integer(nnzmax),
		   	integer(ncol),
		   	ierr = integer(1),
			PACKAGE = "SparseM")
		nnz <- z$ia[nrow+1]-1
		if(z$ierr != 0) stop("insufficient space for sparse matrix multiplication")
		if(length(z$ra)==0){#trap zero matrix
			z$ra <- 0
			z$ja <- as.integer(1)
			z$ia <- as.integer(c(1,rep(2,nrow)))
			}
		z <- new("matrix.csr",ra=z$ra[1:nnz],ja=z$ja[1:nnz],
			ia=z$ia,dimension=as.integer(c(nrow,ncol)))
	   	}
	else{
		if(is.matrix(y)){
			z <- .matmul.matrix.csr(x,as.matrix.csr(y))
			}
		else{
		#matrix-vector multiplication: multiply a sparse csr matrix by a vector
		#A -- csr structure returned from call to function "as.matrix.csr"
		#B -- vector
			nrow <- x@dimension[1]
			ncol <- x@dimension[2]
			if(length(y) != ncol)stop("not conformable for multiplication")
			z <- .Fortran("amux",
		   		as.integer(nrow),
				as.double(y),
			   	y=double(nrow),
			   	as.double(x@ra),
		   		as.integer(x@ja),
		   		as.integer(x@ia),
				PACKAGE = "SparseM")
			z <- z$y
			dim(z) <- c(nrow,1)
			}
		 }
	}
else{
	if(is.matrix(x)){
		z <- .matmul.matrix.csr(as.matrix.csr(x),y)
		}
	else{
	#matrix-vector multiplication: multiply a sparse csr matrix by a vector
	#A -- csr structure returned from call to function "as.matrix.csr"
	#B -- vector
		y <- t(y)
		nrow <- y@dimension[1]
		ncol <- y@dimension[2]
		if(length(x) != ncol)stop("not conformable for multiplication")
		z <- .Fortran("amux",
   			as.integer(nrow),
			as.double(x),
		   	y=double(nrow),
		   	as.double(y@ra),
   			as.integer(y@ja),
   			as.integer(y@ia),
			PACKAGE = "SparseM")
		z <- z$y
		dim(z) <- c(1,nrow)
		}
	}
return(z)
}
#--------------------------------------------------------------------
".kron.matrix.csr" <-
function(X,Y){
	X = as.matrix.coo(X)
	Y = as.matrix.coo(Y)

	la = length(X@ra)
	lb = length(Y@ra)

	ra = rep(Y@ra,la)*rep(X@ra,each=lb)
	ja = as.integer(rep(Y@ja,la)+rep((X@ja-1)*dim(Y)[2],each=lb))
	ia = as.integer(rep(Y@ia,la)+rep((X@ia-1)*dim(Y)[1],each=lb))
	dim = as.integer(dim(X)*dim(Y))

	as.matrix.csr(new("matrix.coo",ra=ra,ia=ia,ja=ja, dimension=dim))
	# /*RSB*/ dim= changed to dimension=
}

#--------------------------------------------------------------------
"chol" <- function(x, ...) UseMethod("chol")
#--------------------------------------------------------------------
"chol.default" <-get("chol", pos=NULL, mode= "function")
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
    class(fit) <- c(if (is.matrix(Y)) "mslm", "slm")
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
    x <- x %*% as(w,"matrix.diag.csr")
    y <- y * w
    fit <- slm.fit.csr(x, y,  ...)
    fit$contrasts <- attr(x, "contrasts")
    fit
}
#--------------------------------------------------------------------
"slm.fit.csr" <-
function (x, y, ...) 
{
#    n <- length(y)
    if(is.matrix(y))
	n <- dim(y)[1]
    else
        n <- length(y)
    p <- x@dimension[2]
    if (n != x@dimension[1]) 
        stop("x and y don't match n")
    chol <- chol(t(x)%*%x)
    xy <- t(x) %*% y
    coef <- backsolve(chol,xy)
    fitted <-  as.matrix(x %*% coef)
    resid <- y - fitted
    list(coefficients=coef, chol=chol, residuals=resid,  fitted=fitted)
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
"residuals.slm" <- function (object, ...){
r <- object$residuals
r
}
#--------------------------------------------------------------------
"summary.mslm" <- function (object, ...)
{
    coef <- coef(object)
    ny <- ncol(coef)
    if (is.null(ny))
        return(NextMethod("summary"))
    effects <- object$effects
    resid <- residuals(object)
    fitted <- fitted(object)
    ynames <- colnames(coef)
    if (is.null(ynames)) {
        lhs <- object$terms[[2]]
        if (mode(lhs) == "call" && lhs[[1]] == "cbind")
            ynames <- as.character(lhs)[-1]
        else ynames <- paste("Y", seq(ny), sep = "")
    }
    value <- vector("list", ny)
    names(value) <- paste("Response", ynames)
    cl <- class(object)
    class(object) <- cl[match("mslm", cl):length(cl)][-1]
    for (i in seq(ny)) {
        object$coefficients <- coef[, i]
        object$residuals <- resid[, i]
        object$fitted.values <- fitted[, i]
        object$effects <- effects[, i]
        object$call$formula[[2]] <- object$terms[[2]] <- as.name(ynames[i])
        value[[i]] <- summary(object, ...)
    }
    class(value) <- "listof"
#    class(value) <- "summary.mslm"
    value
}
#--------------------------------------------------------------------
"summary.slm" <-
function (object, correlation = FALSE, ...)
{
    Chol <- object$chol
    if (is.null(object$terms) || is.null(Chol))
        stop("invalid 'lm' object:  no terms or chol component")
    n <- length(object$residuals)
    p <- object$chol@nrow
    rdf <- n - p
    r <- residuals(object)
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
        rq <- structure(quantile(resid), names = nam)
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
    printCoefmat(x$coef, digits = digits, signif.stars = signif.stars,
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
"print.slm" <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("\nCall:\n", deparse(x$call), "\n\n", sep = "")
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2, 
        quote = FALSE)
    cat("\n")
    invisible(x)
}
"[.matrix.coo" <- function (x, rw = 1:x@dimension[1], cl = 1:x@dimension[2])
{
	nrow <-  x@dimension[1]
	ncol <-  x@dimension[2]
	if (is.matrix.coo(rw)) {
		stop("Indexing by 'matrix.coo' matrices not implemented")
		}
	if (is.matrix.csr(rw)) {
		if(nrow != rw@dimension[1] || ncol != rw@dimension[2])
			stop("Dimension of the indexing matrix is not the same as the matrix being indexed")
		x <- as.matrix.coo(t(as.matrix.csr(x)))
		rw <- as.matrix.coo(t(rw))
		s <- match(paste(rw@ia,rw@ja),paste(x@ia,x@ja))
		ra <- x@ra[s]
		A <- ra
	        }
	else if(is.matrix(rw)){
		if(!all(abs(rw[,1])%in%1:nrow)||!all(abs(rw[,2])%in%1:ncol))
			stop("Subscripts out of bound")
		s <- match(paste(rw[,1],rw[,2]),paste(x@ia,x@ja))
		ra <- x@ra[s]
		ra[is.na(ra)] <- 0
		A <- ra
		}
	else{
		if (is.logical(rw))
			if (length(rw) > nrow)
				stop("logical subscript too long")
			else
				rw <- (1:nrow)[rw]
		if (is.logical(cl))
			if (length(cl) > ncol)
				stop("logical subscript too long")
			else
				cl <- (1:ncol)[cl]
       		if (!all(abs(rw) %in% 1:nrow) || !all(abs(cl) %in% 1:ncol))
			stop("Subscripts out of bound")
	        if (any(rw < 0)) {
			if (!all(rw <= 0))
		                stop("Only 0's may mix with negative subscripts")
			else rw <- setdiff(1:nrow, abs(rw))
		}
		if (any(cl < 0)) {
			if (!all(cl <= 0))
				stop("Only 0's may mix with negative subscripts")
			else cl <- setdiff(1:ncol, abs(cl))
		}
        	if(any(duplicated(rw)) || any(duplicated(cl))){
               		s <- ((x@ja %in% cl) & (x@ia %in% rw))
                	urw <- as.integer(names(table(rw)))
                	ucl <- as.integer(names(table(cl)))
                	ja <- match(x@ja[s],ucl)
                	ia <- match(x@ia[s],urw)
                	dim <- c(length(urw),length(ucl))
                	ra <- x@ra[s]
                	A <- new("matrix.coo",ra=ra,ja=ja,ia=ia,dimension=dim)
			A <- as.matrix.csr(A)

                	#obviously this looping is horrible but is there a better way?
                	#could be fortranized I suppose.
	
                	rn <- table(rw)
                	ri <- as.integer(names(rn))
                	B <- NULL
                	for(i in 1:length(rw)){
				j <- match(rw[i],urw)
				if(is.null(B))
					B <- A[j,]
				else
                               		B <- rbind(B,A[j,])
                        	}
                	A <- B
                	B <- NULL
                	cn <- table(cl)
                	ci <- as.integer(names(cn))
                	for(i in 1:length(cl)){
				j <- match(cl[i],ucl)
				if(is.null(B))
					B <- A[,j]
				else
                               		B <- cbind(B,A[,j])
                        	}
                	A <- as.matrix.coo(B)
                	}
        	else{
                	s <- ((x@ja %in% cl) & (x@ia %in% rw))
                	ja <- match(x@ja[s],cl)
                	ia <- match(x@ia[s],rw)
                	dim <- c(length(rw),length(cl))
                	ra <- x@ra[s]
			if (length(ra) == 0 && length(ja) == 0){ 
				ra = 0
				ja = ia = as.integer(1)
				}
                	A <- new("matrix.coo",ra=ra,ja=ja,ia=ia,dimension=dim)
                	}
	
	}
	return(A)
}
"[<-.matrix.coo" <-
function (x, rw = 1:x@dimension[1], cl = 1:x@dimension[2], value) 
{
    nrow <-  x@dimension[1]
    ncol <-  x@dimension[2]
#    if (!is.matrix.coo(value)) 
#	 stop("replacement matrix must be of matrix.coo class")
    if (is.matrix.coo(rw)) {
         stop("Indexing by 'matrix.coo' matrices not implemented")
	}
    else if (is.matrix.csr(rw)) {
	if(nrow != rw@dimension[1] || ncol != rw@dimension[2])
		stop("Dimension of the indexing matrix is not the same as the matrix being indexed")
	x <- as.matrix.coo(t(as.matrix.csr(x)))
	rw <- as.matrix.coo(t(rw))
	s <- match(paste(rw@ia,rw@ja),paste(x@ia,x@ja))
	len.s <- length(s)
	value <- as.matrix.csr(value)
	len.value <- length(value@ra)
	value.ra <- value@ra
	if (len.s != len.value){
		if (len.value == 1) 
			value.ra <- rep(value@ra,len.s)
		else{
			if (len.value >= len.s){
				value.ra <- value@ra[1:len.s]
				warning("number of items to replace is not a multiple of replacement length")
				}
			else {
				value.ra <- rep(value@ra,len.s%/%len.value+1)[1:len.s]
				warning("number of items to replace is not a multiple of replacement length")
				}
			}
		}
	ra <- c(x@ra[-s],value.ra)
        ia <- as.integer(c(x@ja[-s],rw@ja))
        ja <- as.integer(c(x@ia[-s],rw@ia))
        dim <- rev(x@dimension)
        x <- new("matrix.coo",ra=ra,ja=ja,ia=ia,dimension=dim)
        }
    else if(is.matrix(rw)){
	if(!all(abs(rw[,1])%in%1:nrow)||!all(abs(rw[,2])%in%1:ncol))
		stop("Subscripts out of bound")
	s <- match(paste(rw[,1],rw[,2]),paste(x@ia,x@ja))
#	ra <- c(x@ra[-s],value@ra[!is.na(s)])
#        ja <- as.integer(c(x@ja[-s],rw[,2][!is.na(s)]))
#        ia <- as.integer(c(x@ia[-s],rw[,1][!is.na(s)]))
	value <- as.matrix.csr(value)
	ra <- c(x@ra[-s],value@ra)
        ja <- as.integer(c(x@ja[-s],rw[,2]))
        ia <- as.integer(c(x@ia[-s],rw[,1]))
        dim <- x@dimension
        x <- new("matrix.coo",ra=ra,ja=ja,ia=ia,dimension=dim)
	}
    else {
	if (is.logical(rw))
		if (length(rw) > nrow)
			stop("logical subscript too long")
		else
			rw <- (1:nrow)[rw]
	if (is.logical(cl))
		if (length(cl) > ncol)
			stop("logical subscript too long")
		else
			cl <- (1:ncol)[cl]
        if (!all(abs(rw) %in% 1:nrow) || !all(abs(cl) %in% 1:ncol))
            stop("Subscripts out of bound")
        if (any(rw < 0)) {
            if (!all(rw <= 0))
                stop("Only 0's may mix with negative subscripts")
            else rw <- setdiff(1:nrow, abs(rw))
        }
        if (any(cl < 0)) {
            if (!all(cl <= 0))
                stop("Only 0's may mix with negative subscripts")
            else cl <- setdiff(1:ncol, abs(cl))
        }
        s <- ((x@ja %in% cl) & (x@ia %in% rw))
	value <- as.matrix.coo(as.matrix.csr(value,nrow=length(rw),ncol=length(cl)))
        ra <- c(x@ra[!s],value@ra)
        ja <- as.integer(c(x@ja[!s],cl[value@ja]))
        ia <- as.integer(c(x@ia[!s],rw[value@ia]))
        dim <- x@dimension
        x <- new("matrix.coo",ra=ra,ja=ja,ia=ia,dimension=dim)
        }
return(x)
}
# All the S4 Methods stuff is collected below this point
#require(methods) /*RSB*/ commented out
setClass("matrix.csr",representation(ra="numeric",
	ja="integer",ia="integer", dimension="integer"),
	validity = function(object) {
 		if(!(length(object@dimension) == 2) )
                	return("invalid dimension attribute")
        	else{
               		nrow <- object@dimension[1]
               		ncol <- object@dimension[2]
               	 	}
        	if(!(length(object@ra) ==length(object@ja)))
                	return("ra and ja don't have equal lengths")
        	if(any(object@ja < 1) || any(object@ja > ncol)) 
                	return("ja exceeds dim bounds")
        	if(any(object@ia < 1))
                	return("some elements of ia are <= 0")
		if(any(diff(object@ia)<0))
			return("ia vector not monotone increasing") 
		if(object@ia[length(object@ia)] != length(object@ra)+1)
			return("last element of ia doesn't conform")
		if(length(object@ia) != nrow+1)
			return("ia has wrong number of elements")
        	if(length(object@ra) < 1 || length(object@ra) > prod(object@dimension))
                	return("ra has too few, or too many elements")
		TRUE})
setMethod("initialize", "matrix.csr", 
	function(.Object, ra = 0,  ja = as.integer(1),
		ia = as.integer(c(1,2)),dimension = as.integer(c(1,1))) {
        .Object@ra <- ra
        .Object@ja <- ja
        .Object@ia <- ia
        .Object@dimension <- dimension
	validObject(.Object)
        .Object
       })
setClass("matrix.csc",representation(ra="numeric",ja="integer",ia="integer", dimension="integer"))
setClass("matrix.ssr",representation(ra="numeric",ja="integer",ia="integer", dimension="integer"))
setClass("matrix.ssc",representation(ra="numeric",ja="integer",ia="integer", dimension="integer"))
setClass("matrix.coo",representation(ra="numeric",
	ja="integer",ia="integer", dimension="integer"),
	validity = function(object) {
 		if(!length(object@dimension) == 2 )
                	return("invalid dimension attribute")
        	else{
               		nrow <- object@dimension[1]
               		ncol <- object@dimension[2]
               	 	}
        	if(!(length(object@ra) ==length(object@ja) && length(object@ra) ==length(object@ia)))
                	return("ra,ja,ia don't have equal lengths")
        	if(any(object@ja < 1) || any(object@ja > ncol)) 
                	return("ja exceeds dim bounds")
        	if(any(object@ia < 1) || any(object@ia > nrow))
                	return("ia exceeds dim bounds")
        	if(length(object@ra) < 1 || length(object@ra) > prod(object@dimension))
                	return("ra has too few, or too many elements")
		TRUE})
setMethod("initialize", "matrix.coo", 
	function(.Object, ra = numeric(0),  ja = integer(0),
		ia = integer(0),dimension = integer(0)) {
        .Object@ra <- ra
        .Object@ja <- ja
        .Object@ia <- ia
        .Object@dimension <- dimension
	validObject(.Object)
        .Object
       })
#-------------------------------------------------------------------------
setClass("matrix.csr.chol",representation(nrow="numeric",nnzlindx="numeric",
	nsuper="numeric",lindx="numeric",xlindx="numeric",nnzl="numeric",
	lnz="numeric",xlnz="numeric",invp="numeric",perm="numeric",
	xsuper="numeric",det="numeric",ierr="numeric",time="numeric"))
if(version$major >= 1 && version$minor >= 8.0 || version$major >= 2){
	setClassUnion("numeric or NULL",c("numeric","NULL"))
	setClassUnion("character or NULL",c("character","NULL"))
	}
if(version$major == 1 && version$minor < 8.0){
	setClass("numeric or NULL")
	setIs("numeric","numeric or NULL")
	setIs("NULL","numeric or NULL")
	setClass("character or NULL")
	setIs("character","character or NULL")
	setIs("NULL","character or NULL")
	}
setClass("matrix.csc.hb",representation(ra="numeric",ja="integer",ia="integer", 
	rhs.ra="numeric",guess="numeric or NULL",xexact="numeric or NULL",dimension ="integer",
	rhs.dim="numeric",rhs.mode="character or NULL"))
setClass("matrix.ssc.hb","matrix.csc.hb")
setClass("slm",representation(coefficients="numeric",chol="matrix.csr.chol",
	residuals="numeric",fitted="numeric"))
setClass("mslm","slm")
setClass("summary.slm","slm")
#--------------------------------------------------------------------
setGeneric("as.matrix.csr")
#--------------------------------------------------------------------
setMethod("as.matrix.csr","matrix.csc", function(x, nrow, ncol,eps){
	x <- t(x)
        x@dimension <- as.integer(rev(dim(x)))
        new("matrix.csr",ra = x@ra, ja = x@ja, ia = x@ia, dimension = x@dimension)
        })
#--------------------------------------------------------------------
setMethod("as.matrix.csr","matrix.ssr", function(x, nrow, ncol,eps){.ssr.csr(x)})
#--------------------------------------------------------------------
setMethod("as.matrix.csr","matrix.ssc", function(x, nrow, ncol,eps){.ssr.csr(x)})
#--------------------------------------------------------------------
setMethod("as.matrix.csr","matrix.coo", function(x, nrow, ncol,eps){
#       if (missing(nrow)) nrow <- x@dimension[1]
#       if (missing(ncol)) ncol <- x@dimension[2]
        nrow <- x@dimension[1]
        ncol <- x@dimension[2]
        nnz <- length(x@ra)
        z <- .Fortran("coocsr",
                as.integer(nrow),
                as.integer(nnz),
                as.double(x@ra),
                as.integer(x@ia),
                as.integer(x@ja),
                ao = double(nnz),
                jao = integer(nnz),
                iao = integer(nrow+1),
                PACKAGE = "SparseM")
        nnza <- z$ao[nrow+1]-1
        z <- new("matrix.csr",ra=z$ao,ja=z$jao,ia=z$iao,dimension=x@dimension)
        return(z)
})
#--------------------------------------------------------------------
setGeneric("as.matrix.csc")
#--------------------------------------------------------------------
setMethod("as.matrix.csc","matrix.csr", function(x, nrow, ncol,eps){
	x <- t(x)
        x@dimension <- as.integer(rev(dim(x)))
        new("matrix.csc",ra = x@ra, ja = x@ja, ia = x@ia, dimension = x@dimension)
	})
#--------------------------------------------------------------------
setMethod("as.matrix.csc","matrix.ssr", function(x, nrow, ncol,eps){
	as.matrix.csc(as.matrix.csr(x))})
#--------------------------------------------------------------------
setMethod("as.matrix.csc","matrix.ssc", function(x, nrow, ncol,eps){
	as.matrix.csc(as.matrix.csr(x))})
#--------------------------------------------------------------------
setGeneric("as.matrix.ssr")
#--------------------------------------------------------------------
setMethod("as.matrix.ssr","matrix.csr", function(x, nrow, ncol,eps){.csr.ssr(x)})
#--------------------------------------------------------------------
setMethod("as.matrix.ssr","matrix.csc", function(x, nrow, ncol,eps){
	as.matrix.ssr(as.matrix.csr(x))})
#--------------------------------------------------------------------
setMethod("as.matrix.ssr","matrix.ssc", function(x, nrow, ncol,eps){
	as.matrix.ssr(as.matrix.csr(x))})
#--------------------------------------------------------------------
setGeneric("as.matrix.ssc")
#--------------------------------------------------------------------
setMethod("as.matrix.ssc","matrix.csr", function(x, nrow, ncol,eps){
	x <- as.matrix.csc(x)
        as.matrix.ssc(x)})
#--------------------------------------------------------------------
setMethod("as.matrix.ssc","matrix.csc", function(x, nrow, ncol,eps){.csc.ssc(x)})
#--------------------------------------------------------------------
setMethod("as.matrix.ssc","matrix.ssr", function(x, nrow, ncol,eps){
	as.matrix.ssc(as.matrix.csr(x))})
#--------------------------------------------------------------------
setGeneric("as.matrix.coo")
#--------------------------------------------------------------------
setMethod("as.matrix.coo","matrix.csr", function(x, nrow, ncol,eps){.csr.coo(x)})
#--------------------------------------------------------------------
setMethod("as.matrix","matrix.csr", function(x){
	nrow <- x@dimension[1]
        ncol <- x@dimension[2]
        if(length(x@ra)==1 && is.na(x@ra)){ #trap zero matrix
                dns <- matrix(0,nrow=nrow,ncol=ncol)
                return(dns)
                }
        nan <- is.nan(x@ra)
        infty <- is.infinite(x@ra) & x@ra >0
        ninfty <- is.infinite(x@ra) & x@ra <0
        uniq <- rnorm(3)
        while(any(uniq %in% x@ra[!(nan|infty|ninfty)]))
                uniq <- rnorm(3)
        x@ra[nan] <- uniq[1]
        x@ra[infty] <- uniq[2]
        x@ra[ninfty] <- uniq[3]
        z <- .Fortran("csrdns",
                as.integer(nrow),
                as.integer(ncol),
                as.double(x@ra),
                as.integer(x@ja),
                as.integer(x@ia),
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
	})
setMethod("as.matrix","matrix.ssr",function(x)as.matrix(as.matrix.csr(x)))
setMethod("as.matrix","matrix.ssc",function(x)as.matrix(as.matrix.csr(x)))
setMethod("as.matrix","matrix.coo",function(x)as.matrix(as.matrix.csr(x)))
setMethod("t","matrix.csr",function(x){
	nrow <- x@dimension[1]
        ncol <- x@dimension[2]
        nnz <- x@ia[nrow+1]-1
        z <- .Fortran("csrcsc2", as.integer(nrow), as.integer(ncol),
                as.integer(1), as.integer(1), as.double(x@ra), as.integer(x@ja),
                as.integer(x@ia), ao=double(nnz), jao=integer(nnz), iao=integer(ncol+1),
                PACKAGE = "SparseM")
        dim <- as.integer(rev(x@dimension))
        new("matrix.csr",ra = z$ao, ja = z$jao, ia = z$iao, dimension = dim)
	})
setMethod("t","matrix.csc",function(x) {
    nrow <- x@dimension[1]
    ncol <- x@dimension[2]
    nnz <- x@ia[ncol + 1] - 1
    z <- .Fortran("csrcsc2", as.integer(ncol), as.integer(nrow),
        as.integer(1), as.integer(1), as.double(x@ra), as.integer(x@ja),
        as.integer(x@ia), ao = double(nnz), jao = integer(nnz),
        iao = integer(nrow + 1), PACKAGE = "SparseM")
    dim <- as.integer(rev(x@dimension))
    new("matrix.csc",ra = z$ao, ja = z$jao, ia = z$iao, dimension = dim)
   })
setMethod("t","matrix.coo",function(x) as.matrix.coo(t(as.matrix.csr(x))))
setMethod("dim","matrix.csr",function(x)x@dimension)
setMethod("dim","matrix.csc",function(x)x@dimension)
setMethod("dim","matrix.ssr",function(x)x@dimension)
setMethod("dim","matrix.ssc",function(x)x@dimension)
setMethod("dim","matrix.coo",function(x)x@dimension)
setMethod("diff","matrix.csr", function(x, lag = 1, differences = 1, ...) {
    xlen <- dim(x)[1]
    if (length(lag) > 1 || length(differences) > 1 || lag < 1 ||
        differences < 1)
        stop("`lag' and `differences' must be integers >= 1")
    if (lag * differences >= xlen)
        stop("lag * differences >= nrow")
     r <- x
    i1 <- -1:-lag
        for (i in 1:differences) r <- r[i1,] - r[-nrow(r):-(nrow(r) - lag + 1),]
    r
   })
#setClass("matrix.diag.csr")
setClass("matrix.diag.csr","matrix.csr")
#setIs("matrix.csr","matrix.diag.csr")
setAs("numeric","matrix.diag.csr",function(from){
	#if(!is.numeric(from))stop("non-numeric entries in sparse matrices not allowed")
	if(length(from)==1){
		n <- as.integer(from)
		if(n>0) from  <-  rep(1,n)
		else stop("Sparse identity matrices must have positive, integer dimension") 
		}
	else n <- length(from)
	return(new("matrix.diag.csr", ra = from ,ja = as.integer(1:n), 
		ia = as.integer(1:(n+1)), dimension = as.integer(c(n,n))))
	})
setAs("matrix.csr","matrix.diag.csr",function(from){
        nr <- from@dimension[1]
        nc <- from@dimension[2]
        if( nr!=nc) stop("Resulting 'matrix.diag.csr' matrix has to be square")
        if(!(setequal(from@ja,1:nc)&setequal(from@ia, 1:(nr+1)))) stop("matrix not diagonal")
        new("matrix.diag.csr", ra = from@ra, ja = from@ja, ia = from@ia, dimension = from@dimension)
        })
setMethod("diag","matrix.csr", function (x = 1, nrow, ncol=n){
    if (is.matrix.csr(x) && nargs() == 1) {
        if ((m <- min(dim(x))) == 0)
            return(numeric(0))
        y <- rep(0,m)
        ia <- rep(1:dim(x)[1],diff(x@ia))
        y[x@ja[ia == x@ja]] <- x@ra[ia == x@ja]
        n <- sum(ia == x@ja)
        nms <- dimnames(x)
        if (is.list(nms) && !any(sapply(nms, is.null)) && all((nm <- nms[[1]][1:m]) ==
            nms[[2]][1:m]))
            names(y) <- nm
        return(y)
	}
   else stop("diag method for class matrix.csr doesn't understand nrow and ncol args")
   })
setMethod("diag<-","matrix.csr", function(x,value) {
     dx <- dim(x)
     if (length(dx) != 2)
         stop("only matrix diagonals can be replaced")
     i <- seq(length = min(dx))
     if (length(value) != 1 && length(value) != length(i))
         stop("replacement diagonal has wrong length")
     if (length(value) == 1) value <- rep(value,min(dx))
     x[cbind(i, i)] <- value
     x
     })
setMethod("diag<-","matrix.diag.csr",function(x,value) {
	y <- as(x,"matrix.csr")
	diag(y) <- value
	as(y,"matrix.diag.csr")
	})
setGeneric("det")
setMethod("det","matrix",get("det", pos=NULL, mode= "function"))
setMethod("det","matrix.csr", function(x, ...) det(chol(x))^2)
setMethod("det","matrix.csr.chol", function(x, ...) x@det)
setGeneric("norm",function(x, ...)standardGeneric("norm"))
setMethod("norm","matrix.csr", function(x, type = "sup", ...){
  switch(type,
                sup = max(abs(x@ra)),
                HS = sqrt(sum(x@ra^2)),
                l1 = sum(abs(x@ra))
                )
        }
)
setGeneric("chol")
setMethod("chol","matrix",get("chol", pos=NULL, mode= "function"))
setMethod("chol","matrix.csr", function(x, pivot = FALSE, 
	nsubmax, nnzlmax, tmpmax, eps = .Machine$double.eps, ...){
# Interface for a sparse least squares solver via Ng-Peyton's Cholesky
# factorization
#       x -- csr structure returned from call to function "as.matrix.csr"
#       cachsz -- size of the cache memory -- machine dependent.
# Check that input matrix is symmetric
	if(norm(t(x)-x) > eps) stop("Input matrix to chol() not symmetric")
        cachsz <- 64
        nrow <- x@dimension[1]
        ncol <- x@dimension[2]
        if(nrow!=ncol) stop("Can't perform Cholesky Factorization for Non-square matrix\n")
        nnzdmax <- x@ia[nrow+1]-1
        nnzdsm <- nnzdmax + nrow + 1
        iwmax <- 7*nrow+3
        if(missing(nsubmax)) nsubmax <- nnzdmax
        if(missing(nnzlmax)) nnzlmax <- max(4*nnzdmax,floor(.2*nnzdmax^1.3))
        if(missing(tmpmax)) tmpmax <- 10*nrow
        level <- 8
        z <- .Fortran("chol", nrow = as.integer(nrow), nnzdmax = as.integer(nnzdmax),
                d = as.double(x@ra), jd = as.integer(x@ja), id = as.integer(x@ia),
                nnzdsm = as.integer(nnzdsm), dsub = double(nnzdsm), jdsub = integer(nnzdsm),
                nsub = integer(1), nsubmax = as.integer(nsubmax), lindx = integer(nsubmax),
                xlindx = integer(nrow+1), nsuper = integer(1), nnzlmax = as.integer(nnzlmax),
                lnz = double(nnzlmax), xlnz = integer(nrow+1), invp = integer(nrow),
                perm = integer(nrow), iwmax = as.integer(iwmax), iwork = integer(iwmax),
                colcnt = integer(nrow), snode = integer(nrow), xsuper = integer(nrow+1),
                split = integer(nrow), tmpmax = as.integer(tmpmax), tmpvec = double(tmpmax),
                cachsz = as.integer(cachsz), level = as.integer(level), ierr = integer(1),
                time = double(1), PACKAGE = "SparseM") 
if (z$ierr != 0){
                if(z$ierr == 9) mess <- "singularity problem"
                else if(z$ierr == 4) mess <- "Increase nnzlmax"
                else if(z$ierr == 5) mess <- "Increase nsubmax"
                else if(z$ierr %in% c(8,10)) mess <- "Increase tmpmax"
                else mess <- "insufficient space"
                if(z$ierr == 9) warning(mess) else stop(mess)
                }
        nnzl <- z$xlnz[length(z$xlnz)]-1
        nnzlindx <- z$nsub
	k <- z$xlnz
	R <- z$lnz
	det <- prod(R[k[-length(k)]])
        new("matrix.csr.chol",nrow=z$nrow,nnzlindx=nnzlindx,
                nsuper=z$nsuper,lindx=z$lindx[1:nnzlindx],xlindx=z$xlindx,
                nnzl=as.integer(nnzl),lnz=z$lnz[1:nnzl],xlnz=z$xlnz,invp=z$invp,
                perm=z$perm,xsuper=z$xsuper,det=det,ierr=z$ierr,time=z$time)
	})
setMethod("chol","matrix.csc",function(x,pivot = FALSE, ...)chol(as.matrix.csr(x)))
setMethod("backsolve","matrix.csr.chol", 
function(r, x, k = NULL, upper.tri = NULL, transpose = NULL){
# backsolve for Ng-Peyton's Cholesky factorization
#       Solves linear system A b = x where r is chol(A)
# Input:
#       r -- structure returned by chol.matrix.csr
#       x --  rhs  may be a matrix in dense form
        m <- r@nrow
        if(!is.matrix(x)) x <- as.matrix(x)
        if(nrow(x)!=m)stop("chol not conformable with x")
        p <- ncol(x)
        z <- .Fortran("bckslv", m = as.integer(m), nnzlindx = as.integer(r@nnzlindx),
                as.integer(r@nsuper), as.integer(p), as.integer(r@lindx),
                as.integer(r@xlindx), as.integer(r@nnzl), as.double(r@lnz),
                as.integer(r@xlnz), as.integer(r@invp), as.integer(r@perm),
                as.integer(r@xsuper), double(m), sol = double(m*p), as.double(x),
                time = double(1), PACKAGE = "SparseM")
        z <- matrix(z$sol,nrow=nrow(x),ncol=ncol(x))
        drop(z)
	})
setMethod("solve","matrix.csr", function (a, b, ...) {
    missing.b <- FALSE
    if(!is.matrix.csr(a))stop("a not in csr format")
    nr <- nrow(a)
    nc <- ncol(a)
    if (nc != nr)
         stop("only square systems can be solved")
    a <- chol(a,...)
    if (missing(b)) {
        b <- diag(1, nc)
        missing.b <- TRUE
    }
    else
        if(!is.matrix.csr(b))b <- as.matrix(b)
    z <- backsolve(a,b)
    if (missing.b)
        z <- as.matrix.csr(z)
    z
   })
setGeneric("model.matrix", function(object, ...) # /*RSB*/ changed definition
	standardGeneric("model.matrix")) # /*RSB*/
	
setMethod("model.matrix","matrix.csc.hb", function(object,...){
        object <- new("matrix.csc",ra=object@ra,ja=object@ja,
                ia=object@ia,dimension=object@dimension)
        as.matrix.csr(object)
	})
setMethod("model.matrix","matrix.ssc.hb", function(object, ...){
	object <- new("matrix.ssc",ra=object@ra,ja=object@ja,ia=object@ia,
		dimension=object@dimension)
        as.matrix.csr(object)
	})

setGeneric("model.response", function(data, type) # /*RSB*/ changed definition
	standardGeneric("model.response")) # /*RSB*/
setMethod("model.response","ANY", # /*RSB*/
function(data,type="any"){ # /*RSB*/
	stats:::model.response(data, type="any") # /*RSB*/
	}) # /*RSB*/


setMethod("model.response","matrix.csc.hb",
function(data,type="any"){  
        if(is.null(data@rhs.mode)) stop("Right-hand side doesn't exist")
        if (data@rhs.mode == "F")
                z <- data@rhs.ra
        else{
                z <- new("matrix.csc",ra=data@rhs.ra,ja=data@rhs.ja,ia=data@rhs.ia,
			dimension=data@rhs.dim)
                z <-  as.matrix.csr(z)
                }
        z
	})
setMethod("model.response","matrix.ssc.hb",
function(data,type="any"){
        data <- new("matrix.ssc",ra=data@ra,ja=data@ja,ia=data@ia,
		dimension=data@dimension)
        as.matrix.csr(data)
	})
setMethod("%*%",signature(x="matrix.csr",y="matrix.csr"),.matmul.matrix.csr)
setMethod("%*%",signature(x="matrix.csr",y="matrix"),.matmul.matrix.csr)
setMethod("%*%",signature(x="matrix.csr",y="numeric"),.matmul.matrix.csr)
setMethod("%*%",signature(x="matrix",y="matrix.csr"),.matmul.matrix.csr)
setMethod("%*%",signature(x="numeric",y="matrix.csr"),.matmul.matrix.csr)
if(version$major == 1 && version$minor >= 8.0 || version$major >= 2){
	setMethod("%x%",signature(X="matrix.csr",Y="matrix.csr"),
		function(X,Y) .kron.matrix.csr(X,Y))
	setMethod("%x%",signature(X="matrix.csr",Y="numeric"),
		function(X,Y) .kron.matrix.csr(X,Y))
	setMethod("%x%",signature(X="numeric",Y="matrix.csr"),
		function(X,Y) .kron.matrix.csr(X,Y))
	setMethod("%x%",signature(X="matrix",Y="matrix.csr"),
		function(X,Y) .kron.matrix.csr(X,Y))
	setMethod("%x%",signature(X="matrix.csr",Y="matrix"),
		function(X,Y) .kron.matrix.csr(X,Y))
        }
if(version$major == 1 && version$minor < 8.0 ){
	setMethod("%x%",signature(X="matrix.csr",Y="matrix.csr",FUN = "missing",make.dimnames = "missing"), .kron.matrix.csr)
	setMethod("%x%",signature(X="matrix.csr",Y="numeric",FUN = "missing",make.dimnames = "missing"), .kron.matrix.csr)
	setMethod("%x%",signature(X="numeric",Y="matrix.csr",FUN = "missing",make.dimnames = "missing"), .kron.matrix.csr)
	setMethod("%x%",signature(X="matrix",Y="matrix.csr",FUN = "missing",make.dimnames = "missing"), .kron.matrix.csr)
	setMethod("%x%",signature(X="matrix.csr",Y="matrix",FUN = "missing",make.dimnames = "missing"), .kron.matrix.csr)
        }

setGeneric("image", function(x, ...) standardGeneric("image")) 
	# /*RSB*/  changed definition

setMethod("image","matrix.csr",
function(x,col=c("white","gray"),xlab="column",ylab="row", ...){
	n <- x@dimension[1]
	p <- x@dimension[2]
	z <- matrix(0,n,p)
	column <- x@ja
	row <- rep(n:1,diff(x@ia))
	z[cbind(row,column)] <- 1
	image.default(x=1:p,y=-(n:1),t(z),axes=FALSE, col=col,xlab=xlab,ylab=ylab, ...) # /*RSB*/ changed call to make sure
	axis(1,pretty(1:p), ...)
	axis(2,pretty(-(n:1)),labels=rev(pretty(1:n)), ...)
	box()
	})
#setMethod("summary","slm",summary.slm)
#setMethod("summary","mslm",summary.mslm)
#setMethod("coef","slm",coef.slm)
#setMethod("fitted","slm",fitted.slm)
#setMethod("residuals","slm",residuals.slm)
#setMethod("print","summary.slm",print.summary.slm)
#--------------------------------------------------------------------
