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
		s <- ((x@ja %in% cl) & (x@ia %in% rw))
		ja <- match(x@ja[s],cl)
		ia <- match(x@ia[s],rw)
		dim <- c(length(rw),length(cl))
		ra <- x@ra[s]
		if (length(ra) == 0 && length(ja) == 0){ #trap all zeros returned matrix
			ra = 0
			ja = ia = as.integer(1)
			}
		A <- new("matrix.coo",ra=ra,ja=ja,ia=ia,dimension=dim)
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
