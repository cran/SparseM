require(methods)
setClass("matrix.csr",representation(ra="numeric",
	ja="integer",ia="integer", dimension="integer"),
	validity = function(object) {
 		if(!length(object@dimension) == 2 )
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
			return("ia has wrong number of elments")
        	if(length(object@ra) < 1 || length(object@ra) > nrow*ncol)
                	return("ra has too few, or too many elements")
		TRUE})
setMethod("initialize", "matrix.csr", 
	function(.Object, ra = numeric(0),  ja = integer(0),
		ia = integer(0),dimension = integer(0)) {
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
        	if(length(object@ra) < 1 || length(object@ra) > nrow*ncol)
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
	xsuper="numeric",ierr="numeric",time="numeric"))
setClass("numeric or NULL")
setIs("numeric","numeric or NULL")
setIs("NULL","numeric or NULL")
setClass("character or NULL")
setIs("character","character or NULL")
setIs("NULL","character or NULL")
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
setMethod("as.matrix.csr","matrix.csc",as.matrix.csr.matrix.csc)
setMethod("as.matrix.csr","matrix.ssr",as.matrix.csr.matrix.ssr)
setMethod("as.matrix.csr","matrix.ssc",as.matrix.csr.matrix.ssc)
setMethod("as.matrix.csr","matrix.coo",as.matrix.csr.matrix.coo)
setGeneric("as.matrix.csc")
setMethod("as.matrix.csc","matrix.csr",as.matrix.csc.matrix.csr)
setMethod("as.matrix.csc","matrix.ssr",as.matrix.csc.matrix.ssr)
setMethod("as.matrix.csc","matrix.ssc",as.matrix.csc.matrix.ssc)
setGeneric("as.matrix.ssr")
setMethod("as.matrix.ssr","matrix.csr",as.matrix.ssr.matrix.csr)
setMethod("as.matrix.ssr","matrix.csc",as.matrix.ssr.matrix.csc)
setMethod("as.matrix.ssr","matrix.ssc",as.matrix.ssr.matrix.ssc)
setGeneric("as.matrix.ssc")
setMethod("as.matrix.ssc","matrix.csr",as.matrix.ssc.matrix.csr)
setMethod("as.matrix.ssc","matrix.csc",as.matrix.ssc.matrix.csc)
setMethod("as.matrix.ssc","matrix.ssr",as.matrix.ssc.matrix.ssr)
setGeneric("as.matrix.coo")
setMethod("as.matrix.coo","matrix.csr",as.matrix.coo.matrix.csr)
setMethod("as.matrix","matrix.csr",as.matrix.matrix.csr)
setMethod("as.matrix","matrix.csc",as.matrix.matrix.csc)
setMethod("as.matrix","matrix.ssr",as.matrix.matrix.ssr)
setMethod("as.matrix","matrix.ssc",as.matrix.matrix.ssc)
setMethod("as.matrix","matrix.coo",as.matrix.matrix.coo)
setMethod("t","matrix.csr",t.matrix.csr)
setMethod("t","matrix.csc",t.matrix.csc)
setMethod("t","matrix.coo",function(x) as.matrix.coo(t(as.matrix.csr(x))))
setMethod("dim","matrix.csr",dim.matrix.csr)
setMethod("dim","matrix.csc",dim.matrix.csc)
setMethod("dim","matrix.ssr",dim.matrix.ssr)
setMethod("dim","matrix.ssc",dim.matrix.ssc)
setMethod("dim","matrix.coo",dim.matrix.coo)
setMethod("diff","matrix.csr",diff.matrix.csr)
setGeneric("diag")
setMethod("diag","matrix.csr",diag.matrix.csr)
setMethod("diag<-","matrix.csr",diag.assign.matrix.csr)
setGeneric("chol")
setMethod("chol","matrix",chol.default)
setMethod("chol","matrix.csr",chol.matrix.csr)
setMethod("chol","matrix.csc",chol.matrix.csc)
setMethod("backsolve","matrix.csr.chol",backsolve.matrix.csr.chol)
setMethod("solve","matrix.csr",solve.matrix.csr)
setMethod("model.matrix","matrix.csc.hb",model.matrix.matrix.csc.hb)
setMethod("model.matrix","matrix.ssc.hb",model.matrix.matrix.ssc.hb)
setMethod("model.response",signature(data="matrix.csc.hb",type="ANY"),model.response.matrix.csc.hb)
setMethod("model.response",signature(data="matrix.ssc.hb",type="ANY"),model.response.matrix.csc.hb)
setMethod("%*%",signature(x="matrix.csr",y="matrix.csr"),.matmul.matrix.csr)
setMethod("%*%",signature(x="matrix.csr",y="matrix"),.matmul.matrix.csr)
setMethod("%*%",signature(x="matrix.csr",y="numeric"),.matmul.matrix.csr)
setMethod("%*%",signature(x="matrix",y="matrix.csr"),.matmul.matrix.csr)
setMethod("%*%",signature(x="numeric",y="matrix.csr"),.matmul.matrix.csr)
setMethod("image","matrix.csr",image.matrix.csr)
setMethod("summary","slm",summary.slm)
setMethod("summary","mslm",summary.mslm)
setMethod("coef","slm",coef.slm)
setMethod("fitted","slm",fitted.slm)
setMethod("residuals","slm",residuals.slm)
setMethod("print","summary.slm",print.summary.slm)
#--------------------------------------------------------------------
