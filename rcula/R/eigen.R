
setClass("geigen",
		representation(
				vectors = "gmatrix",
				values  = "gvector"),
		prototype = list(
				vectors = NULL,
				values  = NULL)
) 

setGeneric("vectors",
		function(a)
			standardGeneric("vectors")
)

setMethod("vectors", c("geigen"),
		function(a) {
			return(a@vectors)
		})

setGeneric("values",
		function(a)
			standardGeneric("values")
)

setMethod("values", c("geigen"),
		function(a) {
			return(a@values)
		})

gsymeigen = function (x) 
{
	x <- as.gmatrix(x,dup=FALSE)
	if(x@type>1L)
		type(x)=0L
	x =gnamestrip(x)
	n <- nrow(x)
	if (!n) 
		stop("0 x 0 matrix")
	if (n != ncol(x)) 
		stop("non-square matrix in 'eigen'")
	
	if (sum(!is.finite(x),retgpu=FALSE)>0) 
		stop("infinite or missing values in 'x'")
	
	checkDevice(x@device)
	
	ret=new("geigen", vectors=gdup(x), values=gvector(n))
	z <- .Call("rcula_eigen_symm", ret@vectors, ret@values)

	ord <- order(h(ret@values), decreasing = TRUE)
	ret@values =ret@values[ord]
	ret@vectors=ret@vectors[,ord, drop=FALSE]
	return(ret)
}
