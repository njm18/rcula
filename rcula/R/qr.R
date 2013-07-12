# TODO: Add comment
# 
# Author: nmorris
###############################################################################


setClass("gqr",
		representation(
				qr="gmatrix",
				qraux = "gvector"),
		prototype = list(
				qr=NULL,
				qraux=NULL)
) 


setMethod("vectors", c("geigen"),
		function(a) {
			return(a@vectors)
		})

setMethod("qr", "gmatrix",
		function(x) 
		{	
			checkDevice(x@device)
			if(x@type>1L)
				type(x)=-0L
			checkDevice(x@device)
			res=new("gqr",
					qr=gdup(x),
					qraux=gvector(ncol(x), type=x@type)
			)
			#browser()
			tmp <- .Call("rcula_qr", res@qr, res@qraux@ptr)
			#res@rank=as.integer(sum(res@qraux/max(res@qraux)>10^-6, retgpu=FALSE))
			#if (!is.null(cn <- colnames(x))) 
			#	colnames(res@qr) <- cn[res@pivot]
			if (!is.null(cn <- colnames(x))) 
				colnames(res@qr) <- cn
			return(res)
		}
)



setMethod("qr.coef", "gqr",
		function (qr, y) 
		{
			#browser()
			if(qr@qr@type>1L)
				type(qr@qr)=0L
			if(qr@qraux@type>1L)
				type(qr@qraux)=0L
			if(qr@qr@type!=qr@qraux@type) {
				totype=min(c(qr@qr@type,qr@qraux@type))
				type(qr@qr)=totype
				type(qr@qraux)=totype
			}

			if (class(y)!="gmatrix") 
				y <- as.gmatrix(y)
			else
				y=gdup(y)
			
			if(y@type>1L)
				type(y)=type(qr@qr)
			if(y@type!=qr@qr@type)
				stop("Type mismatch.")
			if(qr@qraux@type!=qr@qr@type)
				stop("Type mismatch.")
			
			checkDevice(c(y@device,qr@qr@device,qr@qraux@device))
			
			n <- as.integer(nrow(qr@qr))
			p <- as.integer(ncol(qr@qr))
			#k <- as.integer(qr$rank)
			ny <- as.integer(ncol(y))
			
			if (p == 0L) 
				error("gmatrix in qr has a dimension of 0")
			#browser()
			#ix <- if (p > n) 
			#			c(seq_len(n), rep(NA, p - n))
			#		else seq_len(p)
			#coef <- gmatrix(NA, nrow = p, ncol = ny)
			#coef[qr$pivot, ] <- .Call("qr_coef_real", qr@qr, qr@qraux@ptr, y)[ix,]
			dummy=.Call("rcula_modqr_coef", qr@qr, qr@qraux, y)
			if (p > n) 
				stop("Columns of qr matrix must be less than or equal to rows.")
			return(y)
		}
)


