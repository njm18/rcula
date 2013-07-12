# TODO: Add comment
# 
# Author: nmorris
###############################################################################


.onLoad <- function(lib, pkg) {
	.Call("rcula_initialize")
}
.onUnload <- function(libname) {
	.Call("rcula_shutdown")
}
