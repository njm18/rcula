The "rcula" Package
=================================================
NOTE: This package is still in the early stages of development.
The rcula package is provides an R interface for the [CULA](http://www.culatools.com/) dense linear algebra library. 
The actual CULA library is not in general free, but an academic liscence could be obtained freely at the time this package was released. This package is meant to extend the abilities of the package [gmatrix](https://github.com/njm18/gmatrix/tree/master/gmatrix). 

Installation Instructions
-------------------------
1. Install the the prequasites: [CUDA Toolkit](https://developer.nvidia.com/cuda-downloads),
 [R](http://cran.r-project.org/) and [gmatrix](https://github.com/njm18/gmatrix/tree/master/gmatrix). 
2. Optain the [CULA](http://www.culatools.com/) dense library and install it.
3. Start R and then install the 'rcula' package with the following commands. Package compilation may take some time.

```
download.file("http://solomon.case.edu/rcula/rcula_0.5.tar.gz", "rcula.tar.gz")
install.packages("rcula.tar.gz", repos = NULL)
file.remove("rcula.tar.gz")
```
	 
Installation Note
-----------------
By default, when compiling, the makefile assumes that
+ The CUDA library files are located in the folder /usr/local/cuda/lib64.
+ The CUDA include path is /usr/local/cuda/include
+ The CULA library files are located in /usr/local/cula/lib64
+ The CULA include path is /usr/local/cula/include
+ The R libraries are located in the folder /usr/include/R.
+ The compute capibility of the target device is 2.0.

If these are incorrect assumptions, the user may set these values and install using the follwing R commands as an example.
First set the environmental variables:

    Sys.setenv(CUDA_INC_PATH="/usr/include/cuda-5.0/include")
    Sys.setenv(CULA_INC_PATH="/usr/include/cula_new/include")
    Sys.setenv(CUDA_LIB_PATH="/usr/include/cuda-5.0/lib64")
    Sys.setenv(CULA_LIB_PATH="/usr/include/cula_new/lib64")
    Sys.setenv(R_INC_PATH="/usr/local/R/R-2.15.0/lib64/R/include")
    Sys.setenv(NVCC_ARCH="-gencode arch=compute_30,code=sm_30")
    
Next install the package as above:

    download.file("http://solomon.case.edu/rcula/rcula_0.5.tar.gz", "rcula.tar.gz")
    install.packages("rcula.tar.gz", repos = NULL)
    file.remove("rcula.tar.gz")
	    
Testing the Installation
-------------------------
We recoment that the user test the installation using the following commands:

    library(rcula)
    culaTest()
    
Please report any errors to the package maintainer.

Getting Started
---------------
Here are some examples of how to use the solve function:

    A=gmatrix(grnorm(400),20,20)
    b=grnorm(20)
    Ainverse = solve(A)
    z=solve(A,b)

Here are some examples of how to use the qr function:

    B=qr(A)
    Ainverse2=solve(B)

Examples of symmetrice eigen decomposition:

    C=crossprod(A)
    Ceig=gsymeigen(C)
    vals = h(values(Ceig))
    vecs = h(vectors(Ceig))