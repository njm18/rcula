#
# GPUMatrix Makefile
# --------------
# by Nathan Morris
#          (Modified from the R/GPU Makefile by M. Kempenaar)
#


#to debug... make sure nvcc_flags includes -g -G and shlib flag includes -d
#            R -d cuda-gdb
#            run --vanilla
#            Run some code to load in the .so (dyn.load)
#            ctr + c to get back to debugger
#            b local.cu:[line#]
#            c to continue
#            s to step
#            d b to delet al breakpoints


# Nvidia CUDA Compiler
NVCC = nvcc #/usr/local/cuda/bin/nvcc


# LIB / INC locations (default, change these using 'export' in the terminal, please refer
# to the manual on the rcula install instructions)
R_INC_PATH    ?= /usr/include/R
CUDA_INC_PATH ?= /usr/local/cuda/include 
CULA_INC_PATH ?= /usr/local/cula/include 
CUDA_LIB_PATH ?= /usr/local/cuda/lib64
CULA_LIB_PATH ?= /usr/local/cula/lib64

CULA_ARCH ?= sm_20 

# R gcc build script
R_CC = R CMD SHLIB


# Output Files:
#   CUDA output
OUT_OBJ =  rcula.o #solve.o eigen.o lu.o
#   R output
OUT_SO = rcula.so

# CUDA Libraries
LIBS = -lcudart -lcublas -lcula_core -lcula_lapack -d -lRlapack

# CUDA Flags (-g -G removed)
NVCC_FLAGS = --shared -Xcompiler "-fPIC"  -arch=$(CULA_ARCH)



# Build rules
build: $(OUT_OBJ)
	$(R_CC) $(OUT_OBJ) -o $(OUT_SO) $(LIBS) -L$(CUDA_LIB_PATH) -L$(CULA_LIB_PATH)
	
	
rcula.o:
	$(NVCC) -c rcula.cu -o rcula.o $(NVCC_FLAGS) -I$(CUDA_INC_PATH) -I$(CULA_INC_PATH) -I$(R_INC_PATH)



clean:
	-rm -f *.o core
	-rm -f *.so core
	-rm -f *.out
	-rm -f *.linkinfo

rebuild: clean build

