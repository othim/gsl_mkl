CC = g++
OPT = -O3

BLAS_LIB = MKL

CFLAGS2 = -Wall -std=c++11 -stdlib=libc++ $(OPT)
LDFLAGS =  -undefined dynamic_lookup -L/Users/toliver/wigxjpf-1.11/lib/  -L/usr/local/Cellar/gcc/11.3.0/lib/gcc/11/  \
		   -L/usr/local/opt/llvm/lib 

ifeq ($(BLAS_LIB), MKL)
LDFLAGS += -L/opt/intel/oneapi/mkl/2021.4.0/lib -lmkl_intel_ilp64 -lmkl_sequential -lmkl_core
endif
ifeq ($(BLAS_LIB), GSL)
LDFLAGS += -lgslcblas
endif

CFLAGS = -Wall -fPIC -std=c++11 -stdlib=libc++ $(OPT) \
		 -I/opt/anaconda3/envs/nn-mwpc-env/include \
		 -I/opt/anaconda3/include/python3.8 -I/Users/toliver/wigxjpf-1.11/inc/ \
		 -I/usr/local/Cellar/gsl/2.7/include/gsl -I/usr/local/opt/llvm/include


main: main.cpp
	$(CC) $(CFLAGS2) -o proj main.cpp $(LDFLAGS) -lm -lgsl

clean:
	rm proj
