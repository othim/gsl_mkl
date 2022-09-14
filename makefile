CC = icpc
OPT = -O3

BLAS_LIB = MKL


LDFLAGS = -L/net/home/toliver/gsl-2.7/lib -L/usr/local/lib -fopenmp
ifeq ($(BLAS_LIB), MKL)

#LDFLAGS  += -L/net/opt/intel/2022.1.2.146/intel/oneapi/mkl/2022.0.2/lib/intel64 \
			-Wl,--no-as-needed -lmkl_intel_ilp64 -lmkl_gnu_thread -lmkl_core -lgomp -lpthread -ldl
LDFLAGS += -qmkl 
endif
ifeq ($(BLAS_LIB), GSL)
LDFLAGS += -lgslcblas
endif

CFLAGS = -Wall -fPIC $(OPT) \
		 -I/net/home/toliver/gsl-2.7/gsl

CFLAGS2 = -Wall -std=c++11 $(OPT)
main: main.o
	$(CC) $(CFLAGS2) -o prog $(LDFLAGS) main.o -lgsl 
	rm main.o

main.o: main.cpp
	$(CC) $(CFLAGS) -c main.cpp
clean:
	rm prog main.o
