#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <iostream>
#include <ctime>


int main () {
    gsl_vector* v = gsl_vector_calloc(5000);
    gsl_matrix* m = gsl_matrix_calloc(5000, 5000);
    
    std::clock_t start,end;
    int N = 300;
    start = std::clock();
    for (int i=0; i < N; i++) {
        double a = (double)i*i;
        gsl_blas_dgemv(CblasNoTrans, a, m, v, 0.0, v);
    }
    end = std::clock();
    std::cout <<  "Time: " << 1e6*(double)(end-start)/((double)N*(double)CLOCKS_PER_SEC)
          << " us" << std::endl; 
    std::cout << "Program ended" << std::endl;
    return 0;
}
