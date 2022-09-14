#include <gsl/gsl_matrix.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_cblas.h>
#include <gsl/gsl_blas.h>
#include <iostream>
#include <ctime>
#include <complex>
#include <string>
#include <vector>

void test_dgemv();
void test_zgemm();
void test_onshell();

int main (int argc, char** argv) {
    std::vector<std::string> options;
    options.push_back("dgemv");
    options.push_back("zgemm");
    options.push_back("onshell");
    if (argc < 2)
    {
        std::cout << "No option. Options:" << std::endl;
        for( auto i : options)
        {
            std::cout << i << std::endl;
        }
        std::cout << "\n ending program" << std::endl;
        return 1;
    }

    std::cout << argv[1] << std::endl;
    if (std::string(argv[1]) == "dgemv")
    {
        test_dgemv();
    } else if (std::string(argv[1]) == "zgemm")
    {
        test_zgemm();
    } else if (std::string(argv[1]) == "onshell")
    {
        test_onshell();
    }
    return 0;
}

void print_time(std::clock_t start,std::clock_t end,int N)
{
    std::cout <<  "Time (ms): " << 1e3*(double)(end-start)/((double)N*(double)CLOCKS_PER_SEC)
          << std::endl; 
    std::cout <<  "Time (us): " << 1e6*(double)(end-start)/((double)N*(double)CLOCKS_PER_SEC)
         << std::endl; 
}

void test_dgemv()
{
    std::cout << "test_dgemv():\n----------------" << std::endl;
    int M = 5000;
    std::cout << "M=" << M << std::endl;
    gsl_vector* v = gsl_vector_calloc(M);
    gsl_matrix* m = gsl_matrix_calloc(M, M);
    

    std::clock_t start,end;
    int N = 300;
    // Compiling with MKL instead of GSL gives a factor 2 speedup for 
    // M = 5000
    start = std::clock();
    for (int i=0; i < N; i++) {
        double a = (double)i*i;
        gsl_blas_dgemv(CblasNoTrans, a, m, v, 0.0, v);
    }
    end = std::clock();
    print_time(start,end,N);
    // Test of complex types
    // ---------------------
    gsl_complex* gsl_z = &gsl_complex_rect(1,2);
    std::cout << GSL_REAL(*gsl_z) << ", " << GSL_IMAG(*gsl_z) << std::endl;
    std::complex<double>* b = (std::complex<double>*)gsl_z;
    std::cout << *b << std::endl;
    // ---------------------
    
    std::cout << "----------------\ntest_dgemv ended" << std::endl;
}

void set_random_matrix(gsl_matrix_complex* M)
{
    for(int i=0; i<M->size1; i++)
    {
        for(int j=0; j<M->size2; j++)
        {
            gsl_complex el = gsl_complex_rect((double)std::rand()/(double)RAND_MAX
                    ,(double)std::rand()/(double)RAND_MAX);
            //std::cout << GSL_REAL(el) << "," << GSL_IMAG(el) << 
            //    std::endl;
            gsl_matrix_complex_set(M,i,j,el);
        }
    }
}

void test_zgemm()
{
    std::cout << "test_zgemm():\n----------------" << std::endl;
    int M = 200;
    std::cout << "M=" << M << std::endl;
    gsl_matrix_complex* m1 = gsl_matrix_complex_calloc(M, M);
    gsl_matrix_complex* m2 = gsl_matrix_complex_calloc(M, M);
    gsl_matrix_complex* res = gsl_matrix_complex_calloc(M, M);
    
    set_random_matrix(m1);
    set_random_matrix(m2);

    std::clock_t start,end;
    int N = 300;
    for(int i=0; i<N; i++)
    {
        gsl_complex alpha = gsl_complex_rect(1,0);
        gsl_complex beta = gsl_complex_rect(0,0);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, m1,m2,beta,res);
    }
    end = std::clock();
    print_time(start,end,N);
    
    std::cout << "----------------\ntest_dgemv ended" << std::endl;
}


void on_shell_mult(gsl_matrix_complex* m1, gsl_matrix_complex* m2, 
        gsl_matrix_complex* res)
{
    /*
     * This function implements the onshell multiplication, meaning 
     * that the matrix multiplixation does not include the last 
     * column/row of the matrix. This is achieved by first multiplying
     * the matrices as usual and then subtracting the error that is 
     * induced. The time lost is negligable compared to just doing
     * an ordinary multiplication.
     *
     * m1, m2, res needs to be distinct matrices. You cannot have eg.
     * m1 <- m1*m2
     *
     * This function is tested againts on_shell_mult_bf() which is a loop
     * brute force verion of the original sum that we want to compute.
     */

    // Multiply as usual
    gsl_complex alpha = gsl_complex_rect(1,0);
    gsl_complex beta = gsl_complex_rect(0,0);
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, alpha, m1,m2,beta,res);

    // Correct for the inclusion of the on-shell elements in the 
    // matrix-matrix product
    
    // Take the last column of m1
    gsl_vector_complex_view m1_vv = gsl_matrix_complex_column(m1,m1->size2-1);
    gsl_vector_complex* v1 = &m1_vv.vector;
    // Take the last row of m2
    gsl_vector_complex_view m2_vv = gsl_matrix_complex_row(m2,m2->size1-1);
    gsl_vector_complex* v2 = &m2_vv.vector;
    
    // Perform the outer product
    gsl_complex minus_one = gsl_complex_rect(-1,0);
    gsl_blas_zgeru(minus_one, v1,v2, res);
}

void on_shell_mult_bf(gsl_matrix_complex* m1, gsl_matrix_complex* m2, 
        gsl_matrix_complex* res)
{
    for(int i=0; i<m1->size1; i++)
    {
        for(int j=0; j<m2->size2; j++)
        {
            gsl_complex el = gsl_complex_rect(0,0);
            for(int k=0; k< m1->size1-1; k++)
            {
                el = gsl_complex_add(el, gsl_complex_mul(gsl_matrix_complex_get(m1,i,k),
                            gsl_matrix_complex_get(m2,k,j)));
            }
            gsl_matrix_complex_set(res,i,j,el);
        }
    }
}

void test_onshell()
{
    std::cout << "test_onshell():\n----------------" << std::endl;
    int M = 100;
    std::cout << "M=" << M << std::endl;
    gsl_matrix_complex* m1 = gsl_matrix_complex_calloc(M, M);
    gsl_matrix_complex* m2 = gsl_matrix_complex_calloc(M, M);
    gsl_matrix_complex* res1 = gsl_matrix_complex_calloc(M, M);
    gsl_matrix_complex* res2 = gsl_matrix_complex_calloc(M, M);
    
    set_random_matrix(m1);
    set_random_matrix(m2);

    std::clock_t start,end;
    int N = 300;
    start = std::clock();
    for(int i=0; i<N; i++)
    {
        on_shell_mult(m1,m2,res1);
    }
    end = std::clock();
    print_time(start,end,N);
    
    start = std::clock();
    for(int i=0; i<N; i++)
    {
        on_shell_mult_bf(m1,m2,res2);
    }
    end = std::clock();
    print_time(start,end,N);

    double max_abs_diff = 0.0;
    for(int i=0; i<M; i++)
    {
        for(int j=0; j<M; j++)
        {
            gsl_complex r1 = gsl_matrix_complex_get(res1,i,j);
            gsl_complex r2 = gsl_matrix_complex_get(res2,i,j);
            gsl_complex diff = gsl_complex_sub(r1,r2);
            /*std::cout << GSL_REAL(r1) << "," << GSL_IMAG(r1) << 
                std::endl;
            std::cout << GSL_REAL(r2) << "," << GSL_IMAG(r2) << 
                std::endl;
            std::cout << GSL_REAL(diff) << "," << GSL_IMAG(diff) << 
                std::endl << std::endl;*/
            max_abs_diff = std::max(max_abs_diff, std::sqrt(GSL_REAL(diff)*GSL_REAL(diff)+
                GSL_IMAG(diff)*GSL_IMAG(diff)));
        }
    }
    std::cout << "Max abs diff: " << max_abs_diff << std::endl;
    std::cout << "----------------\ntest_onshell ended" << std::endl;
}
