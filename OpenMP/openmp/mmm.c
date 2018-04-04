#ifdef __cplusplus
extern "C" {
#endif
    void mmm_( int *threads, int *len,  double *a, double *b, double*c );
#ifdef __cplusplus
    }
#endif

/*  O p e n M P     C O D E  */

#include <omp.h>

void mmm_( int *threads, int *len,  double *a, double *b, double *c ){

    int i, j, k;
    int veclen = *len;

// Set the number of threads to use here

    omp_set_num_threads(*threads);

#pragma omp parallel shared(veclen) private(i,j,k)
{
    for (i=0; i<veclen; i++) {
        for (j=0; j<veclen; j++) {
            *(c+(i*veclen+j)) = 0.0;
            for (k=0;k<veclen;k++){
                *(c+(i*veclen+j)) += *(a+(i*veclen+k)) * *(b+(k*veclen+j)); 
            }
        }
    }
}
}


