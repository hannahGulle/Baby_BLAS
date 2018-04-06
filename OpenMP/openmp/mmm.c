#ifdef __cplusplus
extern "C" {
#endif
    void mmm_( int *threads, int *len,  double *a, double *b, double*c );
#ifdef __cplusplus
    }
#endif


#include <omp.h>
#include <stdio.h>

void mmm_( int *threads, int *len,  double *a, double *b, double *c ){

    int i, j, k;
    int veclen = *len;
    omp_set_num_threads(*threads);

#pragma omp parallel shared(veclen) private(i,j,k)
{
    printf("Thread %d started\n", omp_get_thread_num() );

    for (i=0; i<veclen; i++) {
        for (j=0; j<veclen; j++) {
            *(c+(i*veclen+j)) = 0.0;
            for (k=0;k<veclen;k++){
                *(c+(i*veclen+j)) += *(a+(i*veclen+k)) * *(b+(k*veclen+j)); 
            }
        }
    }

    PAPI_unregister_thread();
    printf("Thread %d finished\n", omp_get_thread_num() );
}


}


