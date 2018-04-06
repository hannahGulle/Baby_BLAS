#include<stdio.h>
#include<omp.h>

#ifdef __cplusplus
extern "C" {
#endif
    void mvv_(  int *num_threads, int *N, double *ma, double *va, double *vr);
#ifdef __cplusplus
    }
#endif

void mvv_( int *num_threads, int *N, double *ma, double *va, double *vr ){

	int n = *N;
	int i, j;

	omp_set_num_threads(*num_threads);

#pragma omp parallel shared(n) private(i,j) 
{
	printf("Thread %d started\n", omp_get_thread_num() );

	#pragma omp for
	for ( j = 1; j < n; j++ ){
		for( i = 1; i < n; i++ ){
			vr[j] += *(ma+i+j*n) * *(va+i);
		}
	}


	PAPI_unregister_thread();
	printf("Thread %d finished\n", omp_get_thread_num() );
}
}
