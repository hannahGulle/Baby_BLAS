#include<stdio.h>
#include<omp.h>

#ifdef __cplusplus
extern "C" {
#endif
    void vvm_( int *num_threads, int *len, double *va, double *vb, double *ma);
#ifdef __cplusplus
    }
#endif

// Computes the tensor product of two vectors 

void  vvm_( int *num_threads, int *len, double *va, double *vb, double *ma){

	int i, j;
	int threads = *num_threads;
	int alength = *len;
	omp_set_num_threads(*num_threads);
	
#pragma omp parallel shared(alength) private(i,j)
{
	printf("Thread %d started\n", omp_get_thread_num() );

	for (i=0; i<alength; i++) {
		for (j=0; j<alength; j++) {
			*(ma+(alength*i)+j) = *(va+i) * *(vb+j);
		}
	}

	PAPI_unregister_thread();
	printf("Thread %d finished\n", omp_get_thread_num() );
}
}
