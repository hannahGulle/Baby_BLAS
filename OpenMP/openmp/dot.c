#include<stdio.h>
#include<omp.h>

#ifdef __cplusplus
extern "C" {
#endif
    double dot_(  int *num_threads, int *N, double *va, double *vb);
#ifdef __cplusplus
    }
#endif

double dot_( int *num_threads, int *N, double *va, double *vb ) {

	int i;
	double sum = 0.0;
	int n = *N;
	int nthreads = *num_threads;

#pragma omp parallel private(i) shared(n)
{
	#pragma omp for reduction(+:sum)
	for ( i = 1; i < n+1; i++ ) {
		sum += va[i] * vb[i];
	}
}
	return sum;
}
