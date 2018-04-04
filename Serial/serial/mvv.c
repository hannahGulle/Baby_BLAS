#include<stdio.h>
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
	for ( j = 1; j < n+1 ; j++ ){
		for( i = 1; i < n+1; i++ ){
			
			*(vr+j) += *(ma+i+j*n) * *(va+i);
		}
		printf("value @ %i is %d\n", j, *(vr+j));
	}
}
