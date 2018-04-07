
/**************************************************************
 *
 * A matrix-matrix multiplication function designed
 * to work optimally over a given number of threads using
 * the POSIX Threads library (pthreads)
 *
 * Prof. Andrew J. Pounds, Ph.D.
 * Mercer University
 * Spring 2018
 *
 **************************************************************/

#ifdef __cplusplus
extern "C" {
#endif
    double dot_( int *threads, int *len,  double *va, double *vb );
#ifdef __cplusplus
}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

/* Function prototypes */

void *dot_thread_worker();

/* Struct for passing data to thread functions */

struct args {
    int N; 
    int startRow;
    int stopRow;
    double *VAptr;
    double *VBptr;
    double dot;
};

double dot_( int *threads, int *len, double *va, double *vb){

    // This function has to break up the data, spawn the processes, gather the results, and 
    // clean up after itself.


    int numThreads = *threads;
    int matrixDimension = *len;
    int *numberOfRows;
    double *dotProds;
    int startRow, stopRow;
    pthread_t *thread_id;
    struct args *thread_args;
    double sum = 0.0;

    // if there are fewer dimensions than threads, do a simple single-threaded matrix multiplication.
    if ( matrixDimension < numThreads ) {
        for (int i=0;i<matrixDimension;i++) {
        	sum += va[i] * vb[i];
	} 
    }

    else { 

        /* 
         *  The parallel work begins.  The process is to first determine how to break
         *  up the matrix work equitably across the threads.  Once this is done the struct is filled with
         *  the information and a thread is started using the information.  Other than the size of the
         *  matrices and the rows to be processed, only the pointers to the memory locations of the matrices
         *  are passed.
         */ 

        // Malloc an array to keep up with thread id's for each thread
        thread_id = (pthread_t *) malloc (numThreads * sizeof(pthread_t));

        // Malloc an array to keep up with how many rows to work on in each thread
        numberOfRows = ( int * ) malloc( numThreads * sizeof(int) );

	// Malloc an array to keep up with the sums of each thread
	dotProds = ( double * ) malloc( numThreads * sizeof(double) );

        // Here we detemine the number of rows over which each thread will work
        for (int i=0; i<numThreads; i++ ){
            *(numberOfRows+i) = matrixDimension / numThreads;
        }
        for (int i=0; i< matrixDimension % numThreads; i++ ){
            *(numberOfRows+i) = *(numberOfRows+i) + 1;
        }

        // Now that we know how many rows each thread will be responsible for computing,
        // malloc memory for the struct data, pack the struct with the thread-specific info
        // on where it is to start and stop processing, and create the thread using this data. 

        stopRow=0;
        for(int i=0; i < numThreads ; i++) {
            {   
                startRow=stopRow;
                stopRow=startRow+*(numberOfRows+i);
                thread_args = ( struct args * )  malloc(sizeof( struct args));
                thread_args->N   = matrixDimension;
                thread_args->startRow = startRow;
                thread_args->stopRow = stopRow; 
                thread_args->VAptr = va;
                thread_args->VBptr = vb;
		thread_args->dot = *(dotProds+i);

                pthread_create( thread_id+i, NULL, &dot_thread_worker, thread_args );
            }
        }
        for(int i=0; i < numThreads ; i++) {
            pthread_join( *(thread_id+i), NULL); 
        }

        free(numberOfRows);
        free(thread_id);
   	
	for( int i = 0; i < numThreads; i++ ){
		sum += dotProds[i];
		printf("Sum at Thread %i\t is %d\n", i, sum);
	}
    }
    return sum;
    // END OF MMM -- the memory pointed to by *C should now contain the product A.B
}


void *dot_thread_worker( struct args *thread_args  ) {

    // This is where the actual work of multiplying each matrix takes place --
    // with each mmm_thread_worker only responsible for a piece of the matrix A.
    
    // It is worth noting here that while the thread workers are reading simultaneously
    // from *A and *B, they never write to same memory locations in *C.

    int i;
    int rowStart, rowStop; 
    double *va, *vb;
    double sum = 0.0;

    // Unpack the thread_args struct into normal variables
    rowStart =  thread_args->startRow;
    rowStop  =  thread_args->stopRow; 
    va       =  thread_args->VAptr;
    vb       =  thread_args->VBptr;

    // Process the rows for which this thread is responsible
    for (i=rowStart;i<rowStop;i++) {
    	sum += va[i] * vb[i];
    } 
    thread_args->dot = sum;

    free(thread_args);
    pthread_exit(NULL);
}




