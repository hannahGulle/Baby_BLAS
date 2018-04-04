
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
    void mmm_( int *threads, int *len,  double *a, double *b, double*c );
#ifdef __cplusplus
}
#endif

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

/* Function prototypes */

void *mmm_thread_worker();

/* Struct for passing data to thread functions */

struct args {
    int N; 
    int startRow;
    int stopRow;
    double *Aptr;
    double *Bptr;
    double *Cptr;
};

void mmm_( int *threads, int *len, double *A, double *B, double *C ){

    // This function has to break up the data, spawn the processes, gather the results, and 
    // clean up after itself.


    int numThreads = *threads;
    int matrixDimension = *len;
    int *numberOfRows;
    int startRow, stopRow;
    pthread_t *thread_id;
    struct args *thread_args;


    // if there are fewer dimensions than threads, do a simple single-threaded matrix multiplication.
    if ( matrixDimension < numThreads ) {
        for (int i=0;i<matrixDimension;i++) {
            for (int j=0;j<matrixDimension;j++) {
                *(C+(i*matrixDimension+j))=0.0;
                for (int  k=0;k<matrixDimension;k++) {
                    *(C+(i*matrixDimension+j)) += *(A+(i*matrixDimension+k)) * *(B+(k*matrixDimension+j));
                }
            }  
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
                thread_args->Aptr = A;
                thread_args->Bptr = B;
                thread_args->Cptr = C;

                pthread_create( thread_id+i, NULL, &mmm_thread_worker, thread_args );
            }
        }
        for(int i=0; i < numThreads ; i++) {
            pthread_join( *(thread_id+i), NULL); 
        }

        free(numberOfRows);
        free(thread_id);
    }

    // END OF MMM -- the memory pointed to by *C should now contain the product A.B
}


void *mmm_thread_worker( struct args *thread_args  ) {

    // This is where the actual work of multiplying each matrix takes place --
    // with each mmm_thread_worker only responsible for a piece of the matrix A.
    
    // It is worth noting here that while the thread workers are reading simultaneously
    // from *A and *B, they never write to same memory locations in *C.

    int i, j, k;
    double val;
    int rowStart, rowStop, N; 
    double *A, *B, *C;

    // Unpack the thread_args struct into normal variables
    N        =  thread_args->N;
    rowStart =  thread_args->startRow;
    rowStop  =  thread_args->stopRow; 
    A        =  thread_args->Aptr;
    B        =  thread_args->Bptr;
    C        =  thread_args->Cptr;

    // Process the rows for which this thread is responsible
    for (i=rowStart;i<rowStop;i++) {
        for (j=0;j<N;j++) {
            *(C+(i*N+j))=0.0;
            for (k=0;k<N;k++) {
                *(C+(i*N+j)) += *(A+(i*N+k)) * *(B+(k*N+j));
            }
        }  
    } 

    free(thread_args);
    pthread_exit(NULL);
}




