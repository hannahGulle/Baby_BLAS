
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
    void ils_( int *threads, int *len,  double *a, double *b, double*x );
#ifdef __cplusplus
}
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

/* Function prototypes */

void dls_( int *threads, int *len, double *a, double *b, double *x );
int zerosAlongDiagonal( int N, double *a );
int converged( int N, double *a, double *b );
void *ils_thread_worker();

/* Struct for passing data to thread functions */

struct args {
    int N; 
    int startRow;
    int stopRow;
    double *Aptr;
    double *Bptr;
    double *xptr;
    double *x0ptr;
    double *sum1;
    double *sum2;
    int threadid;
};

void ils_( int *threads, int *len, double *a, double *b, double *x ){

    // This function has to break up the data, spawn the processes, gather the results, and 
    // clean up after itself.


    int numThreads = *threads;
    int N = *len;
    int *numberOfRows;
    int startRow, stopRow;
    pthread_t *thread_id;
    struct args *thread_args;


    	int i, j, k, iteration;
	double sum1, sum2;
	double *sum1Ptr, *sum2Ptr;
	double ZERO = 0.0;
	int ITERATION_MAX = 2000;
	double *x0;

    // if there are fewer dimensions than threads, do a simple single-threaded matrix multiplication.
    if ( N < numThreads ) {
	    if ( ! zerosAlongDiagonal( N, a ) ) {

        // Do Jacobi Iterative Method to solve Ax=b. 

        // Create a temporary vector to hold initial values and intermediate steps 

        x0 = malloc( N * sizeof(double) );

        // Fill the x0 vector with initial values of zero

        for (i=0;i<N;i++) *(x+i) = 0.0;

        // Fill the x vector with b vector just so the initial convergence test will fail

        for (i=0;i<N;i++) *(x0+i) = *(b+i);

        // If more than N/3 iterations are done, the direct solver is more efficient
        ITERATION_MAX = fmax(ITERATION_MAX, N/3);

        iteration = 0;
        while ( !converged(N,x,x0) && iteration < ITERATION_MAX ) {

            // copy last result to initial values

            for (i=0;i<N;i++) *(x0+i) = *(x+i);

            // start the reduction process  (ref: Golub and van Loan, Chapter 10)

            for (i=0;i<N;i++) { 
                sum1 = 0.0;
                for (j=0;j<i-1;j++) sum1+= *(a+i*N+j)* *(x0+j); 
                sum2 = 0.0; 
                for (j=i+1;j<N;j++) sum2+= *(a+i*N+j)* *(x0+j); 
                *(x+i) = ( *(b+i) - sum1 - sum2 ) / *(a+i*N+i);
            }

            iteration++;

        }

        // the initial value array is no longer needed
        free(x0);

        if ( iteration == ITERATION_MAX) {
            printf(" *** ITERATIVE SOLVER FAILED TO REACH CONVERGENCE AFTER  ***\n");
            printf(" *** %d ITERATIONS, SWITCHING TO DIRECT SOLVER ***\n", iteration);
            dls_( threads, len, a, b, x );
        }

    }

    else {

        printf(" *** FOUND A ZERO ELEMENT ALONG MATRIX DIAGONAL ***\n");
        printf(" ***  SWITCHING TO DIRECT SOLVER FOR PIVOTING   ***\n");
        dls_( threads, len, a, b, x );

    }

    } // end if n < numThreads

    else { 

        /* 
         *  The parallel work begins.  The process is to first determine how to break
         *  up the matrix work equitably across the threads.  Once this is done the struct is filled with
         *  the information and a thread is started using the information.  Other than the size of the
         *  matrices and the rows to be processed, only the pointers to the memory locations of the matrices
         *  are passed.
         */ 
    if ( ! zerosAlongDiagonal( N, a ) ) {

        // Do Jacobi Iterative Method to solve Ax=b. 

        // Create a temporary vector to hold initial values and intermediate steps 

        x0 = malloc( N * sizeof(double) );

        // Fill the x0 vector with initial values of zero

        for (i=0;i<N;i++) *(x+i) = 0.0;

        // Fill the x vector with b vector just so the initial convergence test will fail

        for (i=0;i<N;i++) *(x0+i) = *(b+i);

        // If more than N/3 iterations are done, the direct solver is more efficient
        ITERATION_MAX = fmax(ITERATION_MAX, N/3);

        iteration = 0;
        while ( !converged(N,x,x0) && iteration < ITERATION_MAX ) {

            // copy last result to initial values

            for (i=0;i<N;i++) *(x0+i) = *(x+i);


	sum1Ptr = ( double * ) malloc ( numThreads * sizeof(double));
	sum2Ptr = ( double * ) malloc ( numThreads * sizeof(double));
            // start the reduction process  (ref: Golub and van Loan, Chapter 10)
        // Malloc an array to keep up with thread id's for each thread
        thread_id = (pthread_t *) malloc (numThreads * sizeof(pthread_t));

        // Malloc an array to keep up with how many rows to work on in each thread
        numberOfRows = ( int * ) malloc( numThreads * sizeof(int) );

        // Here we detemine the number of rows over which each thread will work
        for (int i=0; i<numThreads; i++ ){
            *(numberOfRows+i) = N / numThreads;
        }
        for (int i=0; i< N % numThreads; i++ ){
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
                thread_args->N   = N;
                thread_args->startRow = startRow;
                thread_args->stopRow = stopRow; 
                thread_args->Aptr = a;
                thread_args->Bptr = b;
                thread_args->xptr = x;
		thread_args->x0ptr= x0;
                thread_args->sum1 = sum1Ptr;
                thread_args->sum2 = sum2Ptr;
                thread_args->threadid = i;

                pthread_create( thread_id+i, NULL, &ils_thread_worker, thread_args );
            }
        }

        for(int i=0; i < numThreads ; i++) {
            pthread_join( *(thread_id+i), NULL); 
        }

        free(numberOfRows);
        free(thread_id);
        
            iteration++;

        }
        
                // the initial value array is no longer needed
        free(x0);


        if ( iteration == ITERATION_MAX) {
            printf(" *** ITERATIVE SOLVER FAILED TO REACH CONVERGENCE AFTER  ***\n");
            printf(" *** %d ITERATIONS, SWITCHING TO DIRECT SOLVER ***\n", iteration);
            dls_( threads, len, a, b, x );
        }

    }

    else {

        printf(" *** FOUND A ZERO ELEMENT ALONG MATRIX DIAGONAL ***\n");
        printf(" ***  SWITCHING TO DIRECT SOLVER FOR PIVOTING   ***\n");
        dls_( threads, len, a, b, x );

    }

    }

    // END OF MMM -- the memory pointed to by *C should now contain the product A.B
}


void *ils_thread_worker( struct args *thread_args  ) {

    // This is where the actual work of multiplying each matrix takes place --
    // with each mmm_thread_worker only responsible for a piece of the matrix A.
    
    // It is worth noting here that while the thread workers are reading simultaneously
    // from *A and *B, they never write to same memory locations in *C.

    int i, j, k, threadid;
    double val;
    int rowStart, rowStop, N; 
    double *a, *b, *x0, *x, *result1, *result2;
    double sum1, sum2;

    // Unpack the thread_args struct into normal variables
    N        =  thread_args->N;
    rowStart =  thread_args->startRow;
    rowStop  =  thread_args->stopRow; 
    a        =  thread_args->Aptr;
    b        =  thread_args->Bptr;
    x        =  thread_args->xptr;
    x0       =  thread_args->x0ptr;
    result1  =  thread_args->sum1;
    result2  =  thread_args->sum2;
    threadid =  thread_args->threadid;

    // Process the rows for which this thread is responsible
    for(k=rowStart; k<rowStop; k++){
            for (i=0;i<N;i++) { 
                sum1 = 0.0;
                for (j=0;j<i-1;j++) sum1+= *(a+i*N+j)* *(x0+j); 
                sum2 = 0.0; 
                for (j=i+1;j<N;j++) sum2+= *(a+i*N+j)* *(x0+j); 
                *(x+i) = ( *(b+i) - sum1 - sum2 ) / *(a+i*N+i);
            } 
    }

    result1[threadid] = sum1;
    result2[threadid] = sum2;

    free(thread_args);
    pthread_exit(NULL);
}

// Code to check for zeros along the diagonal
int zerosAlongDiagonal ( int N, double *a ) {

    double ZERO;
    int i;
    int foundZero;

    foundZero = 0;
    for (i=0;i<N;i++) { 
        if (!foundZero)  
            foundZero = fabs(*(a+i*N+i)) == ZERO;
    }
    return foundZero;
}

// Code to check for convergence (A. Pounds, 2018)
int converged( int N, double *a, double *b) {

    // Compute the distance between the vectors and see if the 2-Norm is
    // within tolerance

    double const TOL = 5.0e-15;
    double sum, maxb;
    int i;

    // find max in array b for tolerance scaling while computing sum

    maxb=*(b+0); 
    sum = 0.0; 
    for (i=0; i<N; i++) {
        maxb = fmax(maxb,fabs(*(b+i)));
        sum += (*(a+i)-*(b+i))*(*(a+i)-*(b+i));
    }
    sum = sqrt(sum);

    // by dividing by the largest value in the b matrix we effectively
    // scale the 2-Norm so it can achieve machine precision
    return (sum/maxb < TOL);    

}


