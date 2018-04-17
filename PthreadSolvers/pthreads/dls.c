
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
    void dls_( int *threads, int *len,  double *a, double *b, double*x );
#ifdef __cplusplus
}
#endif

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>

int strictlyDiagonallyDominant( int N, double *a ); 
void *dls_thread_worker();

/* Struct for passing data to thread functions */

struct args {
    int N; 
    int k;
    int startRow;
    int stopRow;
    int indexPivot;
    double *Aptr;
    double pivot;
};

void dls_( int *threads, int *len, double *a, double *b, double *x ){

    // This function has to break up the data, spawn the processes, gather the results, and 
    // clean up after itself.


    int numThreads = *threads;
    int N = *len;
    int *numberOfRows;
    int startRow, stopRow;
    pthread_t *thread_id;
    struct args *thread_args;

    int i, j, k, u;
    int singular, iPivot, rows, rows2;
    double pivotMax, tmp, *y;
    double sum;
    double ZERO = 0.0;
    int *p;


    // if there are fewer dimensions than threads, do a simple single-threaded matrix multiplication.
    if ( N < numThreads ) {
	
      // Check A for strict diagonal dominance to see if we can reduce the matrix without 
      // doing any row interchanges.   We could also check for positive definiteness to
      // achieve the same thing.
  
      if ( ! strictlyDiagonallyDominant( N, a ) ) {
  
          // Do Gaussian Elimination with Partial Pivoting 
          //   (modified from Golub and van Load, Chapter 3) 
  
          // Create an array to hold pivot swaps 
  
          p = malloc( (N-1) * sizeof(int) );
  
          for (k=0;k<N-1;k++) *(p+k)=k;
  
          // Search for largest value in the column and swap the 
          // entire row containing that value with the current
          // pivot row.
  
          for (k=0;k<N-1;k++) {
              pivotMax = *(a+k*N+k);
              iPivot = k; 
              for (u=k;u<N;u++) {
                  if ( fabs(*(a+u*N+k)) > fabs(pivotMax) ) {
                      pivotMax = *(a+u*N+k);
                      iPivot = u;
                  }
              }
              // If a greater pivot value was found, swap the rows.
              if ( iPivot != k ) {
                  u = iPivot; 
                  for (j=k;j<N;j++) {
                      tmp = *(a+k*N+j);
                      *(a+k*N+j) = *(a+u*N+j);
                      *(a+u*N+j)=tmp;
                  }
              }
  
              // Now do block reduction
              *(p+k) = iPivot;
              if ( *(a+k*N+k) != ZERO ) {
                  for (rows=k+1;rows<N;rows++) { 
                      *(a+rows*N+k) = *(a+rows*N+k) / *(a+k*N+k);
  
                      for (rows2=k+1;rows2<N;rows2++) { 
                          *(a+rows*N+rows2) = *(a+rows*N+rows2) - 
                              *(a+rows*N+k) * *(a+k*N+rows2) ;
                      }
                  }
              }
  
              else {
  
                  /* Handle the case of a zero pivot element, singular matrix */
  
                  printf( "Element a[%d][%d} = %f\n", k, k, *(a+k*N+k)); 
                  printf( " *** MATRIX A IS SINGULAR *** \n");
                  printf( "    -- EXECUTION HALTED --\n");
                  exit(1);
              }
  
          } // End of Outer For Loop
          // Now that we know we have reduced the matrices, start the 
          // back substitution process to solve for vector x.
  
  
          /* We now need to arrange b so that it has undergone the same 
           * operations as the matrix a.  This will then form
           * the vector y for the solution of Ux=y where U is the 
           * upper-triangular matrix formed in the elimination process
           * above. 
           */
  
          for (k=0; k<N-1; k++ ) {
              // Swap rows x with p(k) 
              tmp = *(b+k);
              *(b+k) = *(b+ *(p+k));
              *(b+ *(p+k)) = tmp;
  
              for (j=k+1;j<N;j++) 
                  *(b+j)= *(b+j) - *(b+k) * *(a+N*j+k);  
          } 
  
          // Now do the backward substitution to get the solution
          // vector x
  
          *(b+N-1) = *(b+N-1) / *(a+N*(N-1)+(N-1));
          for (i=N-2;i>=0;i--){
              tmp = 0.0;
              for (j=i+1;j<N;j++) {
                  tmp = tmp + *(a+i*N+j) * *(b+j);
              }
              *(b+i) = ( *(b+i) - tmp ) / *(a+i*N+i); 
          }
  
          for (i=0;i<N;i++) *(x+i) = *(b+i);
  
          // At this point the solution to the system should be in vector x 
  
          free(p);
      } // end of if diagonally dominant if
  
      else {
  
          // Since we know the matrix is strictly diagonally dominant, verify
          // that none of the pivot elements are equal to zero
  
          singular = 1; 
          i=0;
          while ( i<N  && singular ) {
              singular = *(a+i*N+i) == ZERO;   
              i++;
          }
  
          if ( singular ) {
              printf( " *** MATRIX A IS SINGULAR *** \n");
              printf( "    -- EXECUTION HALTED -- \n");
              exit(1);
          }
  
          // We know at this point that we have a strictly diagonally dominant matrix that is
          // not singular -- so it sould be possible to do the LU factorization.
          //   (modified from Golub and van Loan, Chapter 3)
  
          for (k=0; k<N-1; k++) {
              for (rows=k+1;rows<N;rows++) {
                  *(a+rows*N+k) = *(a+rows*N+k) / *(a+k*N+k);
  
                  for (rows2=k+1;rows2<N;rows2++) { 
                      *(a+rows*N+rows2) = *(a+rows*N+rows2) - 
                          *(a+rows*N+k) * *(a+k*N+rows2) ;
                  }
              }
          }
  
          // At this point the LU factorizaton should be done and we have to do two
          // triangular back substitutions.  The solution to Ax=b is solved by first 
          // solving Ly=b for y and then Ux=y for the solution vector x.
  
          // Solving lower-triangular (Ly=b) first, overwriting b with y
  
          for (k=0; k<N-1; k++ ) {
              for (j=k+1;j<N;j++) 
                  *(b+j)= *(b+j) - *(b+k) * *(a+N*j+k);  
          } 
  
          // Now we can do the backward substitution to get the solution
          // vector x for the upper-triangular system (Ux=y) overwriting y (stored in b)
          // with x
  
          *(b+N-1) = *(b+N-1) / *(a+N*(N-1)+(N-1));
          for (i=N-2;i>=0;i--){
              tmp = 0.0;
              for (j=i+1;j<N;j++) {
                  tmp = tmp + *(a+i*N+j) * *(b+j);
              }
              *(b+i) = ( *(b+i) - tmp ) / *(a+i*N+i); 
          }
  
          for (i=0;i<N;i++) *(x+i) = *(b+i);
  
          // At this point the solution to the system should be in vector x 
  
    } 
    }
    else { //end of if

        /* 
         *  The parallel work begins.  The process is to first determine how to break
         *  up the matrix work equitably across the threads.  Once this is done the struct is filled with
         *  the information and a thread is started using the information.  Other than the size of the
         *  matrices and the rows to be processed, only the pointers to the memory locations of the matrices
         *  are passed.
         */ 
    if ( ! strictlyDiagonallyDominant( N, a ) ) {

        // Do Gaussian Elimination with Partial Pivoting 
        //   (modified from Golub and van Load, Chapter 3) 

        // Create an array to hold pivot swaps 

        p = malloc( (N-1) * sizeof(int) );

        for (k=0;k<N-1;k++) *(p+k)=k;

        // Search for largest value in the column and swap the 
        // entire row containing that value with the current
        // pivot row.

        for (k=0;k<N-1;k++) {
            pivotMax = *(a+k*N+k);
            iPivot = k; 
            
            
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
                  thread_args->pivot = pivotMax;
                  thread_args->indexPivot = iPivot;
  		  thread_args->k = k;
                  pthread_create( thread_id+i, NULL, &dls_thread_worker, thread_args );
              }
          }
          for(int i=0; i < numThreads ; i++) {
              pthread_join( *(thread_id+i), NULL); 
          }
  
          free(numberOfRows);
          free(thread_id);
  
  
          //----------------------------------------------------------------------------
          //--------------- SERIAL CODE BLOCK STARTS HERE ------------------------------
          //----------------------------------------------------------------------------
  
                      // If a greater pivot value was found, swap the rows.
            if ( iPivot != k ) {
                u = iPivot; 
                for (j=k;j<N;j++) {
                    tmp = *(a+k*N+j);
                    *(a+k*N+j) = *(a+u*N+j);
                    *(a+u*N+j)=tmp;
                }
            }

            // Now do block reduction
            *(p+k) = iPivot;
            if ( *(a+k*N+k) != ZERO ) {
                for (rows=k+1;rows<N;rows++) { 
                    *(a+rows*N+k) = *(a+rows*N+k) / *(a+k*N+k);

                    for (rows2=k+1;rows2<N;rows2++) { 
                        *(a+rows*N+rows2) = *(a+rows*N+rows2) - 
                            *(a+rows*N+k) * *(a+k*N+rows2) ;
                    }
                }
            }

            else {

                /* Handle the case of a zero pivot element, singular matrix */

                printf( "Element a[%d][%d} = %f\n", k, k, *(a+k*N+k)); 
                printf( " *** MATRIX A IS SINGULAR *** \n");
                printf( "    -- EXECUTION HALTED --\n");
                exit(1);
            }

        }// end of for with threading
        // Now that we know we have reduced the matrices, start the 
        // back substitution process to solve for vector x.


        /* We now need to arrange b so that it has undergone the same 
         * operations as the matrix a.  This will then form
         * the vector y for the solution of Ux=y where U is the 
         * upper-triangular matrix formed in the elimination process
         * above. 
         */

        for (k=0; k<N-1; k++ ) {
            // Swap rows x with p(k) 
            tmp = *(b+k);
            *(b+k) = *(b+ *(p+k));
            *(b+ *(p+k)) = tmp;

            for (j=k+1;j<N;j++) 
                *(b+j)= *(b+j) - *(b+k) * *(a+N*j+k);  
        } 

        // Now do the backward substitution to get the solution
        // vector x

        *(b+N-1) = *(b+N-1) / *(a+N*(N-1)+(N-1));
        for (i=N-2;i>=0;i--){
            tmp = 0.0;
            for (j=i+1;j<N;j++) {
                tmp = tmp + *(a+i*N+j) * *(b+j);
            }
            *(b+i) = ( *(b+i) - tmp ) / *(a+i*N+i); 
        }

        for (i=0;i<N;i++) *(x+i) = *(b+i);

        // At this point the solution to the system should be in vector x 

        free(p);
    }//end of if diagonally dominant

    else {

        // Since we know the matrix is strictly diagonally dominant, verify
        // that none of the pivot elements are equal to zero

        singular = 1; 
        i=0;
        while ( i<N  && singular ) {
            singular = *(a+i*N+i) == ZERO;   
            i++;
        }

        if ( singular ) {
            printf( " *** MATRIX A IS SINGULAR *** \n");
            printf( "    -- EXECUTION HALTED -- \n");
            exit(1);
        }

        // We know at this point that we have a strictly diagonally dominant matrix that is
        // not singular -- so it sould be possible to do the LU factorization.
        //   (modified from Golub and van Loan, Chapter 3)

        for (k=0; k<N-1; k++) {
            for (rows=k+1;rows<N;rows++) {
                *(a+rows*N+k) = *(a+rows*N+k) / *(a+k*N+k);

                for (rows2=k+1;rows2<N;rows2++) { 
                    *(a+rows*N+rows2) = *(a+rows*N+rows2) - 
                        *(a+rows*N+k) * *(a+k*N+rows2) ;
                }
            }
        }

        // At this point the LU factorizaton should be done and we have to do two
        // triangular back substitutions.  The solution to Ax=b is solved by first 
        // solving Ly=b for y and then Ux=y for the solution vector x.

        // Solving lower-triangular (Ly=b) first, overwriting b with y

        for (k=0; k<N-1; k++ ) {
            for (j=k+1;j<N;j++) 
                *(b+j)= *(b+j) - *(b+k) * *(a+N*j+k);  
        } 

        // Now we can do the backward substitution to get the solution
        // vector x for the upper-triangular system (Ux=y) overwriting y (stored in b)
        // with x

        *(b+N-1) = *(b+N-1) / *(a+N*(N-1)+(N-1));
        for (i=N-2;i>=0;i--){
            tmp = 0.0;
            for (j=i+1;j<N;j++) {
                tmp = tmp + *(a+i*N+j) * *(b+j);
            }
            *(b+i) = ( *(b+i) - tmp ) / *(a+i*N+i); 
        }

        for (i=0;i<N;i++) *(x+i) = *(b+i);

        // At this point the solution to the system should be in vector x 


    }//end if not diagonally dominant

  
          }//end of use threading?
    }//end of function



void *dls_thread_worker( struct args *thread_args  ) {

    int i, j, k;
    double val;
    int rowStart, rowStop, N, iPivot; 
    double *a, max;

    // Unpack the thread_args struct into normal variables
    N        =  thread_args->N;
    k	     =  thread_args->k;
    rowStart =  thread_args->startRow;
    rowStop  =  thread_args->stopRow; 
    a        =  thread_args->Aptr;
    max      =  thread_args->pivot;
    iPivot   =  thread_args->indexPivot;

    // Process the rows for which this thread is responsible
	for(i=rowStart;i<rowStop;i++){
	    for (j=k;j<N;j++) {
        	 if ( fabs(*(a+j*N+k)) > fabs(max) ) {
              		max = *(a+j*N+k);
              		iPivot = j;
    	     	}
    	    }
	}

    free(thread_args);
    pthread_exit(NULL);
}

int strictlyDiagonallyDominant( int N, double *a ) {

    double sum;
    int i, testPassed, row;

    testPassed = 1;
    row = 0;
    sum = 0.0;
    for (row=0;row<N;row++) { 
        if (testPassed) {
            sum = 0.0;
            for (i=0;i<N;i++) sum+=*(a+row*N+i);
            sum-=fabs(*(a+row*N+row)); 
            testPassed = fabs(*(a+row*N+row)) > sum;
        }
    }

    return testPassed;
}




