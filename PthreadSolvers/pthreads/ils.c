/**********************************************************************
 *
 * ITERATIVE LINEAR SOLVER
 *
 * Andrew J. Pounds, Ph.D.
 * Spring 2018
 *
 * Unless otherwise noted, all code and methods belong to the author.
 * Equations for the Jacoby iterative solver we adapted from Golub
 * and van Loan, "Matrix Computations", Johns Hopkins University press,
 * 1996. 
 *
 **********************************************************************/

#ifdef __cplusplus
extern "C" {
#endif
    void ils_( int *threads, int *len,  double *a, double *b, double *x );
#ifdef __cplusplus
}
#endif

#include <math.h>  
#include <stdio.h>  
#include <stdlib.h>

/*  S E R I A L   C O D E  */

// Need the following prototype in case of a zero along the diagonal
void  dls_( int *threads, int *len, double *a, double *b, double *x );

// Prototype for code to check for zeros along the diagonal
int zerosAlongDiagonal ( int N, double *a );

// Prototype for code to check for convergence
int converged( int N, double *a, double *b);

void ils_( int *threads, int *len,  double *a, double *b, double *x ){

    /* in serial code, *threads not used. It is retained here so the code can be called
     * identically to the threaded methods.
     */


    int i, j, k, N, iteration;
    double sum1, sum2;
    double ZERO = 0.0;
    int ITERATION_MAX = 2000;
    double *x0;

    N = *len;

    // The iterative Jacobi method will fail if there are zeros along the diagonal of the
    // matrix.  We could reorder the matrix, but in this case if zeros are found we will switch
    // to our discrete linear solver which will do a form of Gaussian Elimination with partial
    // pivoting to solve the system. 

    // NOTE: we are not modifying matrix A or vector B in any manner in this iterative procedure
    // so that if we fail to achieve convergence we can drop back to the direct solver with the original
    // A matrix and B vector.

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

