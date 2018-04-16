program driver 
#include "f90papi.h"
!#include "omp.h"

!!---------------------------------------------
!! VARIABLE DECLARATIONS
!!---------------------------------------------
integer :: NDIM

real (kind=8) :: wall_start, wall_end
real (kind=8) :: cpu_start, cpu_end
real (kind=8) :: trace

integer :: startval, stopval, stepval, threads
real (kind=8) :: walltime
real (kind=8) :: cputime 
external walltime, cputime

character (len=8) :: carg1, carg2, carg3, carg4

real (kind=8), dimension(:), allocatable :: veca, vecb, vecx
real (kind=8), dimension(:,:), allocatable :: matrixa, matrixb, matrixc
real (kind=8) :: dot, dotProd
external dot

integer, parameter :: NUM_EVENTS = 2
real (kind=4) :: rtime, ptime, mflops, mflips
integer (kind=8) :: flpops, flpins
integer (kind=4) :: check, eventSet
integer (kind=8), dimension (NUM_EVENTS) :: dp_ops
integer :: native

!! RETRIEVE COMMAND LINE ARGUMENTS

call get_command_argument(1, carg1)
call get_command_argument(2, carg2)
call get_command_argument(3, carg3)
call get_command_argument(4, carg4)

! Use Fortran internal files to convert command line arguments to ints
read (carg1,'(i8)') startval
read (carg2,'(i8)') stopval
read (carg3,'(i8)') stepval
read (carg4,'(i8)') threads



!!----------------------------------------------
!! DIMENSION ITERATION LOOP
!!----------------------------------------------
do iter = startval, stopval, stepval

NDIM = iter

allocate ( veca(NDIM), stat=ierr)
allocate ( vecb(NDIM), stat=ierr)
allocate ( vecx(NDIM), stat=ierr)
allocate ( matrixa(NDIM,NDIM), stat=ierr)
allocate ( matrixb(NDIM,NDIM), stat=ierr)
allocate ( matrixc(NDIM,NDIM), stat=ierr)

!! ----------------------------------------------------
!! ACCURACY TESTING BLOCK
!! ----------------------------------------------------

#ifdef MMM
do i = 1, NDIM 
     veca(i) = 1.0
     vecb(i) = 1.0 / sqrt( dble(NDIM))
enddo

matrixa = 0.0
matrixb = 0.0

call vvm(threads, NDIM, veca, vecb, matrixa);
call vvm(threads, NDIM, veca, vecb, matrixb);
#endif

#ifdef VVM
do i = 1, NDIM
    vecb(i) = dble(NDIM)
    veca(i) = 1.0 / dble(NDIM)
enddo
#endif

#ifdef DOT
do i = 1, NDIM
    veca(i) = dble(NDIM)
    vecb(i) = 1.0 / dble(NDIM)
enddo
#endif

#ifdef MVV
do i = 1, NDIM
    veca(i) = 1.0 / dble(NDIM)
    vecx(i) = 0.0
enddo
matrixa = 1.0
#endif


wall_start = walltime()
cpu_start = cputime()


#ifdef MMM
call mmm(threads, NDIM, matrixa, matrixb, matrixc);
#endif

#ifdef VVM
call vvm(threads, NDIM, veca, vecb, matrixc);
#endif

#ifdef DOT
dotProd = dot(threads, NDIM, veca, vecb);
#endif

#ifdef MVV
call mvv(threads, NDIM, matrixa, veca, vecx);
#endif


cpu_end = cputime()
wall_end = walltime()


!! ------------------------------------------------------
!! TRACE BLOCK
!! ------------------------------------------------------

#ifdef MMM
trace = 0.0
do i=1, NDIM 
     trace = trace + matrixc(i,i)
enddo
#endif

#ifdef VVM
trace = 0.0
do i=1, NDIM 
     trace = trace + matrixc(i,i)
enddo
#endif

#ifdef MVV
trace = 0.0
do i=1, NDIM
    trace = trace + vecx(i)
enddo
#endif

!! -----------------------------------------------------
!! END TRACE BLOCK
!! ----------------------------------------------------



!! -----------------------------------------------------
!! RESULTS AND DEALLOCATION BLOCK
!! -----------------------------------------------------

mflops = (2.0/3.0)*dble(NDIM)**3/(wall_end-wall_start)/1.0e6
#ifndef DOT
print *, NDIM, trace, cpu_end-cpu_start, wall_end-wall_start, threads
#else
print *, NDIM, dotProd, cpu_end-cpu_start, wall_end-wall_start, mflops
#endif

if (allocated(matrixa)) deallocate(matrixa)
if (allocated(matrixb)) deallocate(matrixb)
if (allocated(matrixc)) deallocate(matrixc)
if (allocated(veca))    deallocate(veca)
if (allocated(vecb))    deallocate(vecb)
if (allocated(vecx))    deallocate(vecx)

!! -----------------------------------------------------
!! END RESULTS AND DEALLOCATION BLOCK
!! ----------------------------------------------------

enddo
!! END ITERATION LOOP

end program driver 
 
