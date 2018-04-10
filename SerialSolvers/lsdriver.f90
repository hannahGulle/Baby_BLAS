program driver 
#include "f90papi.h"

integer :: NDIM

real (kind=8) :: wall_start, wall_end
real (kind=8) :: cpu_start, cpu_end

integer :: nthreads
real (kind=8) :: walltime
real (kind=8) :: cputime 
external walltime, cputime

real (kind=8), dimension(:), allocatable :: veca, vecb, vecx
real (kind=8), dimension(:,:), allocatable :: matrixa, matrixb, matrixc
logical :: DIAG_DOMINANT, SPARSE_MATRIX
real (kind=8) :: residual

integer,parameter :: NUM_EVENTS = 2
real (kind=4) :: rtime, ptime, mflops, mflips
integer (kind=8) :: flpops, flpins
integer (kind=4) :: check, eventSet
integer (kind=8), dimension (NUM_EVENTS) :: dp_ops


!! --------------------------------------------------------------------------
!! -------------- PAPI BLOCK ------------------------------------------------
!! --------------------------------------------------------------------------
! Initialize Papi Check Variable (Initialize Library)
check = PAPI_VER_CURRENT

call PAPIF_library_init(check);
if ((check .ne. PAPI_VER_CURRENT) .and. (check .gt. 0)) then
    print *, "Papi Library Version Mismatch!"
    call exit()
endif

if (check .lt. 0) then
    print *, "Papi Initialization Error."
    call exit()
endif

call PAPIF_is_initialized(check);
if (check .ne. PAPI_LOW_LEVEL_INITED) then
    print *, "Papi Low Level Initialization Error."
    call exit()
endif

eventSet = PAPI_NULL;

call PAPIF_create_eventset(eventSet, check);
if (check .ne. PAPI_OK) then
    print *, "Could Not Create Papi Event Set."
    call exit()
endif

! Add the events to be counted in each set
call PAPIF_add_event(eventSet, PAPI_VEC_DP, check)
if (check .ne. PAPI_OK) then
    print *, "Could Not Create PAPI_VEC_DP Event."
    call exit()
endif

!! ---------------------------------------------------------------------------
!! -------------- END PAPI BLOCK ---------------------------------------------
!! ---------------------------------------------------------------------------


DIAG_DOMINANT = .false.
SPARSE_MATRIX = .false.

NDIM = 100 
nthreads = 1

#ifdef ILS
print *, "Performing Iterative Solver Accuracy Test"
#else
print *, "Performing Direct Solver Accuracy Test"
#endif

!This portion of code is ONLY used for verifying the accuracy of the code using
!the matrix, vector b, and solution vector x stored on the class website.

!Download the files from theochem using curl (don't store these on anvil!)
!NOTE: for strictly diagonally dominant systems append _dd to last file name, e.g. -- linsolve_a_dd.dat
#ifdef DIAGDOM
call system("curl -s -o linsolve_a.dat --url http://theochem.mercer.edu/csc435/data/linsolve_a_dd.dat")
call system("curl -s -o linsolve_b.dat --url http://theochem.mercer.edu/csc435/data/linsolve_b_dd.dat")
call system("curl -s -o linsolve_x.dat --url http://theochem.mercer.edu/csc435/data/linsolve_x_dd.dat")
#else
call system("curl -s -o linsolve_a.dat --url http://theochem.mercer.edu/csc435/data/linsolve_a.dat")
call system("curl -s -o linsolve_b.dat --url http://theochem.mercer.edu/csc435/data/linsolve_b.dat")
call system("curl -s -o linsolve_x.dat --url http://theochem.mercer.edu/csc435/data/linsolve_x.dat")
#endif

print *, "Files loaded from theochem.mercer.edu"

allocate ( matrixa(NDIM,NDIM), stat=ierr)
allocate ( matrixb(NDIM,NDIM), stat=ierr)
allocate ( matrixc(NDIM,NDIM), stat=ierr)
allocate ( veca(NDIM), stat=ierr)
allocate ( vecb(NDIM), stat=ierr)
allocate ( vecx(NDIM), stat=ierr)

open (unit=5,file="linsolve_a.dat",status="old")
do i = 1, NDIM
  do j = 1, NDIM
     read(5,*) matrixa(j,i)
  enddo
enddo
close(5)
open (unit=5,file="linsolve_b.dat",status="old")
do i = 1, NDIM
   read(5,*) vecb(i)
enddo
close(5)
open (unit=5,file="linsolve_x.dat",status="old")
do i = 1, NDIM
   read(5,*) veca(i)
enddo
close(5)

print *, "Files read into program"

! Delete the files from disk
call system("rm linsolve_a.dat linsolve_b.dat linsolve_x.dat")

print *, "Files deleted from disk."


! If doing ILS or DLS testing, build matrix C explicitly 
! as well as solution vector X and product vector B. You
! can also specify if you want the system to be diagonally 
! dominant.

#ifdef DIAGDOM
DIAG_DOMINANT = .true.
#elif SPARSE
SPARSE_MATRIX = .true.
#endif

call buildLinearSystem( NDIM, matrixa, vecb, vecx,  DIAG_DOMINANT, SPARSE_MATRIX )

! Start the counters in each event set
call PAPIF_start(eventSet, check);
if (check .ne. PAPI_OK) then
    print *, "Could Not Start PAPI_VEC_DP Counter."
    call exit()
endif

wall_start = walltime()
cpu_start = cputime()

! Read set and set eventSet back to zero
call PAPIF_accum(eventSet, dp_ops, check);
if (check .ne. PAPI_OK) then
    print *, "Could Not Accumulate Papi Event Set."
    call exit()
endif

#ifdef DLS
call dls(nthreads, NDIM, matrixa, vecb, vecx)
#else
call ils(nthreads, NDIM, matrixa, vecb, vecx)
#endif

! Read event set
call PAPIF_read(eventSet, dp_ops, check);
if (check .ne. PAPI_OK) then
    print *, "Could Not Read Papi Event Set."
    call exit()
endif

cpu_end = cputime()
wall_end = walltime()

call PAPIF_flops( rtime, ptime, flpops, mflops, check );
call PAPIF_flips( rtime, ptime, flpins, mflips, check );

residual = 0.0
do i=1, NDIM
   residual = max(residual, abs(vecx(i)-dble(i)))
enddo

mflops = (dp_ops(1)/(cpu_end-cpu_start))/1.0e6
mflops2  = (2.0/3.0)*dble(NDIM)**3/ (wall_end-wall_start) / 1.0e6

print *, "Dimension     Residual                  CpuTime                   WallTime                  Est.Flop          PapiFlops" 
print *, NDIM, residual, cpu_end-cpu_start, wall_end-wall_start,  mflops2, mflops

! Free the memory that was allocated based on which version of the program was
! run.

if (allocated(matrixa)) deallocate(matrixa)
if (allocated(matrixb)) deallocate(matrixb)
if (allocated(matrixc)) deallocate(matrixc)
if (allocated(veca)) deallocate(veca)
if (allocated(vecb)) deallocate(vecb)
if (allocated(vecx)) deallocate(vecx)


end program driver 


! Subroutine to build random linear systems for Solvers to Use,
!  (A. Pounds, 2018)
subroutine buildLinearSystem( N, A, B, X,  DIAG_DOMINANT, SPARSE_MATRIX )

integer :: N;
real (kind=8), dimension(N,N) :: A 
real (kind=8), dimension(N) :: B, X
logical :: DIAG_DOMINANT, SPARSE_MATRIX
real (kind=8) :: rowsum

call init_random_seed()
call random_number(A)

do i =1, N
  X(i) = dble(i)
enddo

if (DIAG_DOMINANT ) then
! Force the matrix to be diagonally dominant
    do i=1,N
   rowsum = 0.0
   do j=1, N 
      rowsum=rowsum+abs(A(j,i))
   enddo
   A(i,i) = rowsum-abs(A(i,i))+100.0
 enddo
endif

if (SPARSE_MATRIX) then
!Build sparse matrix with 5 superdiagonals and 5 subdiagonals
do i=1,N
   do j=1, N 
      if (j .lt. (i-5) .or. j .gt. (i+5) ) a(j,i) = 0.0;
   enddo
 enddo
endif

B = matmul(X,A)

end subroutine 



! This is needed for the random number generator
SUBROUTINE init_random_seed() 
  INTEGER :: i, n, clock 
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed 
 
  CALL RANDOM_SEED(size = n) 
  ALLOCATE(seed(n)) 
 
  CALL SYSTEM_CLOCK(COUNT=clock) 

  seed = clock + 37 * (/ (i - 1, i = 1, n) /)
  CALL RANDOM_SEED(PUT = seed)

  DEALLOCATE(seed)
END SUBROUTINE
 
