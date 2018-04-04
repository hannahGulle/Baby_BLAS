program driver 

integer :: NDIM

real (kind=8) :: wall_start, wall_end
real (kind=8) :: cpu_start, cpu_end
real (kind=8) :: trace


integer :: startval, stopval, stepval, nthreads
real (kind=8) :: walltime
real (kind=8) :: cputime 
external walltime, cputime

character (len=8) :: carg1, carg2, carg3, carg4

real (kind=8), dimension(:), allocatable :: veca, vecb
real (kind=8), dimension(:,:), allocatable :: matrixa, matrixb, matrixc

#ifdef ACCURACY_TEST

!This portion of code is ONLY used for verifying the accuracy of the code using
!the matrix and matrix inverse stored on the class website.

!Download the files from theochem using curl (don't store these on anvil!)
call system("curl -s -o matrixa.dat --url http://theochem.mercer.edu/csc435/data/matrixa.dat")
call system("curl -s -o matrixb.dat --url http://theochem.mercer.edu/csc435/data/matrixb.dat")

NDIM = 100  ! The test files are 100x100 double precision matrix and its inverse
nthreads = 2
allocate ( matrixa(NDIM,NDIM), stat=ierr)
allocate ( matrixb(NDIM,NDIM), stat=ierr)
allocate ( matrixc(NDIM,NDIM), stat=ierr)
open (unit=5,file="matrixa.dat",status="old")
do i = 1, NDIM
  do j = 1, NDIM
     read(5,*) matrixa(j,i)
  enddo
enddo
close(5)
open (unit=5,file="matrixb.dat",status="old")
do i = 1, NDIM
  do j = 1, NDIM
     read(5,*) matrixb(j,i)
  enddo
enddo
close(5)

! Delete the files from disk
call system("rm matrixa.dat matrixb.dat")

#else

! Start the normal processing here.  Read the starting, stop, and step values
! as well as the number of threads to use.
! modified to use command line arguments

call get_command_argument(1, carg1)
call get_command_argument(2, carg2)
call get_command_argument(3, carg3)
call get_command_argument(4, carg4)

! Use Fortran internal files to convert command line arguments to ints

read (carg1,'(i8)') startval
read (carg2,'(i8)') stopval
read (carg3,'(i8)') stepval
read (carg4,'(i8)') nthreads 

! Start the outermost loop to run tests of varying matrix size
 
do iter = startval, stopval, stepval

NDIM = iter

allocate ( veca(NDIM), stat=ierr)
allocate ( vecb(NDIM), stat=ierr)
allocate ( matrixa(NDIM,NDIM), stat=ierr)
allocate ( matrixb(NDIM,NDIM), stat=ierr)
allocate ( matrixc(NDIM,NDIM), stat=ierr)

! Build veca and vecb which, their tensor product creates the two matrices 
! to be multiplied.

do i = 1, NDIM 
     veca(i) = 1.0
     vecb(i) = 1.0 / sqrt( dble(NDIM))
enddo

! Zero the matrices using Fortran 90 syntax.
matrixa = 0.0
matrixb = 0.0

call vvm(NDIM, veca, vecb, matrixa)
call vvm(NDIM, veca, vecb, matrixb)

#endif

wall_start = walltime()
cpu_start = cputime()

call mmm(nthreads, NDIM, matrixa, matrixb, matrixc)

cpu_end = cputime()
wall_end = walltime()

trace = 0.0;

do i=1, NDIM 
     trace = trace + matrixc(i,i)
enddo

! Calculate megaflops based on CPU time and Walltime

mflops  = 2*dble(NDIM)**3/ (cpu_end-cpu_start) / 1.0e6
mflops2 = 2*dble(NDIM)**3/ (wall_end-wall_start)/ 1.0e6
 
print *, NDIM, trace, cpu_end-cpu_start, wall_end-wall_start,  mflops, mflops2


! Free the memory that was allocated based on which version of the program was
! run.

deallocate(matrixa)
deallocate(matrixb)
deallocate(matrixc)
#ifndef ACCURACY_TEST
deallocate(veca)
deallocate(vecb)

enddo
#endif

end program driver 
 
