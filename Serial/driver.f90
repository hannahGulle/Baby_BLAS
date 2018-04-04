program driver 
#include "f90papi.h"

!!---------------------------------------------
!! VARIABLE DECLARATIONS
!!---------------------------------------------
integer :: NDIM

real (kind=8) :: wall_start, wall_end
real (kind=8) :: cpu_start, cpu_end
real (kind=8) :: trace

integer :: startval, stopval, stepval
real (kind=8) :: walltime
real (kind=8) :: cputime 
external walltime, cputime

character (len=8) :: carg1, carg2, carg3

real (kind=8), dimension(:), allocatable :: veca, vecb, vecx
real (kind=8), dimension(:,:), allocatable :: matrixa, matrixb, matrixc

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

! Use Fortran internal files to convert command line arguments to ints
read (carg1,'(i8)') startval
read (carg2,'(i8)') stopval
read (carg3,'(i8)') stepval
 



!!----------------------------------------------
!! DIMENSION ITERATION LOOP
!!----------------------------------------------
do iter = startval, stopval, stepval



!! -------------------------------------------
!! PAPI BLOCK
!!---------------------------------------------
! Papi Initialize Check Variable
check = PAPI_VER_CURRENT

! Papi Intitialize Library
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
    print *, "Papi Low Level Initialization Failed."
    call exit()
endif

! Create a Papi Event Set
eventSet = PAPI_NULL;

call PAPIF_create_eventset(eventSet, check)
if (check .ne. PAPI_OK) then
    print *, "Could Not Create Papi Event Set."
    call exit()
endif

! Add the particular events to be counted to each event set
call PAPIF_add_event(eventSet, PAPI_DP_OPS, check)
if (check .ne. PAPI_OK) then
    print *, "Could Not Create PAPI_DP_OPS Event."
    call exit()
endif

! Start the counters in each event set
call PAPIF_start(eventSet, check)
if (check .ne. PAPI_OK) then
    print *, "Could Not Start PAPI_DP_OPS Counter."
    call exit()
endif
!! -----------------------------------------------------
!! END PAPI BLOCK
!! -----------------------------------------------------  



NDIM = iter

allocate ( veca(NDIM), stat=ierr)
allocate ( vecb(NDIM), stat=ierr)
allocate ( matrixa(NDIM,NDIM), stat=ierr)
allocate ( matrixb(NDIM,NDIM), stat=ierr)
allocate ( matrixc(NDIM,NDIM), stat=ierr)

!! ----------------------------------------------------
!! ACCURACY TESTING BLOCK
!! ----------------------------------------------------
do i = 1, NDIM 
     veca(i) = 1.0
     vecb(i) = 1.0 / sqrt( dble(NDIM))
enddo

matrixa = 0.0
matrixb = 0.0

call vvm(NDIM, veca, vecb, matrixa);
call vvm(NDIM, veca, vecb, matrixb);

wall_start = walltime()
cpu_start = cputime()

! Read set and set array back to zero
call PAPIF_accum(eventSet, dp_ops, check);
if (check .ne. PAPI_OK) then
    print *, "Could Not Accumulate Papi Event Set."
    call exit()
endif

call mmm(NDIM, matrixa, matrixb, matrixc);

call PAPIF_read(eventSet, dp_ops, check)
if (check .ne. PAPI_OK) then
    print *, "Could Not Read Papi Event Set."
    call exit()
endif

cpu_end = cputime()
wall_end = walltime()

call PAPIF_flops( rtime, ptime, flpops, mflops, check );
call PAPIF_flips( rtime, ptime, flpins, mflips, check );

!! -----------------------------------------------------
!! END ACCURACY TESTING BLOCK
!! -----------------------------------------------------





!! ------------------------------------------------------
!! TRACE BLOCK
!! ------------------------------------------------------
trace = 0.0;

do i=1, NDIM 
     trace = trace + matrixc(i,i)
enddo
!! -----------------------------------------------------
!! END TRACE BLOCK
!! ----------------------------------------------------



!! -----------------------------------------------------
!! RESULTS AND DEALLOCATION BLOCK
!! -----------------------------------------------------
print *, NDIM, trace, cpu_end-cpu_start, wall_end-wall_start, mflops

if (allocated(matrixa)) deallocate(matrixa)
if (allocated(matrixb)) deallocate(matrixb)
if (allocated(matrixc)) deallocate(matrixc)
if (allocated(veca))    deallocate(veca)
if (allocated(vecb))    deallocate(vecb)
if (allocated(vecx))    deallocate(vecx)

call PAPIF_stop(eventSet, dp_ops, check);
if (check .ne. PAPI_OK) then
    print *, "Could Not Stop Papi Counters."
    call exit()
endif

call PAPIF_cleanup_eventset(eventSet, check);
if (check .ne. PAPI_OK) then
    print *, "Could Not Cleanup Papi Event Set."
    call exit()
endif

call PAPIF_destroy_eventset(eventSet, check);
if (check .ne. PAPI_OK) then
    print *, "Could Not Destroy Papi Event Set."
    call exit()
endif

!! -----------------------------------------------------
!! END RESULTS AND DEALLOCATION BLOCK
!! ----------------------------------------------------

enddo
!! END ITERATION LOOP

end program driver 
 
