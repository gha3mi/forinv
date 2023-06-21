program test3

  use :: kinds                 ! Import the module 'kinds' for precision types
  use :: pinverse, only: pinv  ! Import only the 'pinv' function from the 'pinverse' module
  use :: M_stopwatch           ! Import the M_stopwatch module for timing measurements

  implicit none

  real(rk), dimension(:,:), allocatable :: A, Ainv1, Ainv2  ! Define dynamically allocated matrices A and Ainv
  integer                               :: m, n             ! Define variables for matrix dimensions
  type(watchtype)                       :: w                ! Define a watchtype object for timing measurements

  m = 1000                   ! Set the number of rows for matrix A
  n = 100                    ! Set the number of columns for matrix A
  allocate(A(m,n))           ! Allocate memory for matrix A
  call random_number(A)      ! Fill matrix A with random numbers between 0 and 1

  call create_watch(w)       ! Create a stopwatch object
  call start_watch(w)        ! Start the stopwatch

  Ainv1 = pinv(A*10)         ! Calculate the matrix inverse of A using the 'pinv' function

  call stop_watch(w)         ! Stop the stopwatch
  call print_watch(w)        ! Print the elapsed time measured by the stopwatch
  call destroy_watch(w)      ! Destroy the stopwatch object

  Ainv2 = pinv(A*10)             ! Calculate the matrix inverse of A using the 'pinv' function

  print*,"Ainv1(1,1) =", Ainv1(1,1)
  print*,"Ainv2(1,1) =", Ainv2(1,1)
  print*,"relative error:", norm2(Ainv1 - Ainv2)/norm2(Ainv1)

  deallocate(Ainv1,Ainv2)           ! Deallocate memory for matrix Ainv

  print*, "test3 completed." ! Print a message indicating the completion of the program
  print*, "--------------------------------------------------------------"

end program test3
