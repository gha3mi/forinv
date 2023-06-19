program test2

  ! This Fortran test code demonstrates the usage of the pinv function to calculate the matrix inverse
  ! of a randomly generated matrix. It also measures the execution time using the M_stopwatch module.

  use :: kinds                 ! Import the module 'kinds' for precision types
  use :: pinverse, only: pinv  ! Import only the 'pinv' function from the 'pinverse' module
  use :: M_stopwatch           ! Import the M_stopwatch module for timing measurements

  implicit none

  real(rk), dimension(:,:), allocatable :: A, Ainv  ! Define dynamically allocated matrices A and Ainv
  integer                               :: m, n     ! Define variables for matrix dimensions
  type(watchtype)                       :: w        ! Define a watchtype object for timing measurements

  m = 10000                  ! Set the number of rows for matrix A
  n = 1000                   ! Set the number of columns for matrix A
  allocate(A(m,n))           ! Allocate memory for matrix A
  call random_number(A)      ! Fill matrix A with random numbers between 0 and 1

  call create_watch(w)       ! Create a stopwatch object
  call start_watch(w)        ! Start the stopwatch

  Ainv = pinv(A)             ! Calculate the matrix inverse of A using the 'pinv' function

  call stop_watch(w)         ! Stop the stopwatch
  call print_watch(w)        ! Print the elapsed time measured by the stopwatch
  call destroy_watch(w)      ! Destroy the stopwatch object

  deallocate(Ainv)           ! Deallocate memory for matrix Ainv

  print*, "test2 completed." ! Print a message indicating the completion of the program

end program test2
