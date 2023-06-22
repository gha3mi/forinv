program test3

  ! This Fortran test code demonstrates the usage of the pinv function to calculate the matrix inverse
  ! of a randomly generated matrix.

  use :: kinds                 ! Import the module 'kinds' for precision types
  use :: pinverse, only: pinv  ! Import only the 'pinv' function from the 'pinverse' module
  use :: fortime

  implicit none

  real(rk), dimension(:,:), allocatable :: A, Ainv  ! Define dynamically allocated matrices A and Ainv
  integer                               :: m, n     ! Define variables for matrix dimensions
  type(timer)                           :: w        ! Define a watchtype object for timing measurements

  m = 10000                   ! Set the number of rows for matrix A
  n = 1000                    ! Set the number of columns for matrix A
  allocate(A(m,n))            ! Allocate memory for matrix A
  call random_number(A)       ! Fill matrix A with random numbers between 0 and 1

  call timer_start(w)

  Ainv = pinv(A*10)          ! Calculate the matrix inverse of A using the 'pinv' function

  call timer_stop(w,message=' Elapsed time:')

  deallocate(Ainv)           ! Deallocate memory for matrix Ainv

  print *, "Test 3 passed."

end program test3
