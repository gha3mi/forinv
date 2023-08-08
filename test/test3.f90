program test3

  ! This Fortran test code demonstrates the usage of the inv function to calculate the matrix inverse
  ! of a randomly generated matrix.

  use kinds                 ! Import the module 'kinds' for precision types
  use forinv, only: inv     ! Import only the 'inv' function from the 'forinv' module
  use fortime

  implicit none

  real(rk), dimension(:,:), allocatable :: A, Ainv  ! Define dynamically allocated matrices A and Ainv
  integer                               :: m, n     ! Define variables for matrix dimensions
  type(timer)                           :: w        ! Define a watchtype object for timing measurements

  m = 2000                    ! Set the number of rows for matrix A
  n = 200                     ! Set the number of columns for matrix A
  allocate(A(m,n))            ! Allocate memory for matrix A
  call random_number(A)       ! Fill matrix A with random numbers between 0 and 1

  call timer_start(w)

  Ainv = inv(A*10)          ! Calculate the matrix inverse of A using the 'inv' function

  call timer_stop(w,message=' Elapsed time (2000*200 , gesvd):')

  call timer_start(w)

  Ainv = inv(A*10, method='getrf')          ! Calculate the matrix inverse of A using the 'inv' function and the getrf method

  call timer_stop(w,message=' Elapsed time (2000*200 , getrf):')

  deallocate(Ainv)            ! Deallocate memory for matrix Ainv
  deallocate(A)               ! Deallocate memory for matrix A

  m = 2000                    ! Set the number of rows for matrix A
  n = 1800                    ! Set the number of columns for matrix A
  allocate(A(m,n))            ! Allocate memory for matrix A
  call random_number(A)       ! Fill matrix A with random numbers between 0 and 1

  call timer_start(w)

  Ainv = inv(A*10)           ! Calculate the matrix inverse of A using the 'inv' function

  call timer_stop(w,message=' Elapsed time (2000*1800, gesvd):')

  call timer_start(w)

  Ainv = inv(A*10, method='getrf')          ! Calculate the matrix inverse of A using the 'inv' function and the getrf method

  call timer_stop(w,message=' Elapsed time (2000*1800, getrf):')

  deallocate(Ainv)            ! Deallocate memory for matrix Ainv
  deallocate(A)               ! Deallocate memory for matrix A

  print *, "Test 3 passed."

end program test3
