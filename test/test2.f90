program test2

  use :: kinds                 ! Import the module 'kinds' for precision types
  use :: pinverse, only: pinv  ! Import only the 'pinv' function from the 'pinverse' module

  implicit none

  real(rk), dimension(:,:), allocatable :: A, Ainv1, Ainv2  ! Define dynamically allocated matrices A and Ainv
  integer                               :: m, n             ! Define variables for matrix dimensions
  real(rk)                              :: rel_err

  m = 1000                   ! Set the number of rows for matrix A
  n = 100                    ! Set the number of columns for matrix A
  allocate(A(m,n))           ! Allocate memory for matrix A
  call random_number(A)      ! Fill matrix A with random numbers between 0 and 1

  Ainv1 = pinv(A*10)         ! Calculate the matrix inverse of A using the 'pinv' function
  Ainv2 = pinv(A*10)         ! Calculate the matrix inverse of A using the 'pinv' function

  rel_err = norm2(Ainv1 - Ainv2)/norm2(Ainv1)

  print*,"Ainv1(1,1) =", Ainv1(1,1)
  print*,"Ainv2(1,1) =", Ainv2(1,1)
  print *, "Relative error:", rel_err
  if (rel_err < 1e-13_rk) then
     print *, "Test 2 passed."
     print*,""
  else
     print *, "Test 2 failed!"
     print*,""
  end if

  deallocate(Ainv1,Ainv2)           ! Deallocate memory for matrix Ainv

end program test2
