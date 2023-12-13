program test1

   ! This Fortran test code demonstrates the usage of the inv function to calculate the matrix inverse&
   ! and verifies the results by comparing them with expected values obtained from MATLAB.

   use kinds                ! Import the module 'kinds' for precision types
   use forinv, only: inv    ! Import only the 'inv' function from the 'forinv' module

   implicit none

   real(rk), dimension(4,3) :: A           ! Define a 4x3 matrix A
   real(rk), dimension(3,4) :: Ainv, AinvM ! Define matrices Ainv and AinvM
   real(rk)                 :: rel_err

   !===============================================================================
   ! Initialize matrix A with values
   A(1,1) = 0.814723686393179_rk;  A(1,2) = 0.632359246225410_rk ;  A(1,3) = 0.957506835434298_rk
   A(2,1) = 0.905791937075619_rk;  A(2,2) = 0.0975404049994095_rk;  A(2,3) = 0.964888535199277_rk
   A(3,1) = 0.126986816293506_rk;  A(3,2) = 0.278498218867048_rk ;  A(3,3) = 0.157613081677548_rk
   A(4,1) = 0.913375856139019_rk;  A(4,2) = 0.546881519204984_rk ;  A(4,3) = 0.970592781760616_rk
   !===============================================================================


   !===============================================================================
   ! Define expected matrix AinvM=inv(A) from MATLAB results
   AinvM(1,1) = -11.8748418966254_rk 
   AinvM(1,2) = -1.88508351433968_rk
   AinvM(1,3) = 1.37910840173057_rk 
   AinvM(1,4) = 13.3647936345720_rk
   AinvM(2,1) = -0.232262684677810_rk
   AinvM(2,2) = -1.81026044323168_rk
   AinvM(2,3) = 1.12793815267075_rk 
   AinvM(2,4) = 1.84558847033752_rk
   AinvM(3,1) = 11.2481704396877_rk  
   AinvM(3,2) = 2.87726069020400_rk 
   AinvM(3,3) = -1.70926059704099_rk
   AinvM(3,4) = -12.6490061903496_rk
   !===============================================================================


   !===============================================================================
   ! Calculate the matrix inverse of A using the 'inv' function
   Ainv = inv(A)
   !===============================================================================


   !===============================================================================
   ! Calculate the relative error between MATLAB and Fortran results
   rel_err = norm2(Ainv - AinvM)/norm2(AinvM)

   print *, "Relative error:", rel_err
   if (rel_err < 1e-13_rk) then
      print *, "Test 1 passed."
      print*,""
   else
      print *, "Test 1 failed!"
      print*,""
   end if
   !===============================================================================

end program test1
