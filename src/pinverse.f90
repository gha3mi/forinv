module pinverse

   ! This module provides functions and subroutines
   ! for calculating the singular value decomposition (SVD) and pseudoinverse of a matrix.

   use :: kinds  ! Import the module 'kinds' for precision types

   implicit none
   
   private              ! Make all module variables and procedures private

   public :: svd, pinv  ! Make svd and pinv functions public

   !===============================================================================
   interface svd
      procedure :: svd_rel ! Interface for the svd_rel subroutine
   end interface
   !===============================================================================


   !===============================================================================
   interface pinv
      procedure :: pinverse_rel ! Interface for the pinverse_rel function
   end interface
   !===============================================================================

contains

   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the singular value decomposition (SVD) of a matrix A.
   pure subroutine svd_rel(A, U,S,VT)

      ! Inputs:
      real(rk), dimension(:, :),              intent(in)  :: A    ! Input matrix A  

      ! Outputs:
      real(rk), dimension(:, :), allocatable, intent(out) :: U    ! Left singular vectors
      real(rk), dimension(:, :), allocatable, intent(out) :: VT   ! Right singular vectors
      real(rk), dimension(:),    allocatable, intent(out) :: S    ! Singular values

      ! Local variables
      real(rk), dimension(:),    allocatable              :: work ! Work array
      integer                                             :: m, n, minnm, lwork, info, i, j

      ! External subroutine for calculating the SVD
      interface dgesvd
         pure subroutine dgesvd(jobuf,jobvtf,mf,nf,af,ldaf,sf,uf,lduf,vtf,ldvtf,workf,lworkf,infof)
            use kinds
            character, intent(in)  :: jobuf, jobvtf
            integer,   intent(in)  :: mf, nf, ldaf, lduf, ldvtf, lworkf
            real(rk),  intent(in)  :: Af(ldaf, *)
            real(rk),  intent(out) :: Sf(min(mf, nf))
            real(rk),  intent(out) :: Uf(lduf, *), VTf(ldvtf, *)
            real(rk),  intent(out) :: workf(*)
            integer,   intent(out) :: infof
         end subroutine dgesvd
      end interface

      m = size(A, 1)
      n = size(A, 2)

      minnm = max(1,5*min(n,m))

      allocate(U(m, m), VT(n, n), S(minnm), work(5 * minnm))

      lwork = max(1,3*min(n,m) + max(n,m),5*min(n,m))

      call dgesvd('A', 'A', m, n, A, m, S, U, m, VT, n, work, lwork, info)

      deallocate(work)
   end subroutine svd_rel
   !===============================================================================


   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the pseudoinverse of a matrix A using the SVD.
   pure function pinverse_rel(A) result(Apinv)

      ! Inputs:
      real(rk), dimension(:, :), contiguous, intent(in)  :: A     ! Input matrix A

      ! Outputs:
      real(rk), dimension(:, :), allocatable             :: Apinv ! Pseudoinverse of A
      
      ! Local variables
      real(rk), dimension(:, :), allocatable             :: U     ! Left singular vectors
      real(rk), dimension(:, :), allocatable             :: VT    ! Right singular vectors
      real(rk), dimension(:),    allocatable             :: S     ! Singular values
      integer                                            :: m, n, i, j, irank

      m = size(A, 1)
      n = size(A, 2)

      allocate(Apinv(n, m))

      call svd_rel(A, U,S,VT)

      Apinv = 0.0_rk

      do irank = 1, min(n,m)
         do j = 1, m
            do i = 1, n
               Apinv(i, j) = Apinv(i, j) + VT(irank, i) * U(j, irank) / S(irank)
            end do
         end do
      end do

   end function pinverse_rel
   !===============================================================================

end module pinverse
