module pinverse

   !This module provides functions and subroutines for pseudoinverse calculations.

   use :: kinds
   use :: forsvd, only: svd

   implicit none
   
   private

   public ::  pinv

   !===============================================================================
   interface pinv
      procedure :: pinverse_rel ! Interface for the pinverse_rel function
   end interface
   !===============================================================================

contains

   !===============================================================================
   !> author: Seyed Ali Ghasemi
   !> Calculates the pseudoinverse of a matrix A using the SVD.
   pure function pinverse_rel(A, tol) result(Apinv)

      ! Inputs:
      real(rk), dimension(:, :), contiguous, intent(in)  :: A     ! Input matrix A
      real(rk), intent(in), optional                     :: tol

      ! Outputs:
      real(rk), dimension(size(A,2), size(A,1))          :: Apinv ! Pseudoinverse of A

      ! Local variables
      real(rk), dimension(size(A,1), size(A,1))          :: U    ! Left singular vectors
      real(rk), dimension(size(A,2), size(A,2))          :: VT   ! Right singular vectors
      real(rk), dimension(min(size(A,1), size(A,2)))     :: S    ! Singular values
      integer                                            :: m, n, i, j, irank, rank

      m = size(A, 1)
      n = size(A, 2)

      call svd(A, U,S,VT)

      if (.not. present(tol)) then
         rank = min(m,n)
      else
         rank = count(S > tol)
      end if

      Apinv = 0.0_rk

      do irank = 1, rank
         do j = 1, m
            do i = 1, n
               Apinv(i, j) = Apinv(i, j) + VT(irank, i) * U(j, irank) / S(irank)
            end do
         end do
      end do

   end function pinverse_rel
   !===============================================================================

end module pinverse