      module pinverse

         !This module provides functions and subroutines for pseudoinverse calculations.

         use :: kinds
         use :: forsvd, only: svd

         implicit none

         private

         public :: pinv

         !===============================================================================
         interface pinv
            procedure :: pinv_rel
         end interface
         !===============================================================================

      contains

         !===============================================================================
         !> author: Seyed Ali Ghasemi
         !> Calculates the pseudoinverse of a matrix A.
         pure function pinv_rel(A, method, tol) result(Apinv)

            ! Inputs:
            real(rk),     dimension(:, :), contiguous, intent(in) :: A
            real(rk),     intent(in), optional                    :: tol
            character(*), intent(in), optional                    :: method

            ! Outputs:
            real(rk), dimension(size(A,2), size(A,1))          :: Apinv ! Pseudoinverse of A


            if (present(method)) then
               select case (method)
                case('gesvd','gesdd')
                  Apinv = pinvSVD_rel(A, method, tol)
                case default
                  error stop 'method is not valid.'
               end select
            else
               Apinv = pinvSVD_rel(A, tol=tol)
            end if

         end function pinv_rel
         !===============================================================================


         !===============================================================================
         !> author: Seyed Ali Ghasemi
         !> Calculates the pseudoinverse of a matrix A using the SVD.
         pure function pinvSVD_rel(A, method, tol) result(Apinv)

            ! Inputs:
            real(rk), dimension(:, :), contiguous, intent(in)  :: A     ! Input matrix A
            real(rk),     intent(in), optional                 :: tol
            character(*), intent(in), optional                 :: method

            ! Outputs:
            real(rk), dimension(size(A,2), size(A,1))          :: Apinv ! Pseudoinverse of A

            ! Local variables
            real(rk), dimension(size(A,1), size(A,1))          :: U    ! Left singular vectors
            real(rk), dimension(size(A,2), size(A,2))          :: VT   ! Right singular vectors
            real(rk), dimension(min(size(A,1), size(A,2)))     :: S    ! Singular values
            integer                                            :: m, n, i, j, irank, rank

            m = size(A, 1)
            n = size(A, 2)

            call svd(A, U,S,VT, method)

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

         end function pinvSVD_rel
         !===============================================================================

      end module pinverse
