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
         function pinv_rel(A, method, tol) result(Apinv)

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
                case('getrf')
                  Apinv = pinvLU_rel(A)
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


         !===============================================================================
         !> author: ZUO Zhihua
         !> Calculates the pseudoinverse of a matrix A using the LU decomposition.
         function pinvLU_rel(A) result(Apinv)

            ! Inputs:
            real(rk), dimension(:, :), intent(in)              :: A     ! Input matrix A

            ! Outputs:
            real(rk), dimension(size(A,2), size(A,1))          :: Apinv ! Pseudoinverse of A

            ! Local variables

            if (size(A, 1) == size(A, 2)) then
               Apinv = invLU_rel(A)
            elseif (size(A, 1) > size(A, 2)) then
               Apinv = transpose(A)
               Apinv = gemm(invLU_rel(gemm(Apinv, A)), Apinv)
            else
               Apinv = transpose(A)
               Apinv = gemm(Apinv, invLU_rel(gemm(A, Apinv)))
            end if

         end function pinvLU_rel
         !===============================================================================


         !===============================================================================
         !> author: ZUO Zhihua
         !> Calculates the inverse of a matrix A using the LU decomposition.
         function invLU_rel(A) result(Ainv)

            ! Inputs:
            real(rk), dimension(:, :), intent(in)              :: A     ! Input matrix A

            ! Outputs:
            real(rk), dimension(size(A,2), size(A,1))          :: Ainv  ! Inverse of A

            ! Local variables
            integer                                            :: ipiv(size(A, 1)), info
            real(rk)                                           :: work(size(A, 2))

            ! External subroutine for calculating the inverse of a matrix A using the LU decomposition.
            interface
               subroutine dgetrf(m, n, a, lda, ipiv, info)
                  import rk
                  integer, intent(in) :: m
                  integer, intent(in) :: n
                  integer, intent(in) :: lda
                  integer, intent(out) :: ipiv(*)
                  integer, intent(out) :: info
                  real(rk), intent(inout) :: a(lda, *)
               end subroutine dgetrf
               subroutine dgetri(n, a, lda, ipiv, work, lwork, info)
                  import rk
                  integer, intent(in) :: n
                  integer, intent(in) :: lda
                  integer, intent(in) :: lwork
                  integer, intent(out) :: ipiv(*)
                  integer, intent(out) :: info
                  real(rk), intent(inout) :: a(lda, *)
                  real(rk), intent(out) :: work(*)
               end subroutine dgetri
            end interface

            Ainv = A
            call dgetrf(size(A, 1), size(A, 2), Ainv, size(A, 1), ipiv, info)
            call dgetri(size(A, 2), Ainv, size(A, 1), ipiv, work, size(A, 2), info)

         end function invLU_rel
         !===============================================================================


         !===============================================================================
         !> author: ZUO Zhihua
         !> Calculates the matrix-matrix product of two matrices A and B.
         pure function gemm(A, B) result(C)

            ! Inputs:
            real(rk), dimension(:, :), intent(in)              :: A     ! Input matrix A
            real(rk), dimension(:, :), intent(in)              :: B     ! Input matrix B

            ! Outputs:
            real(rk), dimension(size(A,1), size(B,2))          :: C     ! Matrix-matrix product of A and B

            ! Local variables
            integer                                            :: m, n, k

            ! External subroutine for calculating the matrix-matrix product of two matrices A and B.
            interface
               pure subroutine dgemm(transa, transb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
                  import rk
                  integer, intent(in) :: ldc
                  integer, intent(in) :: ldb
                  integer, intent(in) :: lda
                  character, intent(in) :: transa
                  character, intent(in) :: transb
                  integer, intent(in) :: m
                  integer, intent(in) :: n
                  integer, intent(in) :: k
                  real(rk), intent(in) :: alpha
                  real(rk), intent(in) :: a(lda, *)
                  real(rk), intent(in) :: b(ldb, *)
                  real(rk), intent(in) :: beta
                  real(rk), intent(inout) :: c(ldc, *)
               end subroutine dgemm
            end interface

            m = size(A, 1)
            n = size(B, 2)
            k = size(A, 2)

            call dgemm('N', 'N', m, n, k, 1.0_rk, A, m, B, k, 0.0_rk, C, m)

         end function gemm

      end module pinverse
