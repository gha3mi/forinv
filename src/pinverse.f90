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
#if defined (PURE)
   pure function pinv_rel(A, method, tol) result(Apinv)
#elif defined (IMPURE)
   impure function pinv_rel(A, method, tol) result(Apinv)
#endif

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
#if defined (PURE)
   pure function pinvSVD_rel(A, method, tol) result(Apinv)
#elif defined (IMPURE)
   impure function pinvSVD_rel(A, method, tol) result(Apinv)
#endif

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
   !> author: ZUO Zhihua, Seyed Ali Ghasemi
   !> Calculates the pseudoinverse of a matrix A using the LU decomposition.
#if defined (PURE)
   pure function pinvLU_rel(A) result(Apinv)
#elif defined (IMPURE)
   impure function pinvLU_rel(A) result(Apinv)
#endif

      ! Inputs:
      real(rk), dimension(:, :), contiguous, intent(in) :: A     ! Input matrix A

      ! Outputs:
      real(rk), dimension(size(A,2), size(A,1))         :: Apinv ! Pseudoinverse of A

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
   !> author: ZUO Zhihua, Seyed Ali Ghasemi
   !> Calculates the inverse of a matrix A using the LU decomposition.
#if defined (PURE)
   pure function invLU_rel(A) result(Ainv)
#elif defined (IMPURE)
   impure function invLU_rel(A) result(Ainv)
#endif

      ! Inputs:
      real(rk), dimension(:, :), contiguous, intent(in) :: A     ! Input matrix A

      ! Outputs:
      real(rk), dimension(size(A,2), size(A,1))         :: Ainv  ! Inverse of A

      ! Local variables
      integer                                            :: ipiv(size(A, 1)), info
      real(rk)                                           :: work(size(A, 2))

      ! External subroutine for calculating the inverse of a matrix A using the LU decomposition.
      interface
#if defined (PURE)
         pure subroutine dgetrf(f_m, f_n, f_a, f_lda, f_ipiv, f_info)
#elif defined (IMPURE)
         impure subroutine dgetrf(f_m, f_n, f_a, f_lda, f_ipiv, f_info)
#endif
            import rk
            integer,  intent(in)    :: f_m
            integer,  intent(in)    :: f_n
            integer,  intent(in)    :: f_lda
            integer,  intent(out)   :: f_ipiv(*)
            integer,  intent(out)   :: f_info
            real(rk), intent(inout) :: f_a(f_lda, *)
         end subroutine dgetrf
#if defined (PURE)
         pure subroutine dgetri(f_n, f_a, f_lda, f_ipiv, f_work, f_lwork, f_info)
#elif defined (IMPURE)
         impure subroutine dgetri(f_n, f_a, f_lda, f_ipiv, f_work, f_lwork, f_info)
#endif
            import rk
            integer,  intent(in)    :: f_n
            integer,  intent(in)    :: f_lda
            integer,  intent(in)    :: f_lwork
            integer,  intent(out)   :: f_ipiv(*)
            integer,  intent(out)   :: f_info
            real(rk), intent(inout) :: f_a(f_lda, *)
            real(rk), intent(out)   :: f_work(*)
         end subroutine dgetri
      end interface

      Ainv = A
      call dgetrf(size(A, 1), size(A, 2), Ainv, size(A, 1), ipiv, info)
      call dgetri(size(A, 2), Ainv, size(A, 1), ipiv, work, size(A, 2), info)

   end function invLU_rel
   !===============================================================================


   !===============================================================================
   !> author: ZUO Zhihua, Seyed Ali Ghasemi
   !> Calculates the matrix-matrix product of two matrices A and B.
#if defined (PURE)
   pure function gemm(A, B) result(C)
#elif defined (IMPURE)
   impure function gemm(A, B) result(C)
#endif

      ! Inputs:
      real(rk), dimension(:, :), contiguous, intent(in) :: A     ! Input matrix A
      real(rk), dimension(:, :), contiguous, intent(in) :: B     ! Input matrix B

      ! Outputs:
      real(rk), dimension(size(A,1), size(B,2))          :: C     ! Matrix-matrix product of A and B

      ! Local variables
      integer                                            :: m, n, k

      ! External subroutine for calculating the matrix-matrix product of two matrices A and B.
      interface
#if defined (PURE)
         pure subroutine dgemm(f_transa, f_transb, f_m, f_n, f_k, f_alpha, f_a, f_lda, f_b, f_ldb, f_beta, f_c, f_ldc)
#elif defined (IMPURE)
         impure subroutine dgemm(f_transa, f_transb, f_m, f_n, f_k, f_alpha, f_a, f_lda, f_b, f_ldb, f_beta, f_c, f_ldc)
#endif
            import rk
            integer,   intent(in)    :: f_ldc
            integer,   intent(in)    :: f_ldb
            integer,   intent(in)    :: f_lda
            character, intent(in)    :: f_transa
            character, intent(in)    :: f_transb
            integer,   intent(in)    :: f_m
            integer,   intent(in)    :: f_n
            integer,   intent(in)    :: f_k
            real(rk),  intent(in)    :: f_alpha
            real(rk),  intent(in)    :: f_a(f_lda, *)
            real(rk),  intent(in)    :: f_b(f_ldb, *)
            real(rk),  intent(in)    :: f_beta
            real(rk),  intent(inout) :: f_c(f_ldc, *)
         end subroutine dgemm
      end interface

      m = size(A, 1)
      n = size(B, 2)
      k = size(A, 2)

      call dgemm('N', 'N', m, n, k, 1.0_rk, A, m, B, k, 0.0_rk, C, m)

   end function gemm
   !===============================================================================

end module pinverse
