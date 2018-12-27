! Copyright (c) 2017-2018 Etienne Descamps
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification,
! are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation and/or
!    other materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
! ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
! WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
! IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
! LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
! EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

Module mlf_matrix
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_utils
  Use mlf_rand
  IMPLICIT NONE
  PRIVATE
  Public :: MatrixOuter, SymMatrixEigenDecomposition, InverseSymMatrix, SqrtSymMatrix
  Public :: SymmetrizeMatrix, SymmetrizeMatrixUp, TensorAddVectDown, IdentityMatrix
  Public :: GECovMatrix, randOrthogonal, SymMatrixSelectEigenValues
Contains
  
  ! Easier to use that dsyr or dsyr2
  Function MatrixOuter(V, W) Result(M)
    real(c_double), intent(in) :: V(:)
    real(c_double), intent(in), optional :: W(:)
    real(c_double) :: M(SIZE(V, 1),SIZE(V, 1))
    integer :: N
    N = SIZE(V, 1)
    If(PRESENT(W)) then
      M = MATMUL(RESHAPE(V, (/ N,1 /)), RESHAPE(W, (/ 1,N /)))
    Else
      M = MATMUL(RESHAPE(V, (/ N,1 /)), RESHAPE(V, (/ 1,N /)))
    Endif
  End Function MatrixOuter

  ! Symmetrize matrix (using the lower part) ('L' to full)
  Subroutine SymmetrizeMatrix(M)
    real(c_double), intent(inout) :: M(:,:)
    integer(c_int) :: i, j, N
    N = size(M,1)
    Do i = 1,N
      FORALL (j=(i+1):N) M(i,j) = M(j,i)
    End Do
  End Subroutine SymmetrizeMatrix

  ! Symmetrize matrix (using the upper part) ('U' to full)
  Subroutine SymmetrizeMatrixUp(M)
    real(c_double), intent(inout) :: M(:,:)
    integer(c_int) :: i, j, N
    N = size(M,1)
    Do i = 1,N
      FORALL (j=(i+1):N) M(j,i) = M(i,j)
    End Do
  End Subroutine SymmetrizeMatrixUp
   
  ! Procedure that proceed add to a matrix the matrix product of vect*transpose(vect)
  ! (It is faster to use this method and symmetrized the matrix thereafter)
  Subroutine TensorAddVectDown(N, V, M, w)
    integer(c_int), intent(in) :: N
    integer(c_int) :: i, j
    real(c_double), intent(in) :: V(N), w
    real(c_double), intent(inout) :: M(N,N)
    real(c_double) :: x
    Do i=1,N
      x = w*V(i)
      FORALL (j=i:N) M(j,i) = M(j,i) + x*V(j)
    End Do
  End Subroutine TensorAddVectDown

  ! Utility function for initializing a matrix as the Identity
  Subroutine IdentityMatrix(M)
    real(c_double) :: M(:,:)
    integer(c_int) :: i
    M = 0
    FORALL (i= 1:SIZE(M,1)) M(i,i) = 1
  End Subroutine IdentityMatrix

  ! Generate a covariance matrix and compute the mean of the individuals
  ! If W is provided, W shall be normalized (sum(W) == 1)
  Subroutine GECovMatrix(C, A, Mu, W, minW, Mu0)
    real(c_double), intent(in) :: A(:,:)
    real(c_double), intent(in), optional :: minW, W(:)
    real(c_double), intent(in), target, optional :: Mu0(:)
    real(c_double), intent(out) :: C(:,:)
    real(c_double), intent(out), target :: Mu(:)
    real(c_double), allocatable :: B(:,:)
    real(c_double), pointer :: MuX(:)
    real(c_double) :: minW0, ms
    integer :: i,j,N,M
    external :: dsyrk
    N = SIZE(A,1); M = SIZE(A,2)
    ALLOCATE(B(N,M))
    If(PRESENT(Mu0)) Then
      ! Mu0 is provided, so no Bessel's correction is required
      MuX => Mu0
      ms = 1.0
    Else
      ! Use Bessel's correction for getting an unbiased covariance
      If(PRESENT(W)) Then
        ms = 1d0/(1d0-SUM(W*W))
      Else
        ms = 1d0/REAL(M-1,kind=8)
      Endif
      MuX => Mu
    Endif
    If(PRESENT(W)) Then
      If(PRESENT(minW)) Then
        minW0 = minW
      Else
        ! Remove vector that are beyond the precision limit of the main vector
        minW0 = EPSILON(1d0)*MAXVAL(W)
      Endif
      Mu = MATMUL(A,W)
      j = 0
      Do i=1,M
        if(W(i)<minW0) CYCLE
        j = j+1
        B(:,j) = ms*SQRT(W(i))*(A(:,i)-MuX)
      End Do
    Else
      Mu = Mean(A,2)
      FORALL(i=1:M) B(:,i) = ms*(A(:,i)-MuX)
      j = M
    Endif
    ! Slowest part (in O(j*N²))
    CALL dsyrk('L', 'N', N, j, 1d0, B, N, 0d0, C, N)
    CALL SymmetrizeMatrix(C)
  End Subroutine GECovMatrix

  ! Random orthogonal matrix using the Haar distribution (Stewart 1980)
  ! Inspired by SciPy code (rvs method of _multivariate.py)
  ! (seems that Python is not as clear as FORTRAN for matrix computation)
  ! If Sigma is defined, returns O*diag(Sigma) (useful for MV random sampling)
  Subroutine randOrthogonal(H, Sigma)
    real(c_double), intent(out) :: H(:,:)
    real(c_double), intent(in), optional :: Sigma(:)
    integer :: N, k
    real(c_double), allocatable, target :: D(:), U(:), Mat(:,:)
    real(c_double), pointer :: x(:), Hx(:,:)
    N = SIZE(H, 1)
    ALLOCATE(D(N), U(N), Mat(N,N))
    CALL IdentityMatrix(H)
    D = 1d0
    Do k=1,N-1
      x => U(k:)
      Hx => Mat(k:,k:)
      CALL RandN(x)
      D(k) = SIGN(1d0, x(1))
      x(1) = x(1) - D(k)*NORM2(x)
      ! Householder matrix generated by x
      CALL IdentityMatrix(Hx)
      Hx = Hx- 2d0*MatrixOuter(x)/SUM(x*x)
      ! Most complex operation (O(N*k²))
      H(k:,:)=MATMUL(Hx,H(k:,:))
    End Do
    ! Fix last element of diagonal such as determinant = 1
    D(N) = PRODUCT(D)*(2*MODULO(N,2)-1)
    If(PRESENT(Sigma)) D= D*Sigma
    ! H = H*diag(H) (*diag(Sigma))
    FORALL (k=1:N) H(:,k) = H(:,k)*D(k)
  End Subroutine randOrthogonal
  
  Integer Function SymMatrixEigenDecomposition(C, LD, LB) result(info)
    ! Get orthogonal eigenvectors (LB) and eigenvalue (LD) from the symmetric matrix C
    ! (eigenvalues in ascending order)
    integer :: ND, LWORK
    real(c_double), intent(in) :: C(:,:)
    real(c_double), intent(out) :: LD(:), LB(:,:)
    real(c_double), allocatable :: LE(:), TAU(:), WORK(:,:)
    ND = SIZE(C,1)
    ALLOCATE(LE(ND-1), TAU(ND), WORK(ND,ND))
    LWORK = ND*ND
    LB = C
    ! Use same functions as dsyev
    ! Bidiagonal symmetric transformation
    CALL dsytrd('U', ND, LB, ND, LD, LE, TAU, WORK, LWORK, info)
    If(info /= 0) RETURN
    ! Generate orthogonal matrix using reflectors generated by dsytrd
    CALL dorgtr('U', ND, LB, ND, TAU, WORK, LWORK, info)
    If(info /= 0) RETURN
    ! Compute eigen values and orthogonal matrix
    CALL dsteqr('V', ND, LD, LE, LB, ND, WORK, info)
  End Function SymMatrixEigenDecomposition

  !  
  Integer Function SymMatrixSelectEigenValues(C, SLD, SLB) result(info)
    ! Get the most prominent eigenvectors (SLB) and eigenvalue (SLD) from the symmetric matrix C
    ! (eigenvalues in ascending order)
    real(c_double), intent(in) :: C(:,:)
    real(c_double), intent(out) :: SLD(:), SLB(:,:)
    real(c_double), allocatable :: LD(:), A(:,:), WORK(:,:)
    integer, allocatable :: IWORK(:), ISUPPZ(:)
    integer :: K, ND, nR
    ND = SIZE(C,1); nR = SIZE(SLD)
    ALLOCATE(LD(ND), A(ND,ND), WORK(ND,ND+6), IWORK(10*ND), ISUPPZ(2*nR))
    A = C
    ! Compute eigen values and orthogonal matrix
    CALL dsyevr('V', 'I', 'U', ND, A, ND, 0d0, HUGE(0d0), ND-nR+1, ND, 1d-12, K, LD, SLB, ND, iSUPPZ, &
      WORK, SIZE(WORK), IWORK, SIZE(IWORK), info)
    If(info /= 0) RETURN
    If(K < nR) Then
      info = -1
      RETURN
    Endif
    SLD = LD(1:nR)
  End Function SymMatrixSelectEigenValues

  ! Inverse covariance matrix (using pseudo inverse when singular)
  !   output lnDet = 0.5*log(determinant(2*pi*C))
  Integer Function InverseSymMatrix(C, invC, lnDet, doEigen) Result(info)
    real(c_double), intent(in) :: C(:,:)
    real(c_double), intent(out) :: invC(:,:)
    real(c_double), intent(out) :: lnDet
    logical, intent(in), optional :: doEigen
    logical :: eigen
    integer :: ND, i
    external :: dsytrd, dorgtr, dsteqr, dpotrf, dpotri
    ND = size(C,1)
    eigen = .FALSE.
    If(present(doEigen)) eigen = doEigen
    invC = C
    info = 1
    ! Compute Cholesky's decomposition (if eigen decomposition not requested)
    If(.NOT. eigen) CALL dpotrf('U', ND, invC, ND, info)
    If(info == 0) Then
      ! The fastest method for matrix inversion
      ! At this step invC contains the Cholesky decomposition
      ! So det(invC) = prod invC_ii = sqrt(det(C))
      ! ln(sqrt(det(2*pi*C))) = ND/2*ln(2*pi) + sum(ln(invC_ii))
      lnDet = 0.5d0*ND*log(2.0*mlf_PI)
      Do i=1,ND
        lnDet = lnDet + log(invC(i,i))
      End Do
      ! Compute inverse of C using inversion of the triangular Cholesky's decomposition
      CALL dpotri('U', ND, invC, ND, info)
      CALL SymmetrizeMatrixUp(invC)
    Else
      ! More stable method for inversing matrix
      BLOCK
        ! If the matrix is singular, we use the more stable eigenvalue decomposition
        ! Can occurs when ND is greater or equivalent to dataset size
        ! If lots of problem occurs in this case, you shall reduce the dimensionality (with e.g. PCA)
        real(c_double) :: LD(ND), LE(ND-1), TAU(ND), WORK(ND,ND), valM, penality
        integer :: LWORK
        LWORK = ND*ND
        ! Reinit matrix with C (if Cholesky's decomposition has been iterated)
        If(.NOT. eigen) invC = C
        ! Use same functions as dsyev
        ! Bidiagonal symmetric transformation
        CALL dsytrd('U', ND, invC, ND, LD, LE, TAU, WORK, LWORK, info)
        If(info /= 0) RETURN
        ! Generate orthogonal matrix using reflectors generated by dsytrd
        CALL dorgtr('U', ND, invC, ND, TAU, WORK, LWORK, info)
        If(info /= 0) RETURN
        ! Compute eigen values and orthogonal matrix
        CALL dsteqr('V', ND, LD, LE, invC, ND, WORK, info)
        If(info /= 0) RETURN
        ! invC contains the orthogonal matrix and LD the diagonal of the matrix

        ! Set a penality factor that will be used on null diagonal elements
        ! This factor sets probability of elements outside the hyper plane to a tiny value
        ! (useful in high dimension cases)
        ! We replace all diagonal value below valM = 10e-9*max(LD) by valM
        valM = 10e-9*MAXVAL(LD)
        penality = 1d0/valM ! Inverse of valM
        WHERE (LD <= valM)
          LD = penality
        ELSEWHERE
          LD = 1d0/LD
        END WHERE
        lnDet = 0.5d0*(ND*LOG(2.0*mlf_PI)+SUM(LOG(LD)))
        ! For a symmetric matrix, the pseudo determinant is the product of non null eigenvalues
        ! We compute the pseudo-inverse matrix with a penality value on element that are nul
        ! invC = O * diag(pseudoInv(D)) * O**T
        Do i=1,ND
          work(:,i) = invC(:,i)*LD(i)
        End Do
        invC = MATMUL(work, TRANSPOSE(invC))
      END BLOCK
    Endif
  End Function InverseSymMatrix

  Integer Function SqrtSymMatrix(C, C12) Result(info)
    integer :: ND, i
    real(c_double), intent(in) :: C(:,:)
    real(c_double), intent(out) :: C12(:,:)
    real(c_double), allocatable :: LD(:)
    ND = size(C,1)
    ALLOCATE(LD(ND))
    info = SymMatrixEigenDecomposition(C, LD, C12)
    If(info /= 0) RETURN
    LD = sqrt(LD)
    FORALL (i=1:ND) C12(:, i) = C12(:, i)*LD(i)
  End Function SqrtSymMatrix
End Module mlf_matrix
