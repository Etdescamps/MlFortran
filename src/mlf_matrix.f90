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
  Public :: GECovMatrix, randOrthogonal
Contains
  
  ! Easier to use that dsyr or dsyr2
  function MatrixOuter(V, W) result(M)
    real(c_double), intent(in) :: V(:)
    real(c_double), intent(in), optional :: W(:)
    real(c_double) :: M(size(V, 1),size(V, 1))
    integer :: N
    N = size(V, 1)
    if(present(W)) then
      M = matmul(reshape(V, (/ N,1 /)), reshape(W, (/ 1,N /)))
    else
      M = matmul(reshape(V, (/ N,1 /)), reshape(V, (/ 1,N /)))
    endif
  end function MatrixOuter

  ! Symmetrize matrix (using the lower part) ('L' to full)
  subroutine SymmetrizeMatrix(M)
    real(c_double), intent(inout) :: M(:,:)
    integer(c_int) :: i, j, N
    N = size(M,1)
    do i = 1,N
      forall (j=(i+1):N) M(i,j) = M(j,i)
    end do
  end subroutine SymmetrizeMatrix

  ! Symmetrize matrix (using the upper part) ('U' to full)
  subroutine SymmetrizeMatrixUp(M)
    real(c_double), intent(inout) :: M(:,:)
    integer(c_int) :: i, j, N
    N = size(M,1)
    do i = 1,N
      forall (j=(i+1):N) M(j,i) = M(i,j)
    end do
  end subroutine SymmetrizeMatrixUp
   
  ! Procedure that proceed add to a matrix the matrix product of vect*transpose(vect)
  ! (It is faster to use this method and symmetrized the matrix thereafter)
  subroutine TensorAddVectDown(N, V, M, w)
    integer(c_int), intent(in) :: N
    integer(c_int) :: i, j
    real(c_double), intent(in) :: V(N), w
    real(c_double), intent(inout) :: M(N,N)
    real(c_double) :: x
    do i=1,N
      x = w*V(i)
      forall (j=i:N) M(j,i) = M(j,i) + x*V(j)
    end do
  end subroutine TensorAddVectDown

  ! Utility function for initializing a matrix as the Identity
  subroutine IdentityMatrix(M)
    real(c_double) :: M(:,:)
    integer(c_int) :: i
    M = 0
    forall (i= 1:size(M,1)) M(i,i) = 1
  end subroutine IdentityMatrix

  ! Generate a covariance matrix and compute the mean of the individuals
  ! If W is provided, W shall be normalized (sum(W) == 1)
  subroutine GECovMatrix(C, A, Mu, W, minW, Mu0)
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
    N = size(A,1); M = size(A,2)
    allocate(B(N,M))
    if(present(Mu0)) then
      ! Mu0 is provided, so no Bessel's correction is required
      MuX => Mu0
      ms = 1.0
    else
      ! Use Bessel's correction for getting an unbiased covariance
      if(present(W)) then
        ms = 1d0/(1d0-sum(W*W))
      else
        ms = 1d0/real(M-1,kind=8)
      endif
      MuX => Mu
    endif
    if(present(W)) then
      if(present(minW)) then
        minW0 = minW
      else
        ! Remove vector that are beyond the precision limit of the main vector
        minW0 = epsilon(1d0)*MAXVAL(W)
      endif
      Mu = matmul(A,W)
      j = 0
      do i=1,M
        if(W(i)<minW0) CYCLE
        j = j+1
        B(:,j) = ms*sqrt(W(i))*(A(:,i)-MuX)
      end do
    else
      Mu = Mean(A,2)
      forall(i=1:M) B(:,i) = ms*(A(:,i)-MuX)
      j = M
    endif
    ! Slowest part (in O(j*N²))
    call dsyrk('L', 'N', N, j, 1d0, B, N, 0d0, C, N)
    call SymmetrizeMatrix(C)
  end subroutine GECovMatrix


  ! Random orthogonal matrix using the Haar distribution (Stewart 1980)
  ! Inspired by SciPy code (rvs method of _multivariate.py)
  ! (seems that Python is not as clear as FORTRAN for matrix computation)
  ! If Sigma is defined, returns O*diag(Sigma) (useful for MV random sampling)
  subroutine randOrthogonal(H, Sigma)
    real(c_double), intent(out) :: H(:,:)
    real(c_double), intent(in), optional :: Sigma(:)
    integer :: N, k
    real(c_double), allocatable, target :: D(:), U(:), Mat(:,:)
    real(c_double), pointer :: x(:), Hx(:,:)
    N = size(H, 1)
    allocate(D(N), U(N), Mat(N,N))
    call IdentityMatrix(H)
    D = 1d0
    do k=1,N-1
      x => U(k:)
      Hx => Mat(k:,k:)
      call RandN(x)
      D(k) = sign(1d0, x(1))
      x(1) = x(1) - D(k)*norm2(x)
      ! Householder matrix generated by x
      call IdentityMatrix(Hx)
      Hx = Hx- 2d0*MatrixOuter(x)/sum(x*x)
      ! Most complex operation (O(N*k²))
      H(k:,:)=matmul(Hx,H(k:,:))
    end do
    ! Fix last element of diagonal such as determinant = 1
    D(N) = product(D)*(2*modulo(N,2)-1)
    if(present(Sigma)) D= D*Sigma
    ! H = H*diag(H) (*diag(Sigma))
    forall (k=1:N) H(:,k) = H(:,k)*D(k)
  end subroutine randOrthogonal
  
  integer Function SymMatrixEigenDecomposition(C, LD, LB) result(info)
    ! Get orthogonal eigenvectors (LB) and eigenvalue (LD) of the symmetric matrix C
    ! (eigenvalues in ascending order)
    integer :: ND, LWORK
    real(c_double), intent(in) :: C(:,:)
    real(c_double), intent(out) :: LD(:), LB(:,:)
    real(c_double), allocatable :: LE(:), TAU(:), WORK(:,:)
    ND = size(C,1)
    allocate(LE(ND-1), TAU(ND), WORK(ND,ND))
    LWORK = ND*ND
    LB = C
    ! Use same functions as dsyev
    ! Bidiagonal symmetric transformation
    call dsytrd('U', ND, LB, ND, LD, LE, TAU, WORK, LWORK, info)
    if(info /= 0) RETURN
    ! Generate orthogonal matrix using reflectors generated by dsytrd
    call dorgtr('U', ND, LB, ND, TAU, WORK, LWORK, info)
    if(info /= 0) RETURN
    ! Compute eigen values and orthogonal matrix
    call dsteqr('V', ND, LD, LE, LB, ND, WORK, info)
  End Function SymMatrixEigenDecomposition

  ! Inverse covariance matrix (using pseudo inverse when singular)
  !   output lnDet = 0.5*log(determinant(2*pi*C))
  integer function InverseSymMatrix(C, invC, lnDet, doEigen) result(info)
    real(c_double), intent(in) :: C(:,:)
    real(c_double), intent(out) :: invC(:,:)
    real(c_double), intent(out) :: lnDet
    logical, intent(in), optional :: doEigen
    logical :: eigen
    integer :: ND, i
    external :: dsytrd, dorgtr, dsteqr, dpotrf, dpotri
    ND = size(C,1)
    eigen = .FALSE.
    if(present(doEigen)) eigen = doEigen
    invC = C
    info = 1
    ! Compute Cholesky's decomposition (if eigen decomposition not requested)
    if(.NOT. eigen) call dpotrf('U', ND, invC, ND, info)
    if(info == 0) then
      ! The fastest method for matrix inversion
      ! At this step invC contains the Cholesky decomposition
      ! So det(invC) = prod invC_ii = sqrt(det(C))
      ! ln(sqrt(det(2*pi*C))) = ND/2*ln(2*pi) + sum(ln(invC_ii))
      lnDet = 0.5d0*ND*log(2.0*mlf_PI)
      do i=1,ND
        lnDet = lnDet + log(invC(i,i))
      end do
      ! Compute inverse of C using inversion of the triangular Cholesky's decomposition
      call dpotri('U', ND, invC, ND, info)
      call SymmetrizeMatrixUp(invC)
    else
      ! More stable method for inversing matrix
      block
        ! If the matrix is singular, we use the more stable eigenvalue decomposition
        ! Can occurs when ND is greater or equivalent to dataset size
        ! If lots of problem occurs in this case, you shall reduce the dimensionality (with e.g. PCA)
        real(c_double) :: LD(ND), LE(ND-1), TAU(ND), WORK(ND,ND), valM, penality
        integer :: LWORK
        LWORK = ND*ND
        ! Reinit matrix with C (if Cholesky's decomposition has been iterated)
        if(.NOT. eigen) invC = C
        ! Use same functions as dsyev
        ! Bidiagonal symmetric transformation
        call dsytrd('U', ND, invC, ND, LD, LE, TAU, WORK, LWORK, info)
        if(info /= 0) RETURN
        ! Generate orthogonal matrix using reflectors generated by dsytrd
        call dorgtr('U', ND, invC, ND, TAU, WORK, LWORK, info)
        if(info /= 0) RETURN
        ! Compute eigen values and orthogonal matrix
        call dsteqr('V', ND, LD, LE, invC, ND, WORK, info)
        if(info /= 0) RETURN
        ! invC contains the orthogonal matrix and LD the diagonal of the matrix

        ! Set a penality factor that will be used on null diagonal elements
        ! This factor sets probability of elements outside the hyper plane to a tiny value
        ! (useful in high dimension cases)
        ! We replace all diagonal value below valM = 10e-9*max(LD) by valM
        valM = 10e-9*maxval(LD)
        penality = 1d0/valM ! Inverse of valM
        where (LD <= valM)
          LD = penality
        elsewhere
          LD = 1d0/LD
        end where
        lnDet = 0.5d0*(ND*log(2.0*mlf_PI)+sum(log(LD)))
        ! For a symmetric matrix, the pseudo determinant is the product of non null eigenvalues
        ! We compute the pseudo-inverse matrix with a penality value on element that are nul
        ! invC = O * diag(pseudoInv(D)) * O**T
        do i=1,ND
          work(:,i) = invC(:,i)*LD(i)
        end do
        invC = matmul(work, transpose(invC))
      end block
    endif
  end function InverseSymMatrix

  integer function SqrtSymMatrix(C, C12) result(info)
    integer :: ND, i
    real(c_double), intent(in) :: C(:,:)
    real(c_double), intent(out) :: C12(:,:)
    real(c_double), allocatable :: LD(:)
    ND = size(C,1)
    ALLOCATE(LD(ND))
    info = SymMatrixEigenDecomposition(C, LD, C12)
    if(info /= 0) RETURN
    LD = sqrt(LD)
    forall (i=1:ND) C12(:, i) = C12(:, i)*LD(i)
  end function SqrtSymMatrix
End Module mlf_matrix
