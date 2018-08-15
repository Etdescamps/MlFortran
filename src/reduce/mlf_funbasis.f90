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

Module mlf_funbasis
  Use ieee_arithmetic
  Use iso_c_binding
  Use iso_fortran_env
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_utils
  Use mlf_matrix
  Use mlf_models
  Use mlf_fun_intf
  Use mlf_errors
  IMPLICIT NONE
  PRIVATE

  real(c_double), parameter :: fb_icoeff(4) =  [3d0/8d0, 7d0/6d0, 23d0/24d0, 1d0]
  ! Object handling timer and step evaluations
  Type, Public, extends(mlf_obj_model) :: mlf_algo_funbasis
    class(mlf_basis_fun), pointer :: fun ! Reference function
    real(c_double), pointer :: P(:,:) ! Selected parameter function basis
    real(c_double), pointer :: W(:,:) ! Selected function basis
    real(c_double), pointer :: X(:) ! Selected inputs of X
    real(c_double), pointer :: Vals(:,:) ! Values of the function basis
    real(c_double) :: alpha, x0, xEnd, eA0, eDiff, eAEnd
  Contains
    procedure :: initF => mlf_funbasis_init
    procedure :: getValueBasis => mlf_FunBasisValue
  End Type mlf_algo_funbasis

  Type, Public, extends(mlf_approx_linear) :: mlf_model_funbasis
    class(mlf_algo_funbasis), pointer :: top
  Contains
    procedure :: getProjMult => mlf_FunBasisGetProjection
    procedure :: getProjSingle => mlf_FunBasisGetProjectionSingle
    procedure :: getValueBasis => mlf_model_FunBasisValue
    procedure :: getValueBounds => mlf_model_FunBasisBounds
  End Type mlf_model_funbasis


Contains

  Subroutine mlf_model_FunBasisBounds(this, xMin, xMax)
    class(mlf_model_funbasis), intent(in) :: this
    real(c_double), intent(out), optional :: xMin, xMax
    if(PRESENT(xMin)) xMin = this%top%x0
    if(PRESENT(xMax)) xMax = this%top%xEnd
  End Subroutine mlf_model_FunBasisBounds

  ! Model initilisator
  Subroutine mlf_model_funbasis_init(this, top)
    class(mlf_model), intent(out), allocatable :: this
    class(mlf_algo_funbasis), intent(in), target :: top
    ALLOCATE(mlf_model_funbasis :: this)
    select type(this)
      class is (mlf_model_funbasis)
        this%top => top
    end select
  End Subroutine mlf_model_funbasis_init

  ! C wrapper for init function
  type(c_ptr) Function c_funbasis_init(cfobj, nFPar, alpha, x0, xEnd, cP, nP, sizeBase, nX, nXA, cWP) &
      bind(C, name="mlf_funBasisInit")
    type(c_ptr), value :: cfobj, cP, cWP
    class(mlf_obj), pointer :: fobj
    type(mlf_algo_funbasis), pointer :: x
    class (*), pointer :: obj
    real(c_double), value :: alpha, x0, xEnd
    integer(c_int), value :: nFPar, nP, sizeBase, nX, nXA
    real(c_double), pointer :: WP(:), P(:,:)
    integer :: info
    c_funbasis_init = C_NULL_PTR
    fobj => mlf_getobjfromc(cfobj)
    if(.NOT. associated(fobj)) RETURN
    select type(fobj)
      class is (mlf_basis_fun)
        if(.NOT. C_ASSOCIATED(cP)) RETURN
        call C_F_POINTER(cP, P, [nFPar, nP])
        ALLOCATE(x)
        if(.NOT. C_ASSOCIATED(cWP)) then
          info = x%initF(fobj, alpha, x0, xEnd, P, sizeBase, nX, nXA)
        else
          call C_F_POINTER(cWP, WP, [nP])
          info = x%initF(fobj, alpha, x0, xEnd, P, sizeBase, nX, nXA, WP)
        endif
        obj => x
        c_funbasis_init = c_allocate(obj)
    end select
  End Function c_funbasis_init

  ! Init function for the step object
  integer Function mlf_funbasis_init(this, f, alpha, x0, xEnd, P, sizeBase, nX, nXA, WP, data_handler) result(info)
    class(mlf_algo_funbasis), intent(inout) :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    class(mlf_basis_fun), intent(inout), optional, target :: f
    real(c_double), intent(in), optional :: alpha, x0, xEnd
    real(c_double), intent(in), optional :: WP(:), P(:,:)
    integer, intent(in), optional :: sizeBase, nX, nXA
    real(c_double), allocatable :: C(:,:), LD(:), LB(:,:)
    type(mlf_rsc_numFields) :: numFields
    integer :: i, j, nP, N
    integer(c_int64_t) :: ndP(2), ndW(2), ndV(2), nX0
    numFields = mlf_rsc_numFields(0,0,4)
    info = mlf_arr_init(this, numFields, data_handler)
    if(PRESENT(P) .AND. PRESENT(sizeBase) .AND. PRESENT(nX)) then
      nP = size(P,1); N = size(P,2)
      ndP = int([nP, N], kind=8); ndW = int([N, sizeBase], kind=8)
      if(PRESENT(nXA)) then
        ndV = int([sizeBase, nXA], kind=8)
      else
        ndV = int([sizeBase, nX], kind=8)
      endif
    else if(.NOT. PRESENT(data_handler)) then
      write (error_unit, *) 'mlf_funbasis(init) error: no input parameter nor data_handler'
      info = -1; RETURN
    endif
    info = this%add_rmatrix(numFields, ndP, this%P, C_CHAR_"P", data_handler = data_handler)
    if(info /= 0) RETURN
    ndW(1) = ndP(2)
    info = this%add_rmatrix(numFields, ndW, this%W, C_CHAR_"W", &
      data_handler = data_handler, fixed_dims = [.TRUE., .FALSE.])
    if(info /= 0) RETURN
    ndV(1) = ndW(2)
    info = this%add_rmatrix(numFields, ndV, this%Vals, C_CHAR_"Vals", &
      data_handler = data_handler, fixed_dims = [.TRUE., .FALSE.])
    if(info /= 0) RETURN
    nX0 = ndV(2)
    info = this%add_rarray(numFields, nX0, this%X, C_CHAR_"X", &
      data_handler = data_handler, fixed_dims = [.TRUE.])
    if(info /= 0) RETURN
    if(PRESENT(data_handler)) RETURN
    this%fun => f
    this%P = P
    allocate(C(N,N), LD(N), LB(N,N))
    CALL InitOrDefault(this%xEnd, ieee_value(0d0, ieee_positive_inf), xEnd)
    if(ieee_is_finite(this%xEnd)) then
      CALL InitOrDefault(this%alpha, 0d0, alpha)
    else
      CALL InitOrDefault(this%alpha, 1d0, alpha)
    endif
    CALL InitOrDefault(this%x0, 0d0, x0)
    CALL ComputeFunMatrix(this, C, nX)
    if(PRESENT(WP)) then
      forall(i=1:N,j=1:N) C(i,j) = WP(i)*WP(j)*C(i,j)
    endif
    ! C is the positive definite matrix containing the dot product of each function
    info = SymMatrixEigenDecomposition(C, LD, LB)
    if(info /= 0) RETURN
    ! Choose the (normalised) eigenvector that have the higest eigenvalues
    if(PRESENT(WP)) then
      forall(i=1:sizeBase) this%W(:,i) = LB(:,(N-i+1))/sqrt(LD(N-i+1))*WP(:)
    else
      forall(i=1:sizeBase) this%W(:,i) = LB(:,(N-i+1))/sqrt(LD(N-i+1))
    endif
    CALL ComputeBasisValue(this, int(nX0,4))
    CALL mlf_model_funbasis_init(this%model, this)
  End Function mlf_funbasis_init

  ! Subroutines used for computing a huge number of dot product between function
  ! The goal is to have a restricted basis of function

  ! First subroutine: compute for one particular X all possibles
  ! local dot product
  integer Function ComputeCoeff(fun, x, P, Y, Y0, fact, M) result(info)
    class(mlf_basis_fun), intent(inout), target :: fun
    real(c_double), intent(in) :: P(:,:), x, fact
    real(c_double), intent(out) :: Y(:,:)
    real(c_double), intent(inout) :: Y0(:)
    integer, intent(in) :: M
    info = fun%eval([x], P, Y)
    if(info<0) RETURN
    call ComputeProduct(Y(1,:), Y0, fact, M)    
  End Function ComputeCoeff
  ! First subroutine: compute for one particular X all possibles
  ! local dot products
  Subroutine ComputeProduct(Y, Y0, fact, M)
    real(c_double), intent(in) :: Y(:), fact
    real(c_double), intent(inout) :: Y0(:)
    integer, intent(in) :: M
    integer :: i, j, k
    k = 1
    do i=1,M
      do j=i,M
        Y0(k) = Y0(k) + fact*Y(i)*Y(j)
        k = k+1
      end do
    end do
  End Subroutine ComputeProduct
  ! Second subroutine: compute dot product between the value of a set of function
  ! with a orthonormal function vector of the function basis
  real(c_double) Function ComputeDotProduct(Y, Y0, xEnd) result(F)
    real(c_double), intent(in) :: Y(:), Y0(:), xEnd
    integer :: i, N, l
    F =  0
    N =  size(Y)
    l =  merge(0,1, ieee_is_finite(xEnd))
    do i=1,3
      F = F + fb_icoeff(i+l)*Y(i)*Y0(i)
    end do
    do i=3,N-3
      F = F + Y(i)*Y0(i)
    end do
    do i=N-2,N
      F = F + fb_icoeff(N+1-i)*Y(i)*Y0(i)
    end do
  End Function ComputeDotProduct
  ! Compute the matrix of the dot product between all functions
  ! Main function for determining a function basis
  ! Gives also the error
  Subroutine ComputeFunMatrix(this, C, N)
    class(mlf_algo_funbasis), intent(in) :: this
    real(c_double), intent(out) :: C(:,:)
    integer, intent(in) :: N
    ! constants initialised once
    real(c_double) :: invAlpha, eDiff, difFact, eA0, eAEnd
    real(c_double), allocatable :: Y(:,:), Y0(:)
    integer :: i, j, k, M, P, info
    M = size(C,1)
    P = M*(M+1)/2
    allocate(Y(1,M), Y0(P))
    if(this%alpha == 0) then
      eA0 = this%x0
      eAEnd = this%xEnd
      eDiff = (eAEnd-eA0)/real(N-1, kind=8)
      difFact = eDiff
      Y0 = 0
      ! We made the sum backward for avoiding catastrophic cancellation
      info = ComputeCoeff(this%fun, this%xEnd, this%P, Y, Y0, fb_icoeff(1)*difFact, M)
      info = ComputeCoeff(this%fun, eAEnd-eDiff, this%P, Y, Y0, fb_icoeff(2)*difFact, M)
      info = ComputeCoeff(this%fun, eAEnd-2*eDiff, this%P, Y, Y0, fb_icoeff(3)*difFact, M)
      do i=3,N-4
        info = ComputeCoeff(this%fun, eAEnd-i*eDiff, this%P, Y, Y0, difFact, M)
      end do
      info = ComputeCoeff(this%fun, eA0+2*eDiff, this%P, Y, Y0, fb_icoeff(3)*difFact, M)
      info = ComputeCoeff(this%fun, eA0+eDiff, this%P, Y, Y0, fb_icoeff(2)*difFact, M)
      info = ComputeCoeff(this%fun, this%x0, this%P, Y, Y0, fb_icoeff(1)*difFact, M)
    else
      eA0 = exp(-this%alpha*this%x0)
      invAlpha = 1/this%alpha
      if(ieee_is_finite(this%xEnd)) then
        eAEnd = exp(-this%alpha*this%xEnd)
      else
        eAEnd = 0
      endif
      eDiff = (eA0-eAEnd)/real(N-1, kind=8)
      difFact = eDiff*invAlpha
      Y0 = 0
      ! We made the sum backward for avoiding catastrophic cancellation
      if(eAEnd > 0) then
        info = ComputeCoeff(this%fun, this%xEnd, this%P, Y, Y0, fb_icoeff(1)*difFact, M)
      else
        Y0 = 0
      endif
      info = ComputeCoeff(this%fun, -invAlpha*log(eAEnd+eDiff), this%P, Y, Y0, fb_icoeff(2)*difFact, M)
      info = ComputeCoeff(this%fun, -invAlpha*log(eAEnd+2*eDiff), this%P, Y, Y0, fb_icoeff(3)*difFact, M)
      do i=3,N-4
        info = ComputeCoeff(this%fun, -invAlpha*log(eAEnd+i*eDiff), this%P, Y, Y0, difFact, M)
      end do
      info = ComputeCoeff(this%fun, -invAlpha*log(eA0-2*eDiff), this%P, Y, Y0, fb_icoeff(3)*difFact, M)
      info = ComputeCoeff(this%fun, -invAlpha*log(eA0-eDiff), this%P, Y, Y0, fb_icoeff(2)*difFact, M)
      info = ComputeCoeff(this%fun, this%x0, this%P, Y, Y0, fb_icoeff(1)*difFact, M)
    endif
    k = 1
    do i=1,M
      do j=i,M
        C(i,j) = Y0(k)
        k = k+1
      end do
    end do
    do i=2,M
      forall(j=1:(i-1)) C(i,j) = C(j,i)
    end do
  End Subroutine ComputeFunMatrix

    integer Function mlf_FunBasisGetProjectionSingle(this, Y, W, Aerror) result(info)
    ! Get the projection of the function in the basis of the selected vector
    class(mlf_model_funbasis), intent(in), target :: this
    real(c_double), intent(in) :: Y(:)
    real(c_double), intent(out) :: W(:)
    real(c_double), optional, intent(out) :: Aerror(:)
    real(c_double), allocatable :: F(:,:), P(:)
    integer :: np, nx, i
    real(c_double) :: coeff
    np = size(this%top%W, 2)
    nx = size(this%top%X)
    allocate(F(nx,1), P(nx))
    info = this%top%fun%eval(this%top%X, reshape(Y, [size(Y),1]), F)
    if(info<0) RETURN
    if(this%top%alpha == 0) then
      coeff = this%top%eDiff
    else
      coeff = this%top%eDiff/this%top%alpha
    endif
    do i = 1,np
      P =  this%top%Vals(i,:)
      W(i) = coeff*ComputeDotProduct(F(:,1), P, this%top%xEnd)
      ! It is not essential to remove this part (the function basis is orthonormal)
      ! but the remaining value can be used for estimating the error of the approximation
      F(:,1) = F(:,1)-W(i)*P
    end do
    if(present(Aerror)) then
      F = abs(F)
      Aerror(1) = sum(F)/real(nx)
      Aerror(2) = maxval(F)
    endif
  End Function mlf_FunBasisGetProjectionSingle

  integer Function mlf_FunBasisGetProjection(this, Y, W, Aerror) result(info)
    ! Get the projection of the function in the basis of the selected vector
    class(mlf_model_funbasis), intent(in), target :: this
    real(c_double), intent(in) :: Y(:,:)
    real(c_double), intent(out) :: W(:,:)
    real(c_double), optional, intent(out) :: Aerror(:,:)
    real(c_double), allocatable :: F(:,:), P(:)
    integer :: np, nx, ny, i, j
    real(c_double) :: coeff
    np = size(this%top%W, 2)
    ny = size(Y, 2)
    nx = size(this%top%X)
    allocate(F(nx,ny), P(nx))
    info = this%top%fun%eval(this%top%X, Y, F)
    if(info<0) RETURN
    if(this%top%alpha == 0) then
      coeff = this%top%eDiff
    else
      coeff = this%top%eDiff/this%top%alpha
    endif
    do i = 1,np
      P =  this%top%Vals(i,:)
      do j = 1,ny
        W(i,j) = coeff*ComputeDotProduct(F(:,j), P, this%top%xEnd)
        ! It is not essential to remove this part (the function basis is orthonormal)
        ! but the remaining value can be used for estimating the error of the approximation
        F(:,j) = F(:,j)-W(i,j)*P
      end do
    end do
    if(present(Aerror)) then
      F = abs(F)
      forall (i=1:ny)
        Aerror(i,1) = sum(F(:,i))/real(nx)
        Aerror(i,2) = maxval(F(:,i))
      end forall 
    endif
  End Function mlf_FunBasisGetProjection

  Integer Function mlf_model_FunBasisValue(this, x, Y) result(info)
    class(mlf_model_funbasis), intent(in) :: this
    real(c_double), intent(in) :: x
    real(c_double), intent(out) :: Y(:)
    info = this%top%getValueBasis(x, Y)
  End Function mlf_model_FunBasisValue

  Integer Function mlf_FunBasisValue(this, x, Y) result(info)
    class(mlf_algo_funbasis), intent(in) :: this
    real(c_double), intent(in) :: x
    real(c_double), intent(out) :: Y(:)
    real(c_double) :: t, dX, X0(4)
    integer :: i
    info = 0
    ! The X are dreasingly ordered
    if(this%alpha == 0) then
      i = INT((x-this%x0)/this%eDiff)+1
      if(i < -1 .OR. i > size(this%X)+2) then
        info = 1 ! Interpolation too far to be accurate
        RETURN
      endif
      i = max(2,min(i,size(this%X)-2))
      X0 = (x-this%X(i-1:i+2))/this%eDiff
      ASSOCIATE(Y1 => this%Vals(:,i-1), Y2 => this%Vals(:,i), Y3 => this%Vals(:,i+1), &
          Y4 => this%Vals(:,i+2))
        Y = 1d0/6d0*X0(2)*X0(3)*(X0(1)*Y4-X0(4)*Y1)+0.5d0*X0(1)*X0(4)*(X0(3)*Y2-X0(2)*Y3)
      END ASSOCIATE
    else
      if(x < this%x0) then
        info = 1
        RETURN
      endif
      if(x>=this%X(1)) then
        Y = this%Vals(:,1) &
          * exp(log(this%Vals(:,1)/this%Vals(:,2))/(this%X(1)-this%X(2))*(x-this%X(1)))
        RETURN
      endif
      i = 1+INT((exp(-this%alpha*x)-this%eAEnd)/this%eDiff)
      if(i >= size(this%X, 1)) i = size(this%X, 1)-1
      dX = this%X(i)-this%X(i+1)
      ! We define t as (x-X(i+1))/dX
      !  so when t=0 -> x=X(i+1) and t=1 -> x=X(i)
      t = (x-this%X(i+1))/dX
      ! We get values Y1=F_W(X(i+1)) and Y2=F_W(X(i))
      Y =  (1d0-t)*this%Vals(:,i+1)+t*this%Vals(:,i)
    endif
  End Function mlf_FunBasisValue

  ! Compute the vectors X and V from the structure
  subroutine ComputeBasisValue(this, N)
    class(mlf_algo_funbasis), intent(inout) :: this
    integer, intent(in) :: N
    integer :: M, i, info
    real(c_double) :: eA0, eAEnd, eDiff, invAlpha
    real(c_double), allocatable :: Y(:,:)
    M = size(this%P,2)
    ALLOCATE(Y(1,M))
    if(this%alpha == 0) then
      this%eDiff = (this%xEnd-this%x0)/real(N-1, kind=8)
      forall(i=1:N) this%X(i) = real(i-1, kind=8)*this%eDiff+this%x0 
    else
      eA0 = exp(-this%alpha*this%x0)
      invAlpha = 1d0/this%alpha
      if(ieee_is_finite(this%xEnd)) then
        eAEnd = exp(-this%alpha*this%xEnd)
        eDiff = (eA0-eAEnd)/real(N-1, kind=8)
        this%X(1) = this%xEnd
      else
        eDiff = eA0/real(N, kind=8)
        eAEnd = eDiff
        this%X(1) = -invAlpha*log(eDiff)
      endif
      this%eA0 = eA0
      this%eAEnd = eAEnd
      this%eDiff = eDiff
      this%X(N) = this%x0
      do i = 2,N-1
        this%X(i) = -invAlpha*log(eAEnd+(i-1)*eDiff)
      end do
    endif
    do i = 1,N
      info = this%fun%eval(this%X(i:i), this%P, Y)
      this%Vals(:,i) = matmul(Y(1,:),this%W)
    end do
  end subroutine ComputeBasisValue
End Module mlf_funbasis
