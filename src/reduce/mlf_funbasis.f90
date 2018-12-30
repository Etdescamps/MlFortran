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

  ! Coefficients from Numerical Recipies in C (Press, Teukolsky) Chapter 4.1
  Real(c_double), Parameter :: fb_icoeff(4) =  [3d0/8d0, 7d0/6d0, 23d0/24d0, 1d0]

  ! Object handling timer and step evaluations
  Type, Public, Extends(mlf_obj_model) :: mlf_algo_funbasis
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

  Type, Public, extends(mlf_fast_approx_linear) :: mlf_model_funbasis
    class(mlf_algo_funbasis), pointer :: top
  Contains
    procedure :: getProjMult => mlf_FunBasisGetProjection
    procedure :: getProjSingle => mlf_FunBasisGetProjectionSingle
    procedure :: getValueBasis => mlf_model_FunBasisValue
    procedure :: getValueBounds => mlf_model_FunBasisBounds
    procedure :: getFastProj => mlf_FastProjection
  End Type mlf_model_funbasis
Contains
  Subroutine mlf_model_FunBasisBounds(this, xMin, xMax)
    class(mlf_model_funbasis), intent(in) :: this
    real(c_double), intent(out), optional :: xMin, xMax
    If(PRESENT(xMin)) xMin = this%top%x0
    If(PRESENT(xMax)) xMax = this%top%xEnd
  End Subroutine mlf_model_FunBasisBounds

  ! Model initilisator
  Subroutine mlf_model_funbasis_init(this, top)
    class(mlf_model), intent(out), allocatable :: this
    class(mlf_algo_funbasis), intent(in), target :: top
    ALLOCATE(mlf_model_funbasis :: this)
    Select Type(this)
    Class is (mlf_model_funbasis)
      this%top => top
    End Select
  End Subroutine mlf_model_funbasis_init

  ! C wrapper for init function
  Type(c_ptr) Function c_funbasis_init(cfobj, nFPar, alpha, x0, xEnd, cP, nP, sizeBase, nX, nXA, cWP) &
      Bind(C, name="mlf_funBasisInit")
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
    If(.NOT. ASSOCIATED(fobj)) RETURN
    Select Type(fobj)
    Class is (mlf_basis_fun)
      If(.NOT. C_ASSOCIATED(cP)) RETURN
      CALL C_F_POINTER(cP, P, [nFPar, nP])
      ALLOCATE(x)
      If(.NOT. C_ASSOCIATED(cWP)) Then
        info = x%initF(fobj, alpha, x0, xEnd, P, sizeBase, nX = nX, nXA = nXA)
      Else
        CALL C_F_POINTER(cWP, WP, [nP])
        info = x%initF(fobj, alpha, x0, xEnd, P, sizeBase, nX = nX, nXA = nXA, WP = WP)
      Endif
      obj => x
      c_funbasis_init = c_allocate(obj)
    End Select
  End Function c_funbasis_init

  Integer Function FunBasisInit(this, sizeBase, nX, nXC, WP) Result(info)
    class(mlf_algo_funbasis), intent(inout) :: this
    integer, intent(in) :: sizeBase, nXC, nX
    real(c_double), intent(in), optional :: WP(:)
    real(c_double), allocatable :: C(:,:), LD(:), LB(:,:)
    integer :: i, j, N
    N = SIZE(this%P, 2)
    ALLOCATE(C(N,N), LD(sizeBase), LB(N,sizeBase))
    ! Compute the matrix of the dot product between each function <f_i,f_j>
    ASSOCIATE(f => this%fun)
      Select Type(f)
      Class is (mlf_basis_fun_inv)
        info = ComputeFunMatrixInv(f, this%x0, this%xEnd, this%P, C, nXC)
      Class Default
        info = ComputeFunMatrix(f, this%x0, this%xEnd, this%alpha, this%P, C, nXC)
      End Select
    END ASSOCIATE
    If(info < 0) RETURN
    !
    If(PRESENT(WP)) Then
      FORALL(i=1:N,j=1:N) C(i,j) = WP(i)*WP(j)*C(i,j)
    Endif
    ! C is the positive definite matrix containing the dot product of each function
    info = SymMatrixSelectEigenValues(C, LD, LB)
    If(info /= 0) GOTO 10
    ! Choose the (normalised) eigenvector that have the higest eigenvalues
    If(PRESENT(WP)) Then
      FORALL(i=1:sizeBase) this%W(:,i) = LB(:,sizeBase-i+1)/sqrt(LD(sizeBase-i+1))*WP
    Else
      FORALL(i=1:sizeBase) this%W(:,i) = LB(:,sizeBase-i+1)/sqrt(LD(sizeBase-i+1))
    Endif
    CALL ComputeBasisValue(this, nX)
    RETURN
 10 WRITE (error_unit, *) "mlf_funbasis: error with Eigen decomposition"
    info = -1
  End Function FunBasisInit

  ! Init function for the step object
  Integer Function mlf_funbasis_init(this, f, alpha, x0, xEnd, P, sizeBase, &
      nX, nXA, nXC, WP, data_handler) Result(info)
    class(mlf_algo_funbasis), intent(inout) :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    class(mlf_basis_fun), intent(inout), target :: f
    real(c_double), intent(in), optional :: alpha, x0, xEnd
    real(c_double), intent(in), optional :: WP(:), P(:,:)
    integer, intent(in), optional :: sizeBase, nX, nXC, nXA
    type(mlf_rsc_numFields) :: numFields
    integer :: nP, N, nXC0
    integer(c_int64_t) :: ndP(2), ndW(2), ndV(2), nX0
    N = -1
    numFields = mlf_rsc_numFields(0,0,4)
    info = mlf_arr_init(this, numFields, data_handler)
    If(PRESENT(P) .AND. PRESENT(sizeBase) .AND. PRESENT(nX)) Then
      nP = SIZE(P,1); N = SIZE(P,2)
      ndP = INT([nP, N], kind=8); ndW = INT([N, sizeBase], kind=8)
      If(PRESENT(nXA)) then
        ndV = INT([sizeBase, nXA], kind=8)
      Else
        ndV = INT([sizeBase, nX], kind=8)
      Endif
    Else If(.NOT. PRESENT(data_handler)) Then
      WRITE (error_unit, *) 'mlf_funbasis(init) error: no input parameter nor data_handler'
      info = -1; RETURN
    Endif
    info = this%add_rmatrix(numFields, ndP, this%P, C_CHAR_"P", data_handler = data_handler)
    If(info /= 0) RETURN
    If(N<0) N = INT(ndP(2), kind=4)
    ndW(1) = ndP(2)
    info = this%add_rmatrix(numFields, ndW, this%W, C_CHAR_"W", &
      data_handler = data_handler, fixed_dims = [.TRUE., .FALSE.])
    If(info /= 0) RETURN
    ndV(1) = ndW(2)
    info = this%add_rmatrix(numFields, ndV, this%Vals, C_CHAR_"Vals", &
      data_handler = data_handler, fixed_dims = [.TRUE., .FALSE.])
    If(info /= 0) RETURN
    nX0 = ndV(2)
    info = this%add_rarray(numFields, nX0, this%X, C_CHAR_"X", &
      data_handler = data_handler, fixed_dims = [.TRUE.])
    If(info /= 0) RETURN
    this%fun => f
    If(PRESENT(P)) this%P = P
    If(.NOT. PRESENT(data_handler)) Then
      nXC0 = nX
      If(PRESENT(nXC)) nXC0 = nXC
      If(.NOT. PRESENT(P)) GOTO 10
      CALL InitOrDefault(this%xEnd, ieee_value(0d0, ieee_positive_inf), xEnd)
      If(ieee_is_finite(this%xEnd)) Then
        CALL InitOrDefault(this%alpha, 0d0, alpha)
      Else
        CALL InitOrDefault(this%alpha, 1d0, alpha)
      Endif
      CALL InitOrDefault(this%x0, 0d0, x0)
      info = FunBasisInit(this, sizeBase, INT(nX0, KIND=4), nXC0, WP)
    Endif
    CALL mlf_model_funbasis_init(this%model, this)
    RETURN
 10 WRITE (error_unit, *) 'mlf_funbasis error: missing required parameter'
    info = -1
  End Function mlf_funbasis_init

  ! Subroutines used for computing a huge number of dot product between function
  ! The goal is to have a restricted basis of function

  ! First subroutine: compute for one particular X all possibles
  ! local dot product
  Integer Function ComputeCoeff(fun, x, P, Y, Y0, fact, M) Result(info)
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
  Real(c_double) Function ComputeDotProduct(Y, Y0, xEnd) Result(F)
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
  Integer Function ComputeFunMatrix(fun, x0, xEnd, alpha, Param, C, N) Result(info)
    class(mlf_basis_fun), intent(inout), target :: fun
    real(c_double), intent(in) :: x0, xEnd, alpha, Param(:,:)
    real(c_double), intent(out) :: C(:,:)
    integer, intent(in) :: N
    ! constants initialised once
    real(c_double) :: invAlpha, eDiff, difFact, eA0, eAEnd
    real(c_double), allocatable :: Y(:,:), Y0(:)
    integer :: i, j, k, M, P
    info = 0
    M = SIZE(C,1)
    P = M*(M+1)/2
    ALLOCATE(Y(1,M), Y0(P))
    If(alpha == 0) Then
      eA0 = x0
      eAEnd = xEnd
      eDiff = (eAEnd-eA0)/REAL(N-1, kind=8)
      difFact = eDiff
      Y0 = 0
      ! We made the sum backward for avoiding catastrophic cancellation
      info = ComputeCoeff(fun, xEnd, Param, Y, Y0, fb_icoeff(1)*difFact, M)
      info = ComputeCoeff(fun, eAEnd-eDiff, Param, Y, Y0, fb_icoeff(2)*difFact, M)
      info = ComputeCoeff(fun, eAEnd-2*eDiff, Param, Y, Y0, fb_icoeff(3)*difFact, M)
      Do i=3,N-4
        info = ComputeCoeff(fun, eAEnd-i*eDiff, Param, Y, Y0, difFact, M)
      End Do
      info = ComputeCoeff(fun, eA0+2*eDiff, Param, Y, Y0, fb_icoeff(3)*difFact, M)
      info = ComputeCoeff(fun, eA0+eDiff, Param, Y, Y0, fb_icoeff(2)*difFact, M)
      info = ComputeCoeff(fun, x0, Param, Y, Y0, fb_icoeff(1)*difFact, M)
    Else
      eA0 = EXP(-alpha*x0)
      invAlpha = 1/alpha
      If(ieee_is_finite(xEnd)) Then
        eAEnd = EXP(-alpha*xEnd)
      Else
        eAEnd = 0
      Endif
      eDiff = (eA0-eAEnd)/REAL(N-1, kind=8)
      difFact = eDiff*invAlpha
      Y0 = 0
      ! We made the sum backward for avoiding catastrophic cancellation
      If(eAEnd > 0) Then
        info = ComputeCoeff(fun, xEnd, Param, Y, Y0, fb_icoeff(1)*difFact, M)
      Else
        Y0 = 0
      Endif
      info = ComputeCoeff(fun, -invAlpha*LOG(eAEnd+eDiff), Param, Y, Y0, fb_icoeff(2)*difFact, M)
      info = ComputeCoeff(fun, -invAlpha*LOG(eAEnd+2*eDiff), Param, Y, Y0, fb_icoeff(3)*difFact, M)
      Do i=3,N-4
        info = ComputeCoeff(fun, -invAlpha*LOG(eAEnd+i*eDiff), Param, Y, Y0, difFact, M)
      End Do
      info = ComputeCoeff(fun, -invAlpha*LOG(eA0-2*eDiff), Param, Y, Y0, fb_icoeff(3)*difFact, M)
      info = ComputeCoeff(fun, -invAlpha*LOG(eA0-eDiff), Param, Y, Y0, fb_icoeff(2)*difFact, M)
      info = ComputeCoeff(fun, x0, Param, Y, Y0, fb_icoeff(1)*difFact, M)
    Endif
    k = 1
    Do i=1,M
      Do j=i,M
        C(i,j) = Y0(k)
        k = k+1
      End Do
    End Do
    Do i=2,M
      forall(j=1:(i-1)) C(i,j) = C(j,i)
    End Do
  End Function ComputeFunMatrix

  ! Compute the matrix of the dot product between all functions
  ! Main function for determining a function basis
  ! Gives also the error
  Integer Function ComputeFunMatrixInv(fun, x0, xEnd, Param, C, N) Result(info)
    class(mlf_basis_fun_inv), intent(inout) :: fun
    real(c_double), intent(in) :: x0, xEnd, Param(:,:)
    real(c_double), intent(out) :: C(:,:)
    integer, intent(in) :: N
    ! constants initialised once
    real(c_double) :: S, Fact
    real(c_double), allocatable :: Y(:,:),  II(:), Z(:), Y1(:), P(:,:)
    real(c_double), allocatable :: Y0(:,:), X(:), XV(:)
    integer(c_int), allocatable :: idx(:)
    integer :: i, j, k, l, M, infoI
    M = SIZE(C,1)
    ALLOCATE(Y(2,M), II(M), Z(M), Y1(N), P(SIZE(Param, 1), N), idx(M))
    Y1 = [(REAL(i, KIND=8)/REAL(N+1, KIND=8), i=1,N)]
    info = fun%integral(x0, xEnd, Param, II)
    If(info < 0) RETURN
    info = fun%eval([x0, xEnd], Param, Y)
    If(info < 0) RETURN
    Where(II > 0d0)
      Z = ABS(II)/MAXVAL(ABS(Y), DIM=1)
    ElseWhere
      Z = 0d0
    End Where
    CALL QSortIdx(Z, idx)
    P = Param(:,idx)
    Y = Y(:,idx)
    !$OMP PARALLEL default(shared) PRIVATE(Y0, X, XV, i, j, k, l, infoI, S, Fact)
    ALLOCATE(Y0(N,1), X(N), XV(N))
    !$OMP Do schedule(dynamic, 4)
    Do i = 1, M
      k = idx(i)
      If(II(k) /= 0) Then
        XV = Y1*II(k)
        infoI = fun%inv_integ(x0, XV, P(:,i), X)
        Fact = II(k)/REAL(N+1, KIND=8)
        Do j = i, M
          l = idx(j)
          infoI = fun%eval(X, P(:,j:j), Y0)
          S = (fb_icoeff(1)*SUM(Y(:,j)) + fb_icoeff(2)*(Y0(1,1)+Y0(N,1)) &
            +  fb_icoeff(3)*(Y0(2,1)+Y0(N-1,1)) + SUM(Y0(3:N-2,1)))*Fact
          C(k,l) = S
          C(l,k) = S
        End Do
      Else
        C(k,:) = 0d0
        C(:,k) = 0d0
      Endif
    End Do
    !$OMP END Do
    DEALLOCATE(Y0, X, XV)
    !$OMP END PARALLEL
    info = 0
  End Function ComputeFunMatrixInv

  Integer Function mlf_FunBasisGetProjectionSingle(this, Y, W, Aerror) Result(info)
    ! Get the projection of the function in the basis of the selected vector
    class(mlf_model_funbasis), intent(in), target :: this
    real(c_double), intent(in) :: Y(:)
    real(c_double), intent(out) :: W(:)
    real(c_double), optional, intent(out) :: Aerror(:)
    real(c_double), allocatable :: F(:,:), P(:)
    integer :: np, nx, i
    real(c_double) :: coeff
    np = SIZE(this%top%W, 2)
    nx = SIZE(this%top%X)
    ALLOCATE(F(nx,1), P(nx))
    info = this%top%fun%eval(this%top%X, RESHAPE(Y, [SIZE(Y),1]), F)
    If(info<0) RETURN
    If(this%top%alpha == 0) Then
      coeff = this%top%eDiff
    Else
      coeff = this%top%eDiff/this%top%alpha
    Endif
    Do i = 1,np
      P =  this%top%Vals(i,:)
      W(i) = coeff*ComputeDotProduct(F(:,1), P, this%top%xEnd)
      ! It is not essential to remove this part (the function basis is orthonormal)
      ! but the remaining value can be used for estimating the error of the approximation
      F(:,1) = F(:,1)-W(i)*P
    End Do
    If(PRESENT(Aerror)) Then
      F = ABS(F)
      Aerror(1) = SUM(F)/REAL(nx)
      Aerror(2) = MAXVAL(F)
    Endif
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
    np = SIZE(this%top%W, 2)
    ny = SIZE(Y, 2)
    nx = SIZE(this%top%X)
    ALLOCATE(F(nx,ny), P(nx))
    info = this%top%fun%eval(this%top%X, Y, F)
    If(info<0) RETURN
    If(this%top%alpha == 0) Then
      coeff = this%top%eDiff
    Else
      coeff = this%top%eDiff/this%top%alpha
    Endif
    Do i = 1,np
      P =  this%top%Vals(i,:)
      Do j = 1,ny
        W(i,j) = coeff*ComputeDotProduct(F(:,j), P, this%top%xEnd)
        ! It is not essential to remove this part (the function basis is orthonormal)
        ! but the remaining value can be used for estimating the error of the approximation
        F(:,j) = F(:,j)-W(i,j)*P
      End Do
    End Do
    If(PRESENT(Aerror)) Then
      F = ABS(F)
      FORALL (i=1:ny)
        Aerror(i,1) = SUM(F(:,i))/REAL(nx)
        Aerror(i,2) = MAXVAL(F(:,i))
      END FORALL 
    Endif
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
    integer :: i, sX
    info = -1
    sX = SIZE(this%X)
    ! The X are dreasingly ordered
    If(this%alpha == 0) Then
      i = CEILING((x-this%x0)/this%eDiff)
      If(i < -1) RETURN
      If(i > sX+2) Then
        info = 1 ! Indicate that the function evaluation is inacurate
        Y = this%Vals(:,sX)*exp(this%xEnd-x)
      Endif
      i = MAX(2,MIN(i,sX-2))
      X0 = (x-this%X(i-1:i+2))/this%eDiff
      ASSOCIATE(Y1 => this%Vals(:,i-1), Y2 => this%Vals(:,i), Y3 => this%Vals(:,i+1), &
          Y4 => this%Vals(:,i+2))
        Y = 1d0/6d0*X0(2)*X0(3)*(X0(1)*Y4-X0(4)*Y1)+0.5d0*X0(1)*X0(4)*(X0(3)*Y2-X0(2)*Y3)
      END ASSOCIATE
    Else
      If(x < this%x0) RETURN
      If(x>=this%X(1)) Then
        Y = this%Vals(:,1) &
          * EXP(LOG(this%Vals(:,1)/this%Vals(:,2))/(this%X(1)-this%X(2))*(x-this%X(1)))
        info = 0
        RETURN
      Endif
      i = 1+INT((EXP(-this%alpha*x)-this%eAEnd)/this%eDiff)
      If(i >= sX) i = sX-1
      dX = this%X(i)-this%X(i+1)
      ! We define t as (x-X(i+1))/dX
      !   so when t=0 -> x=X(i+1) and t=1 -> x=X(i)
      t = (x-this%X(i+1))/dX
      ! We get values Y1=F_W(X(i+1)) and Y2=F_W(X(i))
      Y =  (1d0-t)*this%Vals(:,i+1)+t*this%Vals(:,i)
    Endif
    info = 0
  End Function mlf_FunBasisValue

  ! Compute the vectors X and V from the structure
  Subroutine ComputeBasisValue(this, N)
    class(mlf_algo_funbasis), intent(inout) :: this
    integer, intent(in) :: N
    integer :: M, i, info
    real(c_double) :: eA0, eAEnd, eDiff, invAlpha, coeff, Z
    real(c_double), allocatable :: Y(:,:)
    M = SIZE(this%P,2)
    ALLOCATE(Y(1,M))
    If(this%alpha == 0) Then
      this%eDiff = (this%xEnd-this%x0)/REAL(N-1, KIND=8)
      FORALL(i=1:N) this%X(i) = REAL(i-1, KIND=8)*this%eDiff+this%x0
      coeff = this%eDiff
    Else
      eA0 = EXP(-this%alpha*this%x0)
      invAlpha = 1d0/this%alpha
      If(ieee_is_finite(this%xEnd)) Then
        eAEnd = EXP(-this%alpha*this%xEnd)
        eDiff = (eA0-eAEnd)/REAL(N-1, KIND=8)
        this%X(1) = this%xEnd
      Else
        eDiff = eA0/REAL(N, KIND=8)
        eAEnd = eDiff
        this%X(1) = -invAlpha*LOG(eDiff)
      Endif
      this%eA0 = eA0
      this%eAEnd = eAEnd
      this%eDiff = eDiff
      this%X(N) = this%x0
      Do i = 2,N-1
        this%X(i) = -invAlpha*LOG(eAEnd+(i-1)*eDiff)
      End Do
      coeff = this%eDiff*this%alpha
    Endif
    Do i = 1,N
      info = this%fun%eval(this%X(i:i), this%P, Y)
      this%Vals(:,i) = MATMUL(Y(1,:),this%W)
    End Do
    ! Fix the values Vals to ensure that each vector has norm 1
    ! The value of the norm is close to 1 but sufficiently different to add an error
    Do i = 1,SIZE(this%Vals,1)
      Z = coeff*ComputeDotProduct(this%Vals(i,:), this%Vals(i,:), this%xEnd)
      this%Vals(i,:) = this%Vals(i,:)/SQRT(Z)
    End Do
  End Subroutine ComputeBasisValue

  Integer Function mlf_FastProjection(this, P, nX, W) Result(info)
    class(mlf_model_funbasis), intent(in) :: this
    real(c_double), intent(in) :: P(:,:)
    integer, intent(in) :: nX
    real(c_double), intent(out) :: W(:,:)
    real(c_double), allocatable :: X(:), IX(:), IY(:)
    integer :: i, nY
    nY = SIZE(P,2)
    ASSOCIATE(fun => this%top%fun, x0 => this%top%x0, xEnd => this%top%xEnd)
      Select Type(fun)
      Class is (mlf_basis_fun_inv)
        ALLOCATE(X(nX), IY(nY), IX(nX))
        info = fun%integral(x0, xEnd, P, IY)
        If(info < 0) RETURN
        X = [(REAL(i, KIND=8)/REAL(nX+1, KIND=8), i=1,nX)]
        Do i = 1, nY
          If(IX(i) /= 0) Then
            info = fun%inv_integ(x0, X*IY(i), P(:,i), IX)
            If(info < 0) RETURN
            CALL ComputeDotProductAbscissa(this%top%X, this%top%eDiff, this%top%Vals, IX, IY(i), W(:,i))
          Else
            W(:,i) = 0d0
          Endif
        End Do
      Class Default
        info = -1
      End Select
    END ASSOCIATE
  End Function mlf_FastProjection

  ! Compute dot product with list of abscissa (case with mlf_basis_fun_inv)
  Subroutine ComputeDotProductAbscissa(X, dX, Vals, Y, IY, F)
    real(c_double), intent(in) :: X(:), Vals(:,:), dX, IY
    real(c_double), intent(in) :: Y(:) ! Shall not include x0 and xEnd (already included in computation)
    real(c_double), intent(out) :: F(:) ! Dimension: size of base
    real(c_double) :: uY, X0(4), Z(SIZE(F))
    integer :: N, nY, nP, i, j
    N = SIZE(X); nP = SIZE(Vals, 2); nY = SIZE(Y)
    F(:) = fb_icoeff(1)*(Vals(:,1)+Vals(:,nP))
    j = 2; uY = X(2)
    Do i = 1, nY
      If(Y(i) > uY) Then
        j = CEILING((Y(i)-X(1))/dX)
        j = MAX(2, MIN(j, N-2))
        uY = X(j)
      Endif
      X0 = (Y(i)-X(j-1:j+2))/dX
      ASSOCIATE(Y1 => Vals(:,j-1), Y2 => Vals(:,j), Y3 => Vals(:,j+1), Y4 => Vals(:,j+2))
        Z = 1d0/6d0*X0(2)*X0(3)*(X0(1)*Y4-X0(4)*Y1)+0.5d0*X0(1)*X0(4)*(X0(3)*Y2-X0(2)*Y3)
      END ASSOCIATE
      If(i == 1 .OR. i == nY-1) Then
        F = F + fb_icoeff(2)*Z
      Else If(i == 2 .OR. i == nY-2) Then
        F = F + fb_icoeff(3)*Z
      Else
        F = F + Z
      Endif
    End Do
    F = F*IY/REAL(nY+1, KIND = 8)
  End Subroutine ComputeDotProductAbscissa
End Module mlf_funbasis

