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

Module mlf_cmaes
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_step_algo
  Use mlf_optim
  Use mlf_fun_intf
  Use mlf_matrix
  Use mlf_rand
  Use mlf_utils
  Use mlf_errors
  IMPLICIT NONE
  PRIVATE
  
  Public :: mlf_cmaes_objcreate

  ! Common class for CMA-ES and MA-ES
  Type, Public, abstract, extends(mlf_optim_obj) :: mlf_matrix_es_obj
    real(c_double), pointer :: M(:,:) ! M square root matrix of covariance
    real(c_double), pointer :: X0(:) ! Center of the sampling: X(:,i) = D(:,i)+X0(:)
    real(c_double), pointer :: ps(:) ! Isotropically path (that comes from the best mu element of Z)
    real(c_double), pointer :: W(:) ! weight used by the methods on best mu element
    real(c_double), pointer :: alphaCov
    real(c_double), allocatable :: Z(:,:) ! Vectors sampled using Norm(0,Id_N)
    real(c_double), allocatable :: D(:,:) ! Previous vectors multiplied by M (D = MZ)
    real(c_double), allocatable :: Dw(:), DM(:,:)
    real(c_double) :: chiN, mueff, cs, cw, c1, dsigma
  Contains
    procedure :: genX => mlf_matrix_es_obj_genX
    procedure :: stopCond => mlf_matrix_es_obj_stopCond
    procedure :: updateY => mlf_matrix_es_obj_updateY
    procedure :: reinit_X => mlf_matrix_es_obj_reinit_X
    procedure :: reinit => mlf_matrix_es_obj_reinit
    procedure :: updateW => mlf_matrix_es_obj_updateW
    procedure(mlf_matrix_updateF) , deferred :: updateMatrix
  End Type mlf_matrix_es_obj

  Type, Public, extends(mlf_matrix_es_obj) :: mlf_maes_obj
  Contains
    procedure :: init => mlf_maes_init
    procedure :: updateMatrix => mlf_maes_updateMatrix
  End Type mlf_maes_obj


  Type, Public, extends(mlf_matrix_es_obj) :: mlf_cmaes_obj
    real(c_double), pointer :: C(:,:) ! Covariance matrix
    real(c_double), pointer :: pc(:) ! Anisotropic evolution path (that comes from D)
    integer(c_int64_t), pointer :: lastCov, covEvery
    real(c_double) :: cp
  Contains
    procedure :: init => mlf_cmaes_init
    procedure :: reinit_X => mlf_cmaes_reinit_X
    procedure :: updateMatrix => mlf_cmaes_updateMatrix
  End Type mlf_cmaes_obj

  Abstract Interface
    Integer Function mlf_matrix_updateF(this)
      import :: mlf_matrix_es_obj
      class(mlf_matrix_es_obj), intent(inout), target :: this
    End Function mlf_matrix_updateF
  End Interface
Contains
  Integer Function mlf_matrix_es_obj_init(this, numFields, data_handler, params) Result(info)
    class(mlf_matrix_es_obj), intent(inout), target :: this
    class(mlf_step_numFields), intent(inout) :: numFields
    class(mlf_data_handler), intent(inout), optional :: data_handler
    class(mlf_optim_param), intent(inout) :: params
    integer(c_int64_t) :: nD, CND(2), mu
    nD = params%fun%nD
    CND = nD
    If(params%lambda <= 0) Then
      params%lambda = 4 + FLOOR(3.0*LOG(REAL(ND)))
    Endif
    CALL numFields%addFields(nRPar = 1, nRsc = 4)
    info = mlf_optim_init(this, numFields, data_handler, params)
    info = this%add_rmatrix(numFields, CND, this%M, C_CHAR_"M", data_handler = data_handler, fixed_dims = [.TRUE., .TRUE.])
    info = this%add_rarray(numFields, nD, this%X0, C_CHAR_"X0", data_handler = data_handler, fixed_dims = [.TRUE.])
    info = this%add_rarray(numFields, nD, this%ps, C_CHAR_"ps", data_handler = data_handler, fixed_dims = [.TRUE.])
    mu = this%mu
    info = this%add_rarray(numFields, mu, this%W, C_CHAR_"W", data_handler = data_handler)
    this%mu = INT(mu, KIND=c_int)
    CALL this%addRPar(numFields, this%alphaCov, "alphaCov")
    ALLOCATE(this%Z(ND, this%lambda), this%D(ND, this%lambda), this%Dw(ND), this%DM(ND,ND))
  End Function mlf_matrix_es_obj_init

  Subroutine mlf_matrix_es_obj_updateW(this)
    class(mlf_matrix_es_obj), intent(inout), target :: this
    integer :: ND
    ND = SIZE(this%Z, 1)
    ASSOCIATE(W => this%W, mueff => this%mueff, alphacov => this%alphacov, &
        chiN => this%chiN, c1 => this%c1, cs => this%cs, cw => this%cw, &
        dsigma => this%dsigma, nR =>REAL(nD, KIND=8))
      W = 1d0/SUM(W) * W
      mueff = 1d0/(SUM(W * W))
      c1 = alphacov/((nR+1.3d0)**2 + mueff)
      cs = (mueff + 2d0)/(mueff + nR + 5d0)
      cw = MIN(1d0 - c1, alphacov*(mueff+1d0/mueff-2d0)/((nR+2d0)**2 + alphacov*mueff*0.5d0))
      dsigma = 1d0 + cs + 2d0*MAX(0d0, SQRT((mueff-1d0)/(nR+1d0))-1d0)
      chiN = SQRT(nR)*(1d0-1d0/(4d0*nR)+1d0/(21d0*nR*nR))
    END ASSOCIATE
  End Subroutine mlf_matrix_es_obj_updateW 

  Integer Function mlf_matrix_es_obj_reinit(this) Result(info)
    class(mlf_matrix_es_obj), intent(inout), target :: this
    info = this%reinit_X()
  End Function mlf_matrix_es_obj_reinit
 
  Integer Function mlf_matrix_es_obj_reinit_X(this, X0, sigma0, funW, alphacovIn) Result(info)
    class(mlf_matrix_es_obj), intent(inout), target :: this
    class(mlf_weight_fun), optional, intent(in) :: funW
    real(c_double), intent(in), optional :: X0(:), sigma0, alphacovIn
    integer :: lambda, mu, i
    If(PRESENT(sigma0)) this%sigma = sigma0
    info = mlf_optim_reinit(this)
    lambda = SIZE(this%Z, 2); mu = SIZE(this%W)
    CALL IdentityMatrix(this%M)
    this%ps = 0
    If(PRESENT(X0)) Then
      this%X0 = X0
    Else
      CALL this%randomStartPoint(this%X0)
    Endif
    this%alphacov = 2d0
    If(PRESENT(alphacovIn)) this%alphacov=alphacovIn
    If(PRESENT(funW)) Then
      CALL funW%eval(this%W, lambda)
    Else
      this%W = LOG(0.5d0*(lambda+1d0))-[(LOG(REAL(i,KIND=8)), i=1,mu)]
    Endif
    CALL this%updateW()
  End Function mlf_matrix_es_obj_reinit_X

  Integer Function mlf_matrix_es_obj_genX(this, ids) Result(info)
    class(mlf_matrix_es_obj), intent(inout), target :: this
    integer(c_int), intent(in), optional :: ids(:) ! Regenerates selected ids
    integer :: i, j
    info = 0
    If(PRESENT(ids)) Then
      Do i=1,SIZE(ids)
        j = ids(i)
        CALL RandN(this%Z(:,j))
        this%D(:,j) = MATMUL(this%M, this%Z(:,j))
        this%X(:,j) = this%sigma*this%D(:,j)+this%X0
      End Do
    Else
      CALL RandN(this%Z)
      Do i=1,SIZE(this%Z,2)
        this%D(:,i) = MATMUL(this%M, this%Z(:,i))
        this%X(:,i) = this%sigma*this%D(:,i)+this%X0
      End Do
    Endif
  End Function mlf_matrix_es_obj_genX
  
  Integer Function mlf_matrix_es_obj_stopCond(this) Result(info)
    class(mlf_matrix_es_obj), intent(inout), target :: this
    info = mlf_optim_stopCond(this)
  End Function mlf_matrix_es_obj_stopCond

  Integer Function mlf_maes_init(this, funW, alphacovIn, data_handler, params) Result(info)
    class(mlf_maes_obj), intent(inout), target :: this
    class(mlf_optim_param), intent(inout) :: params
    class(mlf_weight_fun), optional, intent(in) :: funW
    real(c_double), intent(in), optional :: alphacovIn
    class(mlf_data_handler), intent(inout), optional :: data_handler
    type(mlf_step_numFields) :: numFields
    CALL numFields%initFields()
    info = mlf_matrix_es_obj_init(this, numFields, data_handler, params)
    If(CheckF(info, 'mlf_maes_init: Error init matrix_es')) RETURN
    If(PRESENT(data_handler)) Then
      CALL this%updateW()
    Else
      info = this%reinit_X(params%X0, funW = funW, alphacovIn = alphacovIn)
    Endif
  End Function mlf_maes_init

  Integer Function mlf_cmaes_init(this, funW, alphacovIn, covEvery, data_handler, params) result(info)
    class(mlf_cmaes_obj), intent(inout), target :: this
    class(mlf_optim_param), intent(inout) :: params
    integer(c_int64_t) :: ND, CND(2)
    class(mlf_weight_fun), optional, intent(in) :: funW
    real(c_double), intent(in), optional :: alphacovIn
    integer, intent(in), optional :: covEvery
    class(mlf_data_handler), intent(inout), optional :: data_handler
    integer :: lambda
    type(mlf_step_numFields) :: numFields
    CALL numFields%initFields(nIPar = 1, nRsc = 2, nIVar =1)
    info = mlf_matrix_es_obj_init(this, numFields, data_handler, params)
    If(CheckF(info, 'mlf_cmaes_init: Error init matrix_es')) RETURN
    lambda = SIZE(this%Z, 2); ND = SIZE(this%Z, 1)
    CND = ND
    info = this%add_rmatrix(numFields, CND, this%C, C_CHAR_"C", &
      data_handler = data_handler, fixed_dims = [.TRUE., .TRUE.])
    info = this%add_rarray(numFields, ND, this%pc, C_CHAR_"pc", &
      data_handler = data_handler, fixed_dims = [.TRUE.])
    CALL this%addIVar(numFields, this%lastCov, "lastCov")
    CALL this%addIPar(numFields, this%covEvery, "covEvery")
    If(PRESENT(data_handler)) Then
      CALL this%updateW()
      this%cp = (this%mueff/ND+4.0D0)/(2.0D0*this%mueff/ND+ND+4d0)
    Else
      this%covEvery = MAX(CEILING(ND/REAL(lambda, KIND=8)),1)
      If(PRESENT(covEvery)) this%covEvery = covEvery
      this%lastCov = 0
      info = this%reinit_X(params%X0, funW = funW, alphacovIn = alphacovIn)
      If(CheckF(info, 'mlf_cmaes_init: Error reinit_X')) RETURN
    Endif
  End Function mlf_cmaes_init

  Function mlf_cmaes_objcreate(ismaes, data_handler, params, funW, alphacovIn, covevery) Result(obj)
    type(mlf_maes_obj), pointer :: mthis
    type(mlf_cmaes_obj), pointer :: cthis
    class(mlf_optim_param), intent(inout) :: params
    class(mlf_obj), pointer :: obj
    class(mlf_weight_fun), optional, intent(in) :: funW
    real(c_double), intent(in), optional :: alphacovIn
    integer, intent(in), optional :: covevery
    class(mlf_data_handler), intent(inout), optional :: data_handler
    logical, intent(in) :: ismaes
    integer :: ret
    If(ismaes) Then
      ALLOCATE(mthis)
      ret = mthis%init(funW, alphacovIn, data_handler, params)
      If(ret<0) Then
        DEALLOCATE(mthis)
      Else
        obj => mthis
      Endif
    Else
      ALLOCATE(cthis)
      ret = cthis%init(funW, alphacovIn, covEvery, data_handler, params)
      If(ret<0) Then
        DEALLOCATE(cthis)
      Else
        obj => cthis
      Endif
    Endif
  End Function mlf_cmaes_objcreate

  Integer Function mlf_cmaes_reinit_X(this, X0, sigma0, funW, alphacovIn) Result(info)
    class(mlf_cmaes_obj), intent(inout), target :: this
    class(mlf_weight_fun), optional, intent(in) :: funW
    real(c_double), intent(in), optional :: X0(:), sigma0, alphacovIn
    integer ND
    info = mlf_matrix_es_obj_reinit_X(this, X0, sigma0, funW, alphacovIn)
    If(CheckF(info, 'mlf_cmaes_reinit_X: Error init matrix_es_obj')) RETURN
    ND = SIZE(this%Z, 1)
    CALL IdentityMatrix(this%C)
    this%pc = 0
    this%cp = (this%mueff/ND+4.0D0)/(2.0D0*this%mueff/ND+ND+4d0)
  End Function mlf_cmaes_reinit_X

  Integer Function mlf_matrix_es_obj_updateY(this, Y) Result(idMin)
    class(mlf_matrix_es_obj), intent(inout), target :: this
    real(c_double), intent(in) :: Y(:,:)
    integer :: lambda, ND, info, i
    lambda = SIZE(this%Z,2); ND = SIZE(this%Z,1)
    ASSOCIATE(mu=> this%mu, idx => this%idx, D => this%D, Z => this%Z, Dw => this%Dw, &
        W => this%W, sigma => this%sigma, mueff => this%mueff, c1 => this%c1, cs => this%cs, &
        cw => this%cw, ps => this%ps, dsigma => this%dsigma, chiN => this%chiN)
      ! Sort idx by the value Y
      CALL QSortIdx(Y, idx)
      PRINT *, "CMA opt sigma:",  sigma
      Do i=1, mu
        PRINT *, REAL(Y(:,idx(i)))
      End Do
      idMin = idx(1)
      ! Select the best values of D and Z
      D(:,:mu) = D(:,idx(:mu))
      Z(:,:mu) = Z(:,idx(:mu))
      ! Compute <d>_w
      Dw = MATMUL(D(:,:mu), W)
      this%X0 = this%X0 + sigma*Dw
      ! Update the value of ps using previous value and <Z>_w
      ps = (1d0-cs)*ps + SQRT(mueff*cs*(2d0-cs))*MATMUL(Z(:,:mu), W)

      info = this%updateMatrix()
      If(CheckF(info, 'mlf_matrix_es_obj_updateY: Error updateMatrix')) Then
        idMin = info
        RETURN
      Endif
      ! Evaluate new sigma
      sigma = sigma*EXP(cs/dsigma*(NORM2(ps)/chiN-1d0))
    END ASSOCIATE
  End Function mlf_matrix_es_obj_updateY

  Integer Function mlf_maes_updateMatrix(this) Result(info)
    class(mlf_maes_obj), intent(inout), target :: this
    integer :: ND, i
    info = 0; ND = SIZE(this%Z,1)
    ASSOCIATE(mu=> this%mu, Z => this%Z, Dw => this%Dw, ps => this%ps, &
        W => this%W, c1 => this%c1, cw => this%cw, DM => this%DM)
      ! Compute DM <d*d^T>_w
      ! The computaion is more precise using a zero matrix
      DM = 0
      Do i = mu,1,-1
        CALL TensorAddVectDown(ND, Z(:,i), DM, 0.5*cw*W(i)) ! Using only lower half is ~40% faster
      End Do
      CALL TensorAddVectDown(ND, ps, DM, 0.5*c1)
      FORALL(i = 1:ND) DM(i,i) = DM(i,i)+(1.0d0-0.5*(c1+cw))
      CALL SymmetrizeMatrix(DM)
      this%M = MATMUL(this%M, DM)
    END ASSOCIATE
  End Function mlf_maes_updateMatrix

  Integer Function mlf_cmaes_updateMatrix(this) result(info)
    class(mlf_cmaes_obj), intent(inout), target :: this
    integer :: ND, i
    info = 0; ND = SIZE(this%Z,1)
    ASSOCIATE(mu=> this%mu, D => this%D, Dw => this%Dw, &
        W => this%W, mueff => this%mueff, c1 => this%c1, C => this%C, &
        cw => this%cw, cp => this%cp, pc => this%pc, DM => this%DM)
      ! Update p using <d>_w
      pc = (1d0-cp)*pc + SQRT(mueff*cp*(2d0-cp))*Dw
      ! Compute DM <d*d^T>_w
      ! The computaion is more precise using a zero matrix
      DM = 0
      Do i = 1,mu
        CALL TensorAddVectDown(ND, D(:,i), DM, cw*W(i)) ! Using only lower half is ~40% faster
      End Do
      CALL TensorAddVectDown(ND, pc, DM, c1) ! Add rank one update
      CALL SymmetrizeMatrix(DM)
      ! Update covariance matrix
      C = (1.0d0-c1-cw)*C + DM
      this%lastCov = this%lastCov + 1
      If(this%lastCov >= this%covEvery) Then
        info = SqrtSymMatrix(C, this%M)
        If(CheckF(info, "mlf_cmaes_updateMatrix: error SqrtSymMatrix")) RETURN
        this%lastCov = 0
      Endif
    END ASSOCIATE
  End Function mlf_cmaes_updateMatrix

End Module mlf_cmaes

