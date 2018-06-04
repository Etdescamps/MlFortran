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

  ! Common class for CMA-ES and MA-ES
  Type, Public, abstract, extends(mlf_optim_obj) :: mlf_matrix_es_obj
    real(c_double), pointer :: M(:,:) ! M square root matrix of covariance
    real(c_double), pointer :: X0(:) ! Center of the sampling: X(:,i) = D(:,i)+X0(:)
    real(c_double), pointer :: ps(:) ! Isotropically path (that comes from the best mu element of Z)
    real(c_double), pointer :: W(:) ! weight used by the methods on best mu element
    real(c_double), pointer :: sigma, alphacov
    real(c_double), allocatable :: Z(:,:) ! Vectors sampled using Norm(0,Id_N)
    real(c_double), allocatable :: D(:,:) ! Previous vectors multiplied by M (D = MZ)
    real(c_double), allocatable :: Dw(:), DM(:,:)
    real(c_double) :: chiN, mueff, cs, cw, c1, dsigma
    integer :: mu
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
    integer Function mlf_matrix_updateF(this)
      import :: mlf_matrix_es_obj
      class(mlf_matrix_es_obj), intent(inout), target :: this
    End Function mlf_matrix_updateF
  End Interface
Contains
  Integer Function mlf_matrix_es_obj_init(this, fun, nipar, nrpar, nrsc, ifields, rfields, &
    lambdaIn, muIn, data_handler) result(info)
    class(mlf_matrix_es_obj), intent(inout), target :: this
    class(mlf_objective_fun), intent(in) :: fun
    integer, intent(inout) :: nrsc
    integer(c_int64_t), intent(inout) :: nipar, nrpar
    character(len=*,kind=c_char) :: ifields, rfields
    integer, intent(in), optional :: lambdaIn, muIn
    class(mlf_data_handler), intent(inout), optional :: data_handler
    integer, parameter :: ni = 0, nr = 2, ns = 4
    integer :: lambda
    integer(c_int64_t) :: nD, CND(2), mu
    nD = fun%nD
    CND = nD
    if(present(lambdaIn)) then
      lambda = lambdaIn
    else
      lambda = 4 + floor(3.0*log(real(ND)))
    end if
    if(present(muIn)) then
      this%mu = muIn
    else
      this%mu = MAX(lambda/2,1)
    endif
    mu = this%mu
    nipar = nipar+ni; nrpar = nrpar+nr; nrsc = nrsc+ns
    info = mlf_optim_init(this, lambda, fun, nipar, nrpar, nrsc, ifields, &
      C_CHAR_"sigma;alphacov;"//rfields, data_handler)
    info = this%add_rmatrix(nrsc, CND, this%M, C_CHAR_"M", data_handler = data_handler, fixed_dims = [.TRUE., .TRUE.])
    info = this%add_rarray(nrsc+1, nD, this%X0, C_CHAR_"X0", data_handler = data_handler, fixed_dims = [.TRUE.])
    info = this%add_rarray(nrsc+2, nD, this%ps, C_CHAR_"ps", data_handler = data_handler, fixed_dims = [.TRUE.])
    info = this%add_rarray(nrsc+3, mu, this%W, C_CHAR_"W", data_handler = data_handler)
    this%mu = int(mu, kind=4)
    this%sigma => this%rpar(nrpar)
    this%alphacov => this%rpar(nrpar+1)
    nipar = nipar+ni; nrpar = nrpar+nr; nrsc = nrsc+ns
    ALLOCATE(this%Z(ND, lambda), this%D(ND, lambda), this%Dw(ND), this%DM(ND,ND))
  End Function mlf_matrix_es_obj_init

  Subroutine mlf_matrix_es_obj_updateW(this)
    class(mlf_matrix_es_obj), intent(inout), target :: this
    integer :: ND
    ND = size(this%Z, 1)
    Associate(W => this%W, mueff => this%mueff, alphacov => this%alphacov, &
        chiN => this%chiN, c1 => this%c1, cs => this%cs, cw => this%cw, &
        dsigma => this%dsigma, nR =>real(nD, kind=8))
      W = 1d0/sum(W) * W
      mueff = 1d0/(sum(W * W))
      c1 = alphacov/((nR+1.3d0)**2 + mueff)
      cs = (mueff + 2d0)/(mueff + nR + 5d0)
      cw = min(1d0 - c1, alphacov*(mueff+1d0/mueff-2d0)/((nR+2d0)**2 + alphacov*mueff*0.5d0))
      dsigma = 1d0 + cs + 2d0*max(0d0, sqrt((mueff-1d0)/(nR+1d0))-1d0)
      chiN = sqrt(nR)*(1d0-1d0/(4d0*nR)+1d0/(21d0*nR*nR))
    End Associate
  End Subroutine mlf_matrix_es_obj_updateW 

  integer Function mlf_matrix_es_obj_reinit(this) result(info)
    class(mlf_matrix_es_obj), intent(inout), target :: this
    info = this%reinit_X()
  End Function mlf_matrix_es_obj_reinit
 
  integer Function mlf_matrix_es_obj_reinit_X(this, X0, sigma0, funW, alphacovIn) result(info)
    class(mlf_matrix_es_obj), intent(inout), target :: this
    class(mlf_weight_fun), optional, intent(in) :: funW
    real(c_double), intent(in), optional :: X0(:), sigma0, alphacovIn
    integer :: lambda, mu, i
    info = mlf_optim_reinit(this)
    lambda = size(this%Z, 2); mu = size(this%W)
    call IdentityMatrix(this%M)
    this%ps = 0
    this%X0 = 0
    if(present(X0)) this%X0 = X0
    this%sigma = 1d0
    if(present(sigma0)) this%sigma = sigma0
    this%alphacov = 2d0
    if(present(alphacovIn)) this%alphacov=alphacovIn
    if(present(funW)) then
      call funW%eval(this%W, lambda)
    else
      this%W = log(0.5d0*(lambda+1d0))-[(log(real(i,kind=8)), i=1,mu)]
    endif
    call this%updateW()
  End Function mlf_matrix_es_obj_reinit_X

  integer Function mlf_matrix_es_obj_genX(this, ids) result(info)
    class(mlf_matrix_es_obj), intent(inout), target :: this
    integer(c_int), intent(in), optional :: ids(:) ! Regenerates selected ids
    integer :: i, j
    info = 0
    if(present(ids)) then
      Do i=1,size(ids)
        j = ids(i)
        call RandN(this%Z(:,j))
        this%D(:,j) = matmul(this%M, this%Z(:,j))
        this%X(:,j) = this%sigma*this%D(:,j)+this%X0
      End Do
    else
      call RandN(this%Z)
      Do i=1,size(this%Z,2)
        this%D(:,i) = matmul(this%M, this%Z(:,i))
        this%X(:,i) = this%sigma*this%D(:,i)+this%X0
      End Do
    endif
  End Function mlf_matrix_es_obj_genX
  
  integer Function mlf_matrix_es_obj_stopCond(this) result(info)
    class(mlf_matrix_es_obj), intent(inout), target :: this
    info = mlf_optim_stopCond(this)
  End Function mlf_matrix_es_obj_stopCond

  Integer Function mlf_maes_init(this, fun, X0, sigma0, funW, lambdaIn, muIn, &
      alphacovIn, data_handler) result(info)
    class(mlf_maes_obj), intent(inout), target :: this
    class(mlf_objective_fun), intent(in) :: fun
    integer(c_int64_t) :: nipar = 0, nrpar = 0
    integer :: nrsc = 0
    class(mlf_weight_fun), optional, intent(in) :: funW
    real(c_double), intent(in), optional :: X0(:), sigma0, alphacovIn
    integer, intent(in), optional :: lambdaIn, muIn
    class(mlf_data_handler), intent(inout), optional :: data_handler
    info = mlf_matrix_es_obj_init(this, fun, nipar, nrpar, nrsc, C_CHAR_"", &
      C_CHAR_"", lambdaIn, muIn, data_handler)
    if(CheckF(info, 'mlf_maes_init: Error init matrix_es')) RETURN
    if(present(data_handler)) then
      call this%updateW()
    else
      info = this%reinit_X(X0, sigma0, funW, alphacovIn)
    endif
  End Function mlf_maes_init

  Integer Function mlf_cmaes_init(this, fun, X0, sigma0, funW, lambdaIn, muIn, &
      alphacovIn, covEvery, data_handler) result(info)
    class(mlf_cmaes_obj), intent(inout), target :: this
    class(mlf_objective_fun), intent(in) :: fun
    integer(c_int64_t) :: nipar = 2, nrpar = 0, ND, CND(2)
    integer :: nrsc = 2
    class(mlf_weight_fun), optional, intent(in) :: funW
    real(c_double), intent(in), optional :: X0(:), sigma0, alphacovIn
    integer, intent(in), optional :: lambdaIn, muIn, covEvery
    class(mlf_data_handler), intent(inout), optional :: data_handler
    integer :: lambda
     
    info = mlf_matrix_es_obj_init(this, fun, nipar, nrpar, nrsc, C_CHAR_"lastCov;covEvery;", &
      C_CHAR_"", lambdaIn, muIn, data_handler)
    if(CheckF(info, 'mlf_cmaes_init: Error init matrix_es')) RETURN
    lambda = size(this%Z, 2); ND = size(this%Z, 1)
    CND = ND
    info = this%add_rmatrix(nrsc, CND, this%C, C_CHAR_"C", data_handler = data_handler, fixed_dims = [.TRUE., .TRUE.])
    info = this%add_rarray(nrsc+1, ND, this%pc, C_CHAR_"pc", data_handler = data_handler, fixed_dims = [.TRUE.])
    this%lastCov => this%ipar(nipar)
    this%covEvery => this%ipar(nipar+1)
    if(present(data_handler)) then
      call this%updateW()
      this%cp = (this%mueff/ND+4.0D0)/(2.0D0*this%mueff/ND+ND+4d0)
    else
      this%covEvery = MAX(CEILING(ND/REAL(lambda, kind=8)),1)
      if(present(covEvery)) this%covEvery = covEvery
      this%lastCov = 0
      info = this%reinit_X(X0, sigma0, funW, alphacovIn)
      if(CheckF(info, 'mlf_cmaes_init: Error reinit_X')) RETURN
    endif
  End Function mlf_cmaes_init

  integer Function mlf_cmaes_reinit_X(this, X0, sigma0, funW, alphacovIn) result(info)
    class(mlf_cmaes_obj), intent(inout), target :: this
    class(mlf_weight_fun), optional, intent(in) :: funW
    real(c_double), intent(in), optional :: X0(:), sigma0, alphacovIn
    integer ND
    info = mlf_matrix_es_obj_reinit_X(this, X0, sigma0, funW, alphacovIn)
    if(CheckF(info, 'mlf_cmaes_reinit_X: Error init matrix_es_obj')) RETURN
    ND = size(this%Z, 1)
    call IdentityMatrix(this%C)
    this%pc = 0
    this%cp = (this%mueff/ND+4.0D0)/(2.0D0*this%mueff/ND+ND+4d0)
  End Function mlf_cmaes_reinit_X

  real(c_double) Function mlf_matrix_es_obj_updateY(this, Y, idMin) result(Ymin)
    class(mlf_matrix_es_obj), intent(inout), target :: this
    real(c_double), intent(in) :: Y(:,:)
    integer, intent(out), optional :: idMin
    integer :: lambda, ND, info
    lambda = size(this%Z,2); ND = size(this%Z,1)
    Associate(mu=> this%mu, idx => this%idx, D => this%D, Z => this%Z, Dw => this%Dw, &
        W => this%W, sigma => this%sigma, mueff => this%mueff, c1 => this%c1, cs => this%cs, &
        cw => this%cw, ps => this%ps, dsigma => this%dsigma, chiN => this%chiN)
      ! Sort idx by the value Y
      call QSortIdx(Y, idx)
      Ymin = Y(1, idx(1))
      if(present(idMin)) idMin = idx(1)
      ! Select the best values of D and Z
      D(:,:mu) = D(:,idx(:mu))
      Z(:,:mu) = Z(:,idx(:mu))
      ! Compute <d>_w
      Dw = matmul(D(:,:mu), W)
      this%X0 = this%X0 + sigma*Dw
      ! Update the value of ps using previous value and <Z>_w
      ps = (1d0-cs)*ps + sqrt(mueff*cs*(2d0-cs))*matmul(Z(:,:mu), W)

      info = this%updateMatrix()
      if(CheckF(info, 'mlf_matrix_es_obj_updateY: Error updateMatrix')) RETURN

      ! Evaluate new sigma
      sigma = sigma*exp(cs/dsigma*(norm2(ps)/chiN-1d0))
    End Associate
  End Function mlf_matrix_es_obj_updateY

  Integer Function mlf_maes_updateMatrix(this) result(info)
    class(mlf_maes_obj), intent(inout), target :: this
    integer :: ND, i
    info = 0; ND = size(this%Z,1)
    Associate(mu=> this%mu, Z => this%Z, Dw => this%Dw, ps => this%ps, &
        W => this%W, c1 => this%c1, cw => this%cw, DM => this%DM)
      ! Compute DM <d*d^T>_w
      ! The computaion is more precise using a zero matrix
      DM = 0
      do i = mu,1,-1
        call TensorAddVectDown(ND, Z(:,i), DM, 0.5*cw*W(i)) ! Using only lower half is ~40% faster
      end do
      call TensorAddVectDown(ND, ps, DM, 0.5*c1)
      forall(i = 1:ND) DM(i,i) = DM(i,i)+(1.0d0-0.5*(c1+cw))
      call SymmetrizeMatrix(DM)
      this%M = matmul(this%M, DM)
    End Associate
  End Function mlf_maes_updateMatrix

  Integer Function mlf_cmaes_updateMatrix(this) result(info)
    class(mlf_cmaes_obj), intent(inout), target :: this
    integer :: ND, i
    info = 0; ND = size(this%Z,1)
    Associate(mu=> this%mu, D => this%D, Dw => this%Dw, &
        W => this%W, mueff => this%mueff, c1 => this%c1, C => this%C, &
        cw => this%cw, cp => this%cp, pc => this%pc, DM => this%DM)
      ! Update p using <d>_w
      pc = (1d0-cp)*pc + sqrt(mueff*cp*(2d0-cp))*Dw
      ! Compute DM <d*d^T>_w
      ! The computaion is more precise using a zero matrix
      DM = 0
      do i = 1,mu
        call TensorAddVectDown(ND, D(:,i), DM, cw*W(i)) ! Using only lower half is ~40% faster
      end do
      call TensorAddVectDown(ND, pc, DM, c1) ! Add rank one update
      call SymmetrizeMatrix(DM)
      ! Update covariance matrix
      C = (1.0d0-c1-cw)*C + DM
      this%lastCov = this%lastCov + 1
      if(this%lastCov >= this%covEvery) then
        info = SqrtSymMatrix(C, this%M)
        if(CheckF(info, "mlf_cmaes_updateMatrix: error SqrtSymMatrix")) RETURN
        this%lastCov = 0
      endif
    End Associate
  End Function mlf_cmaes_updateMatrix

End Module mlf_cmaes

