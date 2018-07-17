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

Module mlf_optim
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_step_algo
  Use mlf_fun_intf
  Use mlf_errors
  Use mlf_utils
  Use mlf_rand
  IMPLICIT NONE
  PRIVATE

  Public :: mlf_optim_init, mlf_optim_reinit, mlf_optim_stopCond

  Type, Public :: mlf_optim_param
    integer(c_int) :: lambda, mu
    integer(c_int64_t) :: nevalFunMax
    real(c_double) :: sigma, targetFun
    real(c_double), pointer :: XMin(:), XMax(:), X0(:)
    class(mlf_objective_fun), pointer :: fun
  Contains
    procedure :: init => mlf_optim_setparams
  End Type mlf_optim_param
  
  Type, Public, abstract, extends(mlf_step_obj) :: mlf_optim_obj
    class(mlf_objective_fun), pointer :: fun
    integer(c_int64_t), pointer :: nevalFun, nevalFunMax
    real(c_double), pointer :: minFun, targetFun, sigma, minX(:)
    real(c_double), allocatable :: X(:,:), Csr(:,:), Y(:,:)
    real(c_double), allocatable :: XInitMin(:), XInitMax(:)
    real(c_double) :: sigma0
    integer(c_int) :: lambda, mu
    integer, allocatable ::  idx(:)
  Contains
    procedure :: stepF => mlf_optim_stepF
    procedure :: stopCond => mlf_optim_stopCond
    procedure :: constraints_optim =>  mlf_optim_constraints
    procedure :: randomStartPoint => mlf_optim_startpoint
    procedure (mlf_genX), deferred :: genX
    procedure (mlf_updateY), deferred :: updateY
  End Type mlf_optim_obj

  Abstract Interface
    integer Function mlf_genX(this, ids)
      Use iso_c_binding
      import :: mlf_optim_obj
      class(mlf_optim_obj), intent(inout), target :: this
      integer(c_int), intent(in), optional :: ids(:) ! Regenerates selected ids
    End Function mlf_genX
    Function mlf_updateY(this, Y, idMin)
      Use iso_c_binding
      import :: mlf_optim_obj
      class(mlf_optim_obj), intent(inout), target :: this
      real(c_double), intent(in) :: Y(:,:)
      integer, intent(out), optional :: idMin
      real(c_double) :: mlf_updateY
    End Function mlf_updateY
  End Interface
Contains
  integer Function mlf_optim_init(this, nipar, nrpar, nrsc, ifields, rfields, data_handler, params) &
      result(info)
    class(mlf_optim_obj), intent(inout), target :: this
    class(mlf_optim_param), intent(inout) :: params
    integer(c_int64_t), intent(inout) :: nipar, nrpar
    integer, intent(inout) :: nrsc
    integer, parameter :: ni = 2, nr = 3, ns = 1
    integer(c_int64_t) :: nD
    character(len=*,kind=c_char) :: ifields, rfields
    class(mlf_data_handler), intent(inout), optional :: data_handler
    info = -1
    if(.NOT. ASSOCIATED(params%fun)) RETURN
    this%fun => params%fun
    nD = this%fun%nD
    nipar = nipar+ni; nrpar = nrpar+nr; nrsc = nrsc+ns
    info = mlf_step_obj_init(this, nipar, nrpar, nrsc, C_CHAR_"nevalFun;nevalFunMax;"//ifields,&
      C_CHAR_"minfun;targetFun;sigma;"//rfields, data_handler)
    if(info < 0) RETURN
    this%nevalFun => this%ipar(nipar+1)
    this%nevalFunMax => this%ipar(nipar+2)
    this%minFun => this%rpar(nrpar+1)
    this%targetFun => this%rpar(nrpar+2)
    this%sigma => this%rpar(nrpar+3)
    this%targetFun = params%targetFun
    this%sigma0 = params%sigma
    this%lambda = params%lambda
    if(params%mu <0) then
      this%mu = MAX(this%lambda/2, MIN(this%lambda, 2))
    else
      this%mu = params%mu
    endif
    this%nevalFunMax = params%nevalFunMax
    info = this%add_rarray(nrsc+1, nD, this%minX, C_CHAR_"minX", &
      data_handler = data_handler, fixed_dims = [.TRUE.])
    if(info < 0) RETURN
    nipar = nipar+ni; nrpar = nrpar+nr; nrsc = nrsc+ns
    ALLOCATE(this%X(ND, params%lambda), this%Y(this%fun%nY,params%lambda), this%idx(params%lambda))
    if(this%fun%nC>0) ALLOCATE(this%Csr(this%fun%nC, params%lambda))
    if(associated(params%XMin)) ALLOCATE(this%XInitMin, source=params%XMin)
    if(associated(params%XMax)) ALLOCATE(this%XInitMin, source=params%XMax)
  End Function mlf_optim_init

  integer Function mlf_optim_reinit(this) result(info)
    class(mlf_optim_obj), intent(inout), target :: this
    info = mlf_step_obj_reinit(this)
    this%minFun = HUGE(this%targetFun)
    this%nevalFun = 0
    this%sigma = this%sigma0
  End Function mlf_optim_reinit

  Subroutine mlf_optim_constraints(this, targetFun, nevalFunMax)
    class(mlf_optim_obj), intent(inout), target :: this
    integer(c_int64_t), optional :: nevalFunMax
    real(c_double), optional :: targetFun
    if(present(targetFun)) this%targetFun = targetFun
    if(present(nevalFunMax)) this%nevalFunMax = nevalFunMax
  End Subroutine mlf_optim_constraints

  integer Function mlf_optim_stopCond(this) result(info)
    class(mlf_optim_obj), intent(inout), target :: this
    if(this%minFun < this%targetFun) then
      info = 1
      RETURN
    endif
    if(this%nevalFun > this%nevalFunMax) then
      info = 1
      RETURN
    endif
    info = 0
  End Function mlf_optim_stopCond

  Integer Function mlf_optim_stepF(this, niter) result(info)
    class(mlf_optim_obj), intent(inout), target :: this
    integer(kind=8), intent(inout), optional :: niter
    integer(kind=8) :: i, niter0=1
    integer :: lambda, nR, nR2, j, idMin
    logical, allocatable :: VU(:)
    real(c_double) :: minv, minc
    minv = huge(minv); minc = huge(minc)
    lambda = size(this%X,2)
    info = -1
    if(present(niter)) niter0 = niter
    ALLOCATE(VU(lambda))
    LIter: Do i=1,niter0
      info = this%stopCond()
      if(info /= 0) EXIT LIter
      info = this%genX()
      if(info<0) EXIT LIter
      if(allocated(this%Csr)) then
        info = this%fun%constraints(this%X, this%Csr)
        if(info<0) EXIT LIter
        VU = any(this%Csr>0, dim = 1)
        nR = count(VU)
        if(nR == lambda .OR. (nR == lambda-1 .AND. nR > 2)) then
          minC = this%updateY(this%Csr)
          CYCLE LIter
        endif
        if(nR > 0) then
          this%idx = [(j, j=1,lambda)]
          this%idx(:nR) = pack(this%idx, VU)
          LConst: Do While(nR>0)
            info = this%genX(this%idx(:nR))
            if(info<0) EXIT LIter
            info = this%fun%constraints(this%X(:,this%idx(:nR)), this%Csr(:,:nR))
            if(info<0) EXIT LIter
            VU(:nR) = any(this%Csr(:,:nR)>0, dim = 1)
            nR2 = count(VU(:nR))
            if(nR2 == 0) EXIT LConst
            this%idx(:nR2) = pack(this%idx(:nR), VU(:nR))
            nR = nR2
          End Do LConst
        endif
      endif
      info = this%fun%eval(this%X, this%Y)
      if(info<0) EXIT LIter
      WHERE(ieee_is_nan(this%Y)) this%Y = huge(0d0) 
      minV = this%updateY(this%Y, idMin)
      if(idMin<0) then
        info = idMin
        EXIT LIter
      endif
      if(minV < this%minFun) then
        this%minFun = minV
        this%minX = this%X(:,idMin)
      endif
      this%nevalFun = this%nevalFun + lambda
      if(present(niter)) niter = i
    End Do LIter
  End Function mlf_optim_stepF

  Subroutine mlf_optim_startpoint(this, X)
    class(mlf_optim_obj), intent(inout), target :: this
    real(c_double), intent(out) :: X(:)
    if(allocated(this%XInitMin) .AND. allocated(this%XInitMax)) then
      call random_number(X)
      X = this%XInitMin+X*this%XInitMax
    else
      call randN(X)
      X = this%sigma0*X
    endif
  End Subroutine mlf_optim_startpoint

  Subroutine mlf_optim_setparams(p, fun, lambda, mu, XMin, XMax, X0, nevalFunMax, sigma, targetFun)
    class(mlf_optim_param), intent(inout) :: p
    class(mlf_objective_fun), target :: fun
    integer(c_int) :: lambda
    integer(c_int), optional :: mu
    integer(c_int64_t), optional :: nevalFunMax
    real(c_double), optional :: sigma, targetFun
    real(c_double), optional, target :: XMin(:), XMax(:), X0(:)
    p%fun => fun
    p%lambda = lambda
    call SetIf(p%mu, -1, mu)
    call SetIf(p%nevalFunMax, HUGE(p%nevalFunMax), nevalFunMax)
    call SetIf(p%sigma, 1d0, sigma)
    call SetIf(p%targetFun, 1d0, targetFun)
    p%XMin => NULL(); p%XMax => NULL(); p%X0 => NULL()
    if(present(XMin)) p%XMin => XMin
    if(present(XMax)) p%XMax => XMax
    if(present(X0)) p%X0 => X0
  End Subroutine mlf_optim_setparams
End module mlf_optim

