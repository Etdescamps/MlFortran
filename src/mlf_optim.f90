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
  
  Type, Public, Abstract, Extends(mlf_step_obj) :: mlf_optim_obj
    class(mlf_objective_fun), pointer :: fun
    integer(c_int64_t), pointer :: nevalFun, nevalFunMax
    real(c_double), pointer :: minFun, targetFun, sigma, X(:,:), minX(:)
    real(c_double), allocatable :: Y(:,:)
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
    Integer Function mlf_genX(this, ids)
      Use iso_c_binding
      import :: mlf_optim_obj
      class(mlf_optim_obj), intent(inout), target :: this
      integer(c_int), intent(in), optional :: ids(:) ! Regenerates selected ids
    End Function mlf_genX

    Integer Function mlf_updateY(this, Y)
      Use iso_c_binding
      import :: mlf_optim_obj
      class(mlf_optim_obj), intent(inout), target :: this
      real(c_double), intent(in) :: Y(:,:)
    End Function mlf_updateY
  End Interface
Contains
  Integer Function mlf_optim_init(this, numFields, data_handler, params) &
      Result(info)
    class(mlf_optim_obj), intent(inout), target :: this
    class(mlf_step_numFields), intent(inout) :: numFields
    class(mlf_data_handler), intent(inout), optional :: data_handler
    class(mlf_optim_param), intent(inout) :: params
    integer(c_int64_t) :: nD, nLambdaND(2), nY
    info = -1
    If(.NOT. ASSOCIATED(params%fun)) RETURN
    this%fun => params%fun
    nD = this%fun%nD
    CALL numFields%addFields(nIPar = 1, nRPar = 1, nRsc = 2, nIVar = 1, nRVar = 2)
    info = mlf_step_obj_init(this, numFields, data_handler)
    If(info < 0) RETURN
    CALL this%addIVar(numFields, this%nevalFun, "nevalFun")
    CALL this%addIPar(numFields, this%nevalFunMax, "nevalFunMax")
    CALL this%addRVar(numFields, this%minFun, "minFun")
    CALL this%addRVar(numFields, this%sigma, "sigma")
    CALL this%addRPar(numFields, this%targetFun, "targetFun")
    this%targetFun = params%targetFun
    this%sigma0 = params%sigma
    this%lambda = params%lambda
    this%nevalFunMax = params%nevalFunMax
    info = this%add_rarray(numFields, nD, this%minX, C_CHAR_"minX", &
      data_handler = data_handler, fixed_dims = [.TRUE.])
    If(info < 0) RETURN
    nLambdaND = [ND, INT(this%lambda, KIND=c_int64_t)]
    info = this%add_rmatrix(numFields, nLambdaND, this%X, C_CHAR_"X", &
      data_handler = data_handler, fixed_dims = [.TRUE., .TRUE.])
    this%lambda = INT(nLambdaND(2), KIND=4)
    If(params%mu <0 .AND. this%lambda > 0) Then
      this%mu = MAX(this%lambda/2, MIN(this%lambda, 2))
    Else
      this%mu = params%mu
    Endif
    nY = this%fun%nY+this%fun%nC
    ALLOCATE(this%Y(nY,params%lambda), this%idx(params%lambda))
    If(ASSOCIATED(params%XMin)) ALLOCATE(this%XInitMin, source=params%XMin)
    If(ASSOCIATED(params%XMax)) ALLOCATE(this%XInitMax, source=params%XMax)
  End Function mlf_optim_init

  Integer Function mlf_optim_reinit(this) Result(info)
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
    If(PRESENT(targetFun)) this%targetFun = targetFun
    If(PRESENT(nevalFunMax)) this%nevalFunMax = nevalFunMax
  End Subroutine mlf_optim_constraints

  Integer Function mlf_optim_stopCond(this) Result(info)
    class(mlf_optim_obj), intent(inout), target :: this
    If(this%minFun < this%targetFun) Then
      info = 1
      RETURN
    Endif
    If(this%nevalFun > this%nevalFunMax) Then
      info = 1
      RETURN
    Endif
    info = 0
  End Function mlf_optim_stopCond

  Integer Function mlf_optim_stepF(this, niter) Result(info)
    class(mlf_optim_obj), intent(inout), target :: this
    integer(kind=8), intent(inout), optional :: niter
    integer(kind=8) :: i, niter0
    integer :: lambda, idMin
    logical, allocatable :: VU(:)
    real(c_double) :: minv, minc
    minv = HUGE(minv); minc = HUGE(minc)
    lambda = SIZE(this%X,2)
    info = 0
    CALL InitOrDefault(niter0, 1_8, niter)
    ALLOCATE(VU(lambda))
    Do i=1,niter0
      info = this%stopCond()
      If(info /= 0) RETURN
      info = this%genX()
      If(info<0) RETURN
      info = this%fun%eval(this%X, this%Y)
      If(info<0) RETURN
      WHERE(ieee_is_nan(this%Y)) this%Y = HUGE(0d0) 
      ! Y contains firstly the constraints and then the objectives
      idMin = this%updateY(this%Y)
      If(idMin<0) Then
        info = idMin
        EXIT
      Endif
      If(this%fun%nC > 0) Then
        If(ANY(this%Y(1:this%fun%nC, idMin) /= 0)) Then
          minV = HUGE(1d0)
        Else
          minV = this%Y(1+this%fun%nC, idMin)
        Endif
      Else
        minV = this%Y(1, idMin)
      Endif
      If(minV < this%minFun) Then
        this%minFun = minV
        this%minX = this%X(:,idMin)
      Endif
      this%nevalFun = this%nevalFun + lambda
      If(PRESENT(niter)) niter = i
    End Do
    If(info > 0) info = 0
  End Function mlf_optim_stepF

  Subroutine mlf_optim_startpoint(this, X)
    class(mlf_optim_obj), intent(inout), target :: this
    real(c_double), intent(out) :: X(:)
    If(ALLOCATED(this%XInitMin) .AND. ALLOCATED(this%XInitMax)) Then
      CALL RANDOM_NUMBER(X)
      X = this%XInitMin+X*this%XInitMax
    Else
      CALL randN(X)
      X = this%sigma0*X
    Endif
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
    CALL SetIf(p%mu, -1, mu)
    CALL SetIf(p%nevalFunMax, HUGE(p%nevalFunMax), nevalFunMax)
    CALL SetIf(p%sigma, 1d0, sigma)
    CALL SetIf(p%targetFun, 1d0, targetFun)
    p%XMin => NULL(); p%XMax => NULL(); p%X0 => NULL()
    If(PRESENT(XMin)) p%XMin => XMin
    If(PRESENT(XMax)) p%XMax => XMax
    If(PRESENT(X0)) p%X0 => X0
  End Subroutine mlf_optim_setparams
End module mlf_optim

