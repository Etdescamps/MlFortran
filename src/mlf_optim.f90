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
  IMPLICIT NONE
  PRIVATE

  Public :: mlf_optim_init, mlf_optim_reinit, mlf_optim_stopCond
  
  Type, Public, abstract, extends(mlf_step_obj) :: mlf_optim_obj
    class(mlf_objective_fun), allocatable :: fun
    integer(c_int64_t), pointer :: nevalFun, nevalFunMax
    real(c_double), pointer :: minFun, targetFun, minX(:)
    real(c_double), allocatable :: X(:,:), Csr(:,:), Y(:,:)
    integer, allocatable ::  idx(:)
  Contains
    procedure :: stepF => mlf_optim_stepF
    procedure :: stopCond => mlf_optim_stopCond
    procedure :: constraints_optim =>  mlf_optim_constraints
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
  integer Function mlf_optim_init(this, lambda, fun, nipar, nrpar, nrsc, ifields, rfields, data_handler) &
      result(info)
    class(mlf_optim_obj), intent(inout), target :: this
    class(mlf_objective_fun), intent(in) :: fun
    integer(c_int64_t), intent(inout) :: nipar, nrpar
    integer, intent(inout) :: nrsc
    integer, intent(in) :: lambda
    integer, parameter :: ni = 2, nr = 2, ns = 1
    integer(c_int64_t) :: nD
    character(len=*,kind=c_char) :: ifields, rfields
    class(mlf_data_handler), intent(inout), optional :: data_handler
    nD = fun%nD
    nipar = nipar+ni; nrpar = nrpar+nr; nrsc = nrsc+ns
    info = mlf_step_obj_init(this, nipar, nrpar, nrsc, C_CHAR_"nevalFun;nevalFunMax;"//ifields,&
      C_CHAR_"minfun;targetFun;"//rfields, data_handler)
    if(info < 0) RETURN
    this%nevalFun => this%ipar(nipar)
    this%nevalFunMax => this%ipar(nipar+1)
    this%minFun => this%rpar(nrpar)
    this%targetFun => this%rpar(nrpar+1)
    this%fun = fun
    info = this%add_rarray(nrsc, nD, this%minX, C_CHAR_"minX", &
      data_handler = data_handler, fixed_dims = [.TRUE.])
    if(info < 0) RETURN
    nipar = nipar+ni; nrpar = nrpar+nr; nrsc = nrsc+ns
    ALLOCATE(this%X(fun%ND, lambda), this%Y(fun%nY,lambda), this%idx(lambda))
    if(fun%nC>0) ALLOCATE(this%Csr(fun%nC, lambda))
  End Function mlf_optim_init

  integer Function mlf_optim_reinit(this) result(info)
    class(mlf_optim_obj), intent(inout), target :: this
    info = mlf_step_obj_reinit(this)
    this%targetFun = -HUGE(this%targetFun)
    this%minFun = HUGE(this%targetFun)
    this%nevalFun = 0
    this%nevalFunMax = HUGE(this%nevalFunMax)
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
    End Do LIter
    if(present(niter)) niter = i
  End Function mlf_optim_stepF
End module mlf_optim

