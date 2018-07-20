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

Module mlf_step_algo
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_models
  Use mlf_errors
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: mlf_step_obj_init, mlf_step_obj_reinit

  Type, Public, extends(mlf_rsc_numFields) :: mlf_step_numFields
    integer(c_int64_t) :: nIVar
    integer(c_int64_t) :: nRVar
  Contains
    procedure :: addFields => mlf_step_addFields
    procedure :: initFields => mlf_step_initFields
  End Type mlf_step_numFields

  ! Object handling timer and step evaluations
  Type, Public, extends(mlf_obj_model), abstract :: mlf_step_obj
    integer(kind=8) :: start_time
    real :: cpu_start
    integer :: idRVar, idIVar
    real(c_double), pointer :: rVar(:)
    integer(c_int64_t), pointer :: niter, niterMax, iVar(:)
    real(c_double), pointer :: cpuTime, realTime, cpuMax, realMax
  Contains
    procedure :: start_timer => mlf_step_start_timer
    procedure :: stop_timer => mlf_step_stop_timer
    procedure :: step =>  mlf_step_f
    procedure :: constraints_step =>  mlf_step_constraints
    procedure :: reinit => mlf_step_obj_reinit
    procedure :: addRVar => mlf_step_obj_addRVar
    procedure :: addIVar => mlf_step_obj_addIVar
    procedure(mlf_step_fun), deferred :: stepF
  End Type mlf_step_obj

  Abstract Interface
    Integer Function mlf_step_fun(this, niter)
      import :: mlf_step_obj
      class(mlf_step_obj), intent(inout), target :: this
      integer(kind=8), intent(inout), optional :: niter
    End Function mlf_step_fun
  End Interface
Contains

  Subroutine mlf_step_addFields(this, nIPar, nRPar, nRsc, nIVar, nRVar)
    class(mlf_step_numFields), intent(inout) :: this
    integer, intent(in), optional :: nIPar, nRPar, nRsc, nIVar, nRVar
    if(present(nIPar)) this%nIPar = this%nIPar+nIPar
    if(present(nRPar)) this%nRPar = this%nRPar+nRPar
    if(present(nRsc)) this%nRsc = this%nRsc+nRsc
    if(present(nIVar)) this%nIVar = this%nIVar+nIVar
    if(present(nRVar)) this%nRVar = this%nRVar+nRVar
  End Subroutine mlf_step_addFields

  Subroutine mlf_step_initFields(this, nIPar, nRPar, nRsc, nIVar, nRVar)
    class(mlf_step_numFields), intent(inout) :: this
    integer, intent(in), optional :: nIPar, nRPar, nRsc, nIVar, nRVar
    this%nIPar = 0; this%nRPar = 0; this%nRsc = 0
    this%nIVar = 0; this%nRVar = 0
    if(present(nIPar)) this%nIPar = nIPar
    if(present(nRPar)) this%nRPar = nRPar
    if(present(nRsc)) this%nRsc = nRsc
    if(present(nIVar)) this%nIVar = nIVar
    if(present(nRVar)) this%nRVar = nRVar
  End Subroutine mlf_step_initFields

  Subroutine mlf_step_obj_addRVar(this, numFields, pnt, field)
    class(mlf_step_obj), intent(inout), target :: this
    class(mlf_step_numFields), intent(inout) :: numFields
    real(c_double), intent(out), pointer :: pnt
    character(len=*,kind=c_char), optional :: field
    numFields%nRVar = numFields%nRVar+1
    pnt => this%rVar(numFields%nRVar)
    if(present(field)) call this%v(this%idRVar)%addField(field)
  End Subroutine mlf_step_obj_addRVar

  Subroutine mlf_step_obj_addIVar(this, numFields, pnt, field)
    class(mlf_step_obj), intent(inout), target :: this
    class(mlf_step_numFields), intent(inout) :: numFields
    integer(c_int64_t), intent(out), pointer :: pnt
    character(len=*,kind=c_char), optional :: field
    numFields%nIVar = numFields%nIVar+1
    pnt => this%iVar(numFields%nIVar)
    if(present(field)) call this%v(this%idIVar)%addField(field)
  End Subroutine mlf_step_obj_addIVar
  
  ! Init function for the step object
  integer Function mlf_step_obj_init(this, numFields, data_handler) result(info)
    class(mlf_step_obj), intent(inout), target :: this
    class(mlf_step_numFields), intent(inout) :: numFields
    class(mlf_data_handler), intent(inout), optional :: data_handler
    call numFields%addFields(nIPar = 1, nRPar = 2, nRsc = 2, nIVar = 1, nRVar = 2)
    info = mlf_arr_init(this, numFields, data_handler)
    if(info < 0) RETURN
    info = this%add_i64array(numFields, numFields%niVar, this%iVar, C_CHAR_"iVar", &
      data_handler = data_handler, fixed_dims = [.TRUE.])
    if(CheckF(info, "Error creating iVar")) RETURN
    this%idIVar = numFields%nRsc
    info = this%add_rarray(numFields, numFields%nrVar, this%rVar, C_CHAR_"rVar", &
      data_handler = data_handler, fixed_dims = [.TRUE.])
    if(CheckF(info, "Error creating rVar")) RETURN
    this%idRVar = numFields%nRsc
    numFields%nIVar = 0; numFields%nRVar = 0 ! Reinit these fields to zero
    call this%addIVar(numFields, this%niter, "niter")
    call this%addIPar(numFields, this%niterMax, "niterMax")
    call this%addRVar(numFields, this%cpuTime, "cpuTime")
    call this%addRVar(numFields, this%realTime, "realTime")
    call this%addRPar(numFields, this%cpuMax, "cpuMax")
    call this%addRPar(numFields, this%realMax, "realMax")
  End Function mlf_step_obj_init

  integer Function mlf_step_obj_reinit(this) result(info)
    class(mlf_step_obj), intent(inout), target :: this
    this%niter = 0
    this%nitermax = huge(this%nitermax)
    this%cputime = 0d0
    this%realtime = 0d0
    this%realmax = huge(0d0)
    this%cpumax = huge(0d0)
    info = 0
  End Function mlf_step_obj_reinit

  Subroutine mlf_step_constraints(this, niterMax, realMax, cpuMax)
    class(mlf_step_obj), intent(inout), target :: this
    integer(c_int64_t), optional :: niterMax
    real(c_double), optional :: realMax, cpuMax
    if(present(niterMax)) this%nitermax = niterMax
    if(present(realMax)) this%realmax = realMax
    if(present(cpuMax)) this%cpumax = cpuMax
  End Subroutine mlf_step_constraints

  real(c_double) Function mlf_step_stop_timer(this, niter) result(t)
    class(mlf_step_obj), intent(inout), target :: this
    integer(kind=8) :: eltime, clock_rate
    integer(kind=8), optional :: niter
    real :: cputime
    call system_clock(eltime, clock_rate)
    call cpu_time(cputime)
    t = real(eltime-this%start_time,8)/real(clock_rate,8)
    this%realtime = this%realtime + t
    this%cputime = this%cputime + real(cputime-this%cpu_start, 8)
    if(present(niter)) then
      this%niter = this%niter+niter
    else
      this%niter = this%niter+1
    endif
  End Function mlf_step_stop_timer
  
  Subroutine mlf_step_start_timer(this)
    class(mlf_step_obj), intent(inout), target :: this
    call cpu_time(this%cpu_start)
    call system_clock(this%start_time)
  End Subroutine mlf_step_start_timer

  integer Function mlf_step_f(this, dt, niter) result(info)
    class(mlf_step_obj), intent(inout), target :: this
    real(c_double), intent(out), optional :: dt
    integer(kind=8), intent(inout), optional :: niter
    real(c_double) :: dt0
    call this%start_timer()
    info = this%stepF(niter)
    dt0 = this%stop_timer(niter)
    if(present(dt)) dt0 = dt
    if(info < 0) RETURN
    if(this%cputime > this%cpumax .OR. this%realtime > this%realmax) then
      info = 1
      return
    endif
  End Function mlf_step_f

  ! C interface to step function
  integer(c_int) Function mlf_step_c(cptr, dt, niter) result(info) bind(C, name="mlf_step")
    type(c_ptr), value :: cptr
    real(c_double), intent(out) :: dt
    integer(c_int64_t), intent(inout) :: niter
    integer(kind=8) :: niterX
    class(mlf_obj), pointer :: obj
    info = -1
    niterX = niter
    obj => mlf_getobjfromc(cptr)
    if(.NOT. associated(obj)) RETURN
    select type(obj)
      class is (mlf_step_obj)
        info = obj%step(dt, niterX)
        niter = niterX
    end select
  End Function mlf_step_c
End Module mlf_step_algo
