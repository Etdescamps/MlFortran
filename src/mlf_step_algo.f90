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
  Integer, Parameter, Public :: mlf_STEP_STOP = 1, mlf_STEP_OK = 0

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
    If(present(nIPar)) this%nIPar = this%nIPar+nIPar
    If(present(nRPar)) this%nRPar = this%nRPar+nRPar
    If(present(nRsc)) this%nRsc = this%nRsc+nRsc
    If(present(nIVar)) this%nIVar = this%nIVar+nIVar
    If(present(nRVar)) this%nRVar = this%nRVar+nRVar
  End Subroutine mlf_step_addFields

  Subroutine mlf_step_initFields(this, nIPar, nRPar, nRsc, nIVar, nRVar)
    class(mlf_step_numFields), intent(inout) :: this
    integer, intent(in), optional :: nIPar, nRPar, nRsc, nIVar, nRVar
    this%nIPar = 0; this%nRPar = 0; this%nRsc = 0
    this%nIVar = 0; this%nRVar = 0
    If(present(nIPar)) this%nIPar = nIPar
    If(present(nRPar)) this%nRPar = nRPar
    If(present(nRsc)) this%nRsc = nRsc
    If(present(nIVar)) this%nIVar = nIVar
    If(present(nRVar)) this%nRVar = nRVar
  End Subroutine mlf_step_initFields

  Subroutine mlf_step_obj_addRVar(this, numFields, pnt, field)
    class(mlf_step_obj), intent(inout), target :: this
    class(mlf_step_numFields), intent(inout) :: numFields
    real(c_double), intent(out), pointer :: pnt
    character(len=*,kind=c_char), optional :: field
    numFields%nRVar = numFields%nRVar+1
    pnt => this%rVar(numFields%nRVar)
    If(PRESENT(field)) CALL this%v(this%idRVar)%addField(field)
  End Subroutine mlf_step_obj_addRVar

  Subroutine mlf_step_obj_addIVar(this, numFields, pnt, field)
    class(mlf_step_obj), intent(inout), target :: this
    class(mlf_step_numFields), intent(inout) :: numFields
    integer(c_int64_t), intent(out), pointer :: pnt
    character(len=*,kind=c_char), optional :: field
    numFields%nIVar = numFields%nIVar+1
    pnt => this%iVar(numFields%nIVar)
    If(PRESENT(field)) CALL this%v(this%idIVar)%addField(field)
  End Subroutine mlf_step_obj_addIVar
  
  ! Init function for the step object
  Integer Function mlf_step_obj_init(this, numFields, data_handler) Result(info)
    class(mlf_step_obj), intent(inout), target :: this
    class(mlf_step_numFields), intent(inout) :: numFields
    class(mlf_data_handler), intent(inout), optional :: data_handler
    CALL numFields%addFields(nIPar = 1, nRPar = 2, nRsc = 2, nIVar = 1, nRVar = 2)
    info = mlf_arr_init(this, numFields, data_handler)
    If(info < 0) RETURN
    info = this%add_i64array(numFields, numFields%niVar, this%iVar, C_CHAR_"iVar", &
      data_handler = data_handler, fixed_dims = [.TRUE.])
    If(CheckF(info, "Error creating iVar")) RETURN
    this%idIVar = numFields%nRsc
    info = this%add_rarray(numFields, numFields%nrVar, this%rVar, C_CHAR_"rVar", &
      data_handler = data_handler, fixed_dims = [.TRUE.])
    If(CheckF(info, "Error creating rVar")) RETURN
    this%idRVar = numFields%nRsc
    numFields%nIVar = 0; numFields%nRVar = 0 ! Reinit these fields to zero
    CALL this%addIVar(numFields, this%niter, "niter")
    CALL this%addIPar(numFields, this%niterMax, "niterMax")
    CALL this%addRVar(numFields, this%cpuTime, "cpuTime")
    CALL this%addRVar(numFields, this%realTime, "realTime")
    CALL this%addRPar(numFields, this%cpuMax, "cpuMax")
    CALL this%addRPar(numFields, this%realMax, "realMax")
  End Function mlf_step_obj_init

  Integer Function mlf_step_obj_reinit(this) Result(info)
    class(mlf_step_obj), intent(inout), target :: this
    this%niter = 0
    this%nitermax = HUGE(this%nitermax)
    this%cputime = 0d0
    this%realtime = 0d0
    this%realmax = HUGE(0d0)
    this%cpumax = HUGE(0d0)
    info = 0
  End Function mlf_step_obj_reinit

  Subroutine mlf_step_constraints(this, niterMax, realMax, cpuMax)
    class(mlf_step_obj), intent(inout), target :: this
    integer(c_int64_t), optional :: niterMax
    real(c_double), optional :: realMax, cpuMax
    If(PRESENT(niterMax)) this%nitermax = niterMax
    If(PRESENT(realMax)) this%realmax = realMax
    If(PRESENT(cpuMax)) this%cpumax = cpuMax
  End Subroutine mlf_step_constraints

  Real(c_double) Function mlf_step_stop_timer(this, niter) Result(t)
    class(mlf_step_obj), intent(inout), target :: this
    integer(kind=8) :: eltime, clock_rate
    integer(kind=8), optional :: niter
    real :: cputime
    CALL SYSTEM_CLOCK(eltime, clock_rate)
    CALL CPU_TIME(cputime)
    t = REAL(eltime-this%start_time,8)/REAL(clock_rate,8)
    this%realtime = this%realtime + t
    this%cputime = this%cputime + REAL(cputime-this%cpu_start, 8)
    If(PRESENT(niter)) Then
      this%niter = this%niter+niter
    Else
      this%niter = this%niter+1
    Endif
  End Function mlf_step_stop_timer
  
  Subroutine mlf_step_start_timer(this)
    class(mlf_step_obj), intent(inout), target :: this
    CALL CPU_TIME(this%cpu_start)
    CALL SYSTEM_CLOCK(this%start_time)
  End Subroutine mlf_step_start_timer

  Integer Function mlf_step_f(this, dt, niter) result(info)
    class(mlf_step_obj), intent(inout), target :: this
    real(c_double), intent(out), optional :: dt
    integer(kind=8), intent(inout), optional :: niter
    real(c_double) :: dt0
    CALL this%start_timer()
    info = this%stepF(niter)
    dt0 = this%stop_timer(niter)
    If(PRESENT(dt)) dt0 = dt
    If(info /= mlf_STEP_OK) RETURN
    If(this%cputime > this%cpumax .OR. this%realtime > this%realmax) Then
      info = mlf_STEP_STOP
      RETURN
    Endif
  End Function mlf_step_f

  ! C interface to step function
  Integer(c_int) Function mlf_step_c(cptr, dt, niter) Result(info) bind(C, name="mlf_step")
    type(c_ptr), value :: cptr
    real(c_double), intent(out) :: dt
    integer(c_int64_t), intent(inout) :: niter
    integer(kind=8) :: niterX
    class(mlf_obj), pointer :: obj
    info = -1
    niterX = niter
    obj => mlf_getobjfromc(cptr)
    If(.NOT. ASSOCIATED(obj)) RETURN
    SELECT TYPE(obj)
      CLASS is (mlf_step_obj)
        info = obj%step(dt, niterX)
        niter = niterX
    END SELECT
  End Function mlf_step_c
End Module mlf_step_algo
