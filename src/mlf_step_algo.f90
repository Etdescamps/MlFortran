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
  IMPLICIT NONE
  PRIVATE

  PUBLIC :: mlf_step_obj_init, mlf_step_obj_reinit

  ! Object handling timer and step evaluations
  Type, Public, extends(mlf_obj_model), abstract :: mlf_step_obj
    integer(kind=8) :: start_time
    real :: cpu_start
    integer(c_int64_t), pointer :: niter, nitermax
    real(c_double), pointer :: cputime, realtime, cpumax, realmax
  Contains
    procedure :: start_timer => mlf_step_start_timer
    procedure :: stop_timer => mlf_step_stop_timer
    procedure :: step =>  mlf_step_f
    procedure :: constraints_step =>  mlf_step_constraints
    procedure :: reinit => mlf_step_obj_reinit
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
  ! Init function for the step object
  integer Function mlf_step_obj_init(this, nipar, nrpar, nrsc, ifields, rfields, data_handler) result(info)
    class(mlf_step_obj), intent(inout), target :: this
    integer(c_int64_t), intent(inout) :: nipar, nrpar
    integer, intent(inout) :: nrsc
    integer, parameter :: ni = 2, nr = 4
    character(len=*,kind=c_char) :: ifields, rfields
    class(mlf_data_handler), intent(inout), optional :: data_handler
    nipar = nipar+ni; nrpar = nrpar+nr
    info = mlf_arr_init(this, nipar, nrpar, nrsc, C_CHAR_"niter;nitermax;"//ifields, &
      C_CHAR_"cputime;realtime;cpumax;realmax;"//rfields, data_handler)
    if(info < 0) RETURN
    this%niter => this%ipar(nipar+1)
    this%nitermax => this%ipar(nipar+2)
    this%cputime => this%rpar(nrpar+1)
    this%realtime => this%rpar(nrpar+2)
    this%cpumax => this%rpar(nrpar+3)
    this%realmax => this%rpar(nrpar+4)
    nipar = nipar+ni; nrpar = nrpar+nr
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
