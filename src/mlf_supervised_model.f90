! Copyright (c) 2018-2019 Etienne Descamps
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

Module mlf_supervised_model
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_errors
  Use mlf_rand
  Use mlf_rsc_array
  Use mlf_models
  Use mlf_step_algo
  Use mlf_fun_intf
  Use mlf_utils
  IMPLICIT NONE
  PRIVATE
 
  Type, Public :: mlf_model_parameters
  End Type mlf_model_parameters

  Type, Public, Abstract, Extends(mlf_model_parameters) :: mlf_model_real_parameters
  Contains
    procedure(mlf_real_get_nparameters), deferred :: getNParameters
    procedure(mlf_real_parameters_set), deferred :: set
  End Type mlf_model_real_parameters

  Type, Public :: mlf_model_experiment
  End Type mlf_model_experiment

  Type, Public :: mlf_model_result
  End Type mlf_model_result

  Type, Public, Abstract, Extends(mlf_model_result) :: mlf_result_vectReal
  Contains
    procedure(mlf_vectResults), deferred :: vectResults
    procedure(mlf_real_get_noutput), deferred :: getNOutput
  End Type mlf_result_vectReal

  Type, Public, Abstract, Extends(mlf_step_obj) :: mlf_experience_model
  Contains
    procedure(mlf_model_setupModel), deferred :: setupModel
    procedure :: getResults => mlf_experience_dummy_results
    procedure :: getNCstr => mlf_experience_dummy_getNCstr
  End Type mlf_experience_model

  Type, Public :: mlf_data_point
    class(mlf_model_experiment), allocatable :: experiment
    class(mlf_model_result), allocatable :: results
    real(c_double) :: weight = 1d0
  End Type mlf_data_point

  Type, Public, Abstract, Extends(mlf_obj) :: mlf_experience_runner
    real(c_double), allocatable :: Z(:,:)
  Contains
    procedure(mlf_experience_run), deferred :: runExp
  End Type mlf_experience_runner

  Type, Public, Extends(mlf_experience_runner) :: mlf_local_experience_runner
    class(mlf_experience_model), allocatable :: model
    class(mlf_model_real_parameters), allocatable :: params
    class(mlf_result_vectReal), allocatable :: results
  Contains
    procedure :: runExp => mlf_local_runExp
  End Type mlf_local_experience_runner

  Type, Public, Abstract :: mlf_local_initializer
  Contains
    procedure(mlf_local_init_runner), deferred :: init_runner
  End Type mlf_local_initializer

  Abstract Interface
    Integer Function mlf_vectResults(this, Y)
      Use iso_c_binding
      import :: mlf_result_vectReal
      class(mlf_result_vectReal), intent(in), target :: this
      real(c_double), intent(out) :: Y(:)
    End Function mlf_vectResults

    Integer Function mlf_experience_run(this, experiment, X, Y, Cstr)
      Use iso_c_binding
      import :: mlf_experience_runner, mlf_model_experiment
      class(mlf_experience_runner), intent(inout), target :: this
      class(mlf_model_experiment), intent(in), target :: experiment
      real(c_double), intent(in) :: X(:)
      real(c_double), intent(out) :: Y(:), Cstr(:)
    End Function mlf_experience_run

    Subroutine mlf_real_get_nparameters(this, nPar, nCstr)
      import :: mlf_model_real_parameters
      class(mlf_model_real_parameters), intent(in), target :: this
      integer(8), intent(out) :: nPar, nCstr
    End Subroutine mlf_real_get_nparameters

    Integer(8) Function mlf_real_get_noutput(this)
      import :: mlf_result_vectReal
      class(mlf_result_vectReal), intent(in), target :: this
    End Function mlf_real_get_noutput

    Integer Function mlf_real_parameters_set(this, X, Cstr)
      Use iso_c_binding
      import :: mlf_model_real_parameters
      class(mlf_model_real_parameters), intent(inout), target :: this
      real(c_double), intent(in) :: X(:)
      real(c_double), intent(out), optional :: Cstr(:)
    End Function mlf_real_parameters_set

    Integer Function mlf_model_setupModel(this, param, experiment)
      import :: mlf_experience_model, mlf_model_experiment, mlf_model_parameters
      class(mlf_experience_model), intent(inout), target :: this
      class(mlf_model_parameters), intent(in) :: param
      class(mlf_model_experiment), intent(in) :: experiment
    End Function mlf_model_setupModel

    Integer Function mlf_local_init_runner(this, runner)
      import :: mlf_local_initializer, mlf_local_experience_runner
      class(mlf_local_initializer), intent(inout) :: this
      class(mlf_local_experience_runner), intent(inout) :: runner
    End Function mlf_local_init_runner
  End Interface
Contains
  Integer Function mlf_experience_dummy_results(this, results, Cstr) Result(info)
    class(mlf_experience_model), intent(inout), target :: this
    class(mlf_model_result), intent(inout) :: results
    real(c_double), intent(out), optional :: Cstr(:)
    info = -1 ! You shall reimplement this function to get results
  End Function mlf_experience_dummy_results

  Integer(8) Function mlf_experience_dummy_getNCstr(this) Result(N)
    class(mlf_experience_model), intent(in) :: this
    N = 0_8
  End Function mlf_experience_dummy_getNCstr

  Integer Function mlf_local_runExp(this, experiment, X, Y, Cstr) Result(info)
    class(mlf_local_experience_runner), intent(inout), target :: this
    class(mlf_model_experiment), intent(in), target :: experiment
    real(c_double), intent(in) :: X(:)
    real(c_double), intent(out) :: Y(:), Cstr(:)
    integer(8) :: iMax, nCstr, nIn
    iMax = HUGE(iMax) ! Stop when finished
    Cstr = 0
    CALL this%params%getNParameters(nCstr, nIn)
    info = this%params%set(X, Cstr(1:nCstr))
    If(info /= 0) Then
      Cstr(nCstr+1:) = HUGE(1d0)
      RETURN
    Endif
    info = this%model%reinit()
    If(info < 0) RETURN
    info = this%model%setupModel(this%params, experiment)
    If(info < 0) RETURN
    info = this%model%step(NITER = iMax)
    If(info < 0) RETURN
    info = this%model%getResults(this%results, Cstr(nCstr+1:))
    If(info < 0) RETURN
    info = this%results%vectResults(Y)
  End Function mlf_local_runExp
End Module mlf_supervised_model

