! Copyright (c) 2018 Etienne Descamps
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
  End Type mlf_result_vectReal

  Type, Public, Abstract, Extends(mlf_step_obj) :: mlf_experience_model
  Contains
    procedure(mlf_model_setExperiment), deferred :: setExperiment
    procedure(mlf_model_setParameters), deferred :: setParameters
    procedure :: getResults => mlf_experience_dummy_results
  End Type mlf_experience_model

  Type, Public :: mlf_data_point
    class(mlf_model_experiment), allocatable :: experiment
    class(mlf_model_result), allocatable :: results
    real(c_double) :: weight
  End Type mlf_data_point

  Type, Public, Abstract, Extends(mlf_obj) :: mlf_experience_runner
  Contains
    procedure(mlf_experience_run), deferred :: runExp
  End Type mlf_experience_runner

  Type, Public, Extends(mlf_experience_runner) :: mlf_local_experience_runner
    class(mlf_experience_model), pointer :: model
    class(mlf_model_real_parameters), allocatable :: params
    class(mlf_result_vectReal), allocatable :: results
  Contains
    procedure :: runExp => mlf_local_runExp
  End Type mlf_local_experience_runner

  !Type, Public, Abstract, Extends(mlf_objective_fun) :: mlf_supervised_deterministic
  !  type(mlf_data_point), allocatable :: dataSet(:)
  !  type(mlf_local_experience_runner), allocatable :: runners(:)
  !  integer(8), allocatable :: selectedIds(:)
  !Contains
  !  procedure :: eval => mlf_supervised_deterministic_eval
  !End Type mlf_supervised_deterministic

  Abstract Interface
    Integer Function mlf_vectResults(this, Y)
      Use iso_c_binding
      import :: mlf_result_vectReal
      class(mlf_result_vectReal), intent(inout), target :: this
      real(c_double), intent(out) :: Y(:)
    End Function mlf_vectResults

    Integer Function mlf_experience_run(this, experiment, X, Y)
      Use iso_c_binding
      import :: mlf_experience_runner, mlf_model_experiment
      class(mlf_experience_runner), intent(inout), target :: this
      class(mlf_model_experiment), intent(in), target :: experiment
      real(c_double), intent(in) :: X(:)
      real(c_double), intent(out) :: Y(:)
    End Function mlf_experience_run

    Integer(8) Function mlf_real_get_nparameters(this)
      import :: mlf_model_real_parameters
      class(mlf_model_real_parameters), intent(inout), target :: this
    End Function mlf_real_get_nparameters

    Integer Function mlf_real_parameters_set(this, X)
      Use iso_c_binding
      import :: mlf_model_real_parameters
      class(mlf_model_real_parameters), intent(inout), target :: this
      real(c_double), intent(in) :: X(:)
    End Function mlf_real_parameters_set

    Integer Function mlf_model_setExperiment(this, experiment)
      import :: mlf_experience_model, mlf_model_experiment
      class(mlf_experience_model), intent(inout), target :: this
      class(mlf_model_experiment), intent(in) :: experiment
    End Function mlf_model_setExperiment

    Integer Function mlf_model_setParameters(this, param)
      import :: mlf_experience_model, mlf_model_parameters
      class(mlf_experience_model), intent(inout), target :: this
      class(mlf_model_parameters), intent(in) :: param
    End Function mlf_model_setParameters
  End Interface

Contains
  Integer Function mlf_experience_dummy_results(this, results) Result(info)
    class(mlf_experience_model), intent(inout), target :: this
    class(mlf_model_result), intent(inout) :: results
    info = -1 ! You shall reimplement this function to get results
  End Function mlf_experience_dummy_results

  Integer Function mlf_local_runExp(this, experiment, X, Y) Result(info)
    class(mlf_local_experience_runner), intent(inout), target :: this
    class(mlf_model_experiment), intent(in), target :: experiment
    real(c_double), intent(in) :: X(:)
    real(c_double), intent(out) :: Y(:)
    integer(8) :: iMax
    iMax = HUGE(iMax) ! Stop when finished
    info = this%params%set(X)
    If(info < 0) RETURN
    info = this%model%reinit()
    If(info < 0) RETURN
    info = this%model%setParameters(this%params)
    If(info < 0) RETURN
    info = this%model%setExperiment(experiment)
    If(info < 0) RETURN
    info = this%model%step(NITER = iMax)
    If(info < 0) RETURN
    info = this%model%getResults(this%results)
    If(info < 0) RETURN
    info = this%results%vectResults(Y)
  End Function mlf_local_runExp
End Module mlf_supervised_model

