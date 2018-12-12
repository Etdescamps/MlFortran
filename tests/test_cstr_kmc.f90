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

Program test_cstr_kmc
  Use ieee_arithmetic
  Use iso_c_binding
  Use test_cstr_kmc_model
  Use mlf_hybrid_kmc
  Use mlf_hdf5
  Use mlf_ode45
  Use mlf_step_algo
  Use mlf_fun_intf
  Use mlf_intf
  integer :: info
  type(mlf_hdf5_file) :: h5f
  real(c_double) :: P(5)
  type(test_cstr_kmc_experiment) :: modelExp
  P = [1d0, 10d0, 5d0, 1d0, 1d0]
  modelExp = test_cstr_kmc_experiment(0d0, 0d0, 10d0, 0d0, 1d0)
  info = mlf_init()
  info = h5f%createFile("kmc_cstr.h5")
  info = eval(h5f, "r1000", P, 1000d0, modelExp)
  info = eval(h5f, "r5000", P, 5000d0, modelExp)
  info = eval(h5f, "r20000", P, 20000d0, modelExp)
  info = eval(h5f, "r200000", P, 200000d0, modelExp)
  CALL h5f%finalize()
  info = mlf_quit()
Contains
  Integer Function eval(h5f, rname, Param, Volume, modelExp) Result(info)
    class(mlf_hdf5_file), intent(inout) :: h5f
    character(len=*), intent(in) :: rname
    real(c_double), intent(in) :: Param(5), Volume
    type(test_cstr_kmc_experiment), intent(in) :: modelExp
    integer(c_int64_t) :: NIndiv(2)
    type(model_cstr_kmc) :: model
    type(test_cstr_kmc_parameters) :: p
    integer(8) :: i, N
    real(c_double), allocatable :: points(:,:)
    p%Volume = Volume
    info = p%set(Param)
    info = model%init()
    info = model%setParameters(p)
    info = model%setExperiment(modelExp)
    N = INT(10*Volume, KIND = 8)
    ALLOCATE(points(7, N))
    Do i = 1,N
      info = model%step()
      points(1,i)   = model%ode%t
      points(2:3,i) = REAL(model%NIndiv, KIND=8)/Volume
      points(4:7,i) = model%ode%X(1:4)
      If(info < 0 .OR. info == mlf_ODE_StopTime) EXIT
    End Do
    i = MIN(i, N)
    info = h5f%pushData(points(:,1:i), rname)
    CALL model%finalize()
  End Function eval
End

