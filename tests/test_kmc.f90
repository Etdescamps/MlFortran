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

Program test_kmc
  Use ieee_arithmetic
  Use iso_c_binding
  Use test_kmc_model
  Use mlf_hdf5
  Use mlf_ode45
  Use mlf_step_algo
  Use mlf_fun_intf
  Use mlf_intf
  integer :: info
  type(mlf_hdf5_file) :: h5f
  info = mlf_init()
  info = h5f%createFile("kmc.h5")
  info = eval(h5f, "r200000", 10.0d0, 1d0, 200000d0, [0.4d0, 0.001d0, 0d0])
  info = eval(h5f, "r20000", 10.0d0, 1d0, 20000d0, [0.4d0, 0.001d0, 0d0])
  info = eval(h5f, "r5000", 10.0d0, 1d0, 5000d0, [0.4d0, 0.001d0, 0d0])
  info = eval(h5f, "r1000", 10.0d0, 1d0, 1000d0, [0.4d0, 0.001d0, 0d0])
  info = eval_ode(h5f, "rODE", 10.0d0, 1d0, [0.4d0, 0.001d0, 0d0])
  CALL h5f%finalize()
  info = mlf_quit()
Contains
  Integer Function eval_ode(h5f, rname, Alpha, Beta, CIndiv) Result(info)
    class(mlf_hdf5_file), intent(inout) :: h5f
    character(len=*), intent(in) :: rname
    real(c_double), intent(in) :: Alpha, Beta
    real(c_double), intent(in) :: CIndiv(3)
    real(c_double), allocatable :: points(:,:)
    type(kmc_ode) :: fun
    type(mlf_ode45_obj) :: ode
    integer, parameter :: NStepMax = 100000
    integer :: i
    fun%Alpha = Alpha; fun%Beta = Beta
    info = ode%init(fun, CIndiv, atoli = 1d-6, rtoli = 1d-6)
    ALLOCATE(points(4, NStepMax))
    Do i = 1,NStepMax
      info = ode%step()
      points(1,i)   = ode%t
      points(2:4,i) = ode%X
      If(info /= 0) EXIT
    End Do
    i = MIN(i, NStepMax)
    info = h5f%pushData(points(:,1:i), rname)
  End Function eval_ode

  Integer Function eval(h5f, rname, Alpha, Beta, Volume, CIndiv) Result(info)
    class(mlf_hdf5_file), intent(inout) :: h5f
    character(len=*), intent(in) :: rname
    real(c_double), intent(in) :: Alpha, Beta, Volume
    real(c_double), intent(in) :: CIndiv(3)
    integer(c_int64_t) :: NIndiv(3)
    type(kmc_model) :: model
    integer(8) :: i, N
    real(c_double), allocatable :: points(:,:)
    NIndiv = INT(CIndiv*Volume, KIND=8)
    info = model%init(Alpha, Beta, Volume, NIndiv)
    N = SUM(NIndiv)*2
    ALLOCATE(points(6, N))
    Do i = 1,N
      info = model%step()
      points(1,i)   = model%t
      points(2:4,i) = REAL(model%NIndiv, KIND=8)
      points(5:6,i) = model%Rates
      If(info /= 0) EXIT
    End Do
    i = MIN(i, N)
    info = h5f%pushData(points(:,1:i), rname)
    CALL model%finalize()
  End Function eval
End

