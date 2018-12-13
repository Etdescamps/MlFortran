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

Program test_inverse_cdf
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_inverse_cdf
  Use mlf_rand
  Use mlf_histogram
  Use mlf_hdf5
  Use mlf_intf
  type(mlf_inverseCDFSampler) :: sampler
  type(mlf_histogram_int) :: hist
  real(c_double) :: X(16), Y(16), hX(999)
  real(c_double) :: oX(1000), oY(1000), r
  integer :: info, i
  type(mlf_hdf5_file) :: h5f
  info = mlf_init()
  info = h5f%createFile("histogram.h5")
  X = [((i-1)/15d0, i=1,16)]
  Y = X*X
  info = sampler%init(X, Y)
  hX = [(i/1000d0, i=1,999)]
  info = hist%init(hX)
  If(info < 0) GOTO 1
  Do i = 1, 50000000
    CALL hist%addPoint(sampler%random(), 1_8)
    !CALL RANDOM_NUMBER(r)
    !CALL hist%addPoint(r, 1_8)
  End Do
  info = hist%get(oY, oX)
  info = h5f%pushData(oX, "oX")
  info = h5f%pushData(oY, "oY")
1 CALL h5f%finalize()
  info = mlf_quit()
End

