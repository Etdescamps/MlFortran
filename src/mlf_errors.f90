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

Module mlf_errors
  Use iso_fortran_env
  Use iso_c_binding
  IMPLICIT NONE
  PRIVATE

  Public :: CheckF, CheckNZ
  
  Interface CheckF
    module procedure mlf_checkF
    module procedure mlf_checkF_dims
  End Interface CheckF
  
  Interface CheckNZ
    module procedure mlf_checkNZ
    module procedure mlf_checkNZ_dims
  End Interface CheckNZ
Contains

  Logical Function mlf_checkNZ(info, msg) result(x)
    integer, intent(inout) :: info
    character(len=*), intent(in) :: msg
    x = .FALSE.
    if(info == 0) RETURN
    info = -1
    write (error_unit, *) msg
    x = .TRUE.
  End Function mlf_checkNZ

  Logical Function mlf_checkF(info, msg) result(x)
    integer, intent(in) :: info
    character(len=*), intent(in) :: msg
    x = .FALSE.
    if(info >= 0) RETURN
    write (error_unit, *) msg
    x = .TRUE.
  End Function mlf_checkF

  Logical Function mlf_checkF_dims(info, msg, dims) result(x)
    integer, intent(in) :: info
    character(len=*), intent(in) :: msg
    integer(kind=8), intent(in) :: dims(:)
    x = .FALSE.
    if(info >= 0) RETURN
    write (error_unit, *) msg, dims
    x = .TRUE.
  End Function mlf_checkF_dims

  Logical Function mlf_checkNZ_dims(info, msg, dims) result(x)
    integer, intent(inout) :: info
    character(len=*), intent(in) :: msg
    integer(kind=8), intent(in) :: dims(:)
    x = .FALSE.
    if(info == 0) RETURN
    info = -1
    write (error_unit, *) msg, dims
    x = .TRUE.
  End Function mlf_checkNZ_dims
End Module mlf_errors

