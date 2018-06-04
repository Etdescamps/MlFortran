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

module test_common
  use ieee_arithmetic
  use iso_c_binding
  implicit none

contains
  integer function GetIntParameter(k) result(n)
    integer, intent(in) :: k
    integer :: larg, stat
    character(len=64) :: arg
    larg = 64
    call GET_COMMAND_ARGUMENT(k, arg, larg)
    read(arg,*,iostat=stat) n
    if(stat /= 0) then
      print *,"Error invalid argument ", k, " : ", arg
      stop
    endif
  end function GetIntParameter
  real(c_double) function GetRealParameter(k) result(val)
    integer, intent(in) :: k
    integer :: larg, stat
    character(len=64) :: arg
    larg = 64
    call GET_COMMAND_ARGUMENT(k, arg, larg)
    read(arg,*,iostat=stat) val
    if(stat /= 0) then
      print *,"Error invalid argument ", k, " : ", arg
      stop
    endif
  end function GetRealParameter
end module
