! Copyright (c) 2017-2019 Etienne Descamps
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

Program test_sumtree
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_sumtree
  Use mlf_binomial
  IMPLICIT NONE
  integer :: r
  r = t_sumtree(100, 4, 2, 1000000_8)
  r = t_sumtree(10000, 4, 4, 1000000000_8)
Contains
  Integer Function t_sumtree(nIndiv, nVect, nAct, nTests) Result(info)
    integer, intent(in) :: nIndiv, nVect, nAct
    integer(8), intent(in) :: nTests
    real(c_double) :: vSel(nVect), p, sumP, val, v0
    type(mlf_sumtree_double) :: T
    integer(8) :: i
    integer :: j, k
    integer, allocatable :: cnt(:,:)
    real(c_double), allocatable :: Table(:,:,:)
    info = -1
    If(T%init(nVect, nAct) /= 0) RETURN
    ALLOCATE(Table(nVect, nAct, nIndiv))
    Do j = 1, nIndiv
      CALL RANDOM_NUMBER(Table(:,:,j))
      Do k = 1,nAct
        Table(:,k,j) = k*j*Table(:,k,j)
      End Do
      info = T%updateIndiv(j, Table(:,:,j))
      If(info /= 0) RETURN
    End Do
    ALLOCATE(cnt(nAct, nIndiv))
    cnt = 0
    CALL RANDOM_NUMBER(vSel)
    Do i = 1, nTests
      k = T%getAction(vSel, val)
      j = T%getIndiv(vSel, k, val) 
      cnt(k,j) = cnt(k,j) + 1
    End Do
    sumP = T%getTop(vSel)
    val = 1d0
    Do j = 1, nIndiv
      Do k = 1, nAct
        p = DOT_PRODUCT(Table(:,k,j), vSel)/sumP
        v0 = BinomialLikelyhood(p, nTests, INT(cnt(k,j),8))
        If(val > v0) Then
          val = v0
        Endif
      End Do
    End Do
    PRINT *, "Most unlikely val: ", val
    info = 0
  End Function t_sumtree
End Program test_sumtree

