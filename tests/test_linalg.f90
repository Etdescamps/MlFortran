! Copyright (c) 2019 Etienne Descamps
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

Program test_linalg
  Use mlf_linalg
  CALL test_real()
  CALL test_cmplx()
Contains
  Subroutine test_cmplx()
    complex(kind = 8) :: A(3,3), B(3), C(3)
    integer :: ier, ip(3), i
    A(1,:) = CMPLX([1, 0, -1], [-1, 1, 0], 8)
    A(2,:) = CMPLX([1, 1, 0], [0, -1, 0], 8)
    A(3,:) = CMPLX([1, -1, 1], [1, 0, -1], 8)
    C = CMPLX([2, 5, -3], [-2, 3, 1])
    PRINT *, "C: ", C
    B = MATMUL(A, C)
    PRINT *, "B: ", B
    ier = DecompMatrix(A, ip)
    PRINT *, ier, ip
    Do i = 1, 4
      PRINT *, A(i,:)
    End Do
    CALL SolveMatrix(A, B, ip)
    PRINT *, "B: ", B
  End Subroutine test_cmplx

  Subroutine test_real()
    real(8) :: A(4,4), B(4), C(4)
    integer :: ier, ip(4), i
    A(1,:) = [3, 2, 3, 4]
    A(2,:) = [4, 4, 3, 2]
    A(3,:) = [1, 4, 4, 3]
    A(4,:) = [2, 3, 1, 1]
    C = [1, -1, 2, -3]
    PRINT *, "C: ", C
    B = MATMUL(A, C)
    PRINT *, "B: ", B
    ier = DecompMatrix(A, ip)
    PRINT *, ier, ip
    Do i = 1, 4
      PRINT *, A(i,:)
    End Do
    CALL SolveMatrix(A, B, ip)
    PRINT *, "B: ", B
  End Subroutine test_real
End

