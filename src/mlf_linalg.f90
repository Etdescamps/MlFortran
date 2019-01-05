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

! This module contains a more readable version (in Fortran 2008) 
! of the linear equation solver implemented as ACM's Algorithm 423

Module mlf_linalg
  Use ieee_arithmetic
  Use iso_c_binding  
  IMPLICIT NONE
  PRIVATE

  Public :: DecompMatrix, SolveMatrix

  Interface DecompMatrix
    module procedure DecompMatrixDouble
    module procedure DecompMatrixComplex
  End Interface DecompMatrix

  Interface SolveMatrix
    module procedure SolveMatrixDouble
    module procedure SolveMatrixComplex
  End Interface SolveMatrix
Contains
   Integer Function DecompMatrixDouble(A, IP) Result(K)
    real(8), intent(inout) :: A(:,:)
    integer, intent(out) :: IP(:)
    integer :: N, M, J
    real(8) :: T
    N = SIZE(A,2)
    IP(N) = 1
    Do K = 1, N-1
      M = K - 1 + MAXLOC(ABS(A(K:,K)), DIM = 1)
      IP(K) = M
      If(M /= K) Then
        IP(N) = -IP(N)
        T = A(M,K); A(M,K) = A(K,K); A(K,K) = T ! Swap A(K,K) and A(M,K)
      Endif
      If(A(K,K) == 0d0) Then
        IP(N) = 0; RETURN
      Endif
      A(K+1:,K) = -A(K+1:,K)/A(K,K)
      Do J = K+1, N
        T = A(M,J); A(M,J) = A(K,J); A(K,J) = T ! Swap A(M,J) and A(K,J)
        If(A(K,J) /= 0d0) A(K+1:,J) = A(K+1:,J)+A(K+1:,K)*A(K,J)
      End Do
    End Do
    If(A(N,N) == 0d0) Then
      IP(N) = 0
    Else
      K = 0
    Endif
  End Function DecompMatrixDouble

  Subroutine SolveMatrixDouble(A, B, IP)
    real(8), intent(in) :: A(:,:)
    real(8), intent(inout) :: B(:)
    integer, intent(in) :: IP(:)
    integer :: K, N, M
    real(8) :: T
    N = SIZE(A,2)
    Do K = 1, N-1
      M = IP(K)
      T = B(M); B(M) = B(K); B(K) = T ! Swap B(M) and B(K)
      B(K+1:) = B(K+1:) + A(K+1:,K)*B(K)
    End Do
    Do K = N, 2, -1
      B(K) = B(K)/A(K,K)
      B(:K-1) = B(:K-1) - A(:K-1, K)*B(K)
    End Do
    B(1) = B(1)/A(1,1)
  End Subroutine SolveMatrixDouble

  Integer Function DecompMatrixComplex(A, IP) Result(K)
    complex(KIND = 8), intent(inout) :: A(:,:)
    integer, intent(out) :: IP(:)
    integer :: N, M, J
    complex(KIND = 8) :: T
    N = SIZE(A,2)
    IP(N) = 1
    Do K = 1, N-1
      M = K - 1 + MAXLOC(ABS(A(K:,K)), DIM = 1)
      IP(K) = M
      If(M /= K) Then
        IP(N) = -IP(N)
        T = A(M,K); A(M,K) = A(K,K); A(K,K) = T ! Swap A(K,K) and A(M,K)
      Endif
      If(ABS(A(K,K)) == 0d0) Then
        IP(N) = 0; RETURN
      Endif
      A(K+1:,K) = -A(K+1:,K)/A(K,K)
      Do J = K+1, N
        T = A(M,J); A(M,J) = A(K,J); A(K,J) = T ! Swap A(M,J) and A(K,J)
        If(A(K,J) /= 0d0) A(K+1:,J) = A(K+1:,J)+A(K+1:,K)*A(K,J)
      End Do
    End Do
    If(ABS(A(N,N)) == 0d0) Then
      IP(N) = 0
    Else
      K = 0
    Endif
  End Function DecompMatrixComplex

  Subroutine SolveMatrixComplex(A, B, IP)
    complex(KIND = 8), intent(in) :: A(:,:)
    complex(KIND = 8), intent(inout) :: B(:)
    integer, intent(in) :: IP(:)
    integer :: K, N, M
    complex(KIND = 8) :: T
    N = SIZE(A,2)
    Do K = 1, N-1
      M = IP(K)
      T = B(M); B(M) = B(K); B(K) = T ! Swap B(M) and B(K)
      B(K+1:) = B(K+1:) + A(K+1:,K)*B(K)
    End Do
    Do K = N, 2, -1
      B(K) = B(K)/A(K,K)
      B(:K-1) = B(:K-1) - A(:K-1, K)*B(K)
    End Do
    B(1) = B(1)/A(1,1)
  End Subroutine SolveMatrixComplex


End Module mlf_linalg

