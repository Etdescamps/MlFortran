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


! Module relative to gaussian models
Module mlf_gaussian
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_utils
  Use mlf_matrix
  IMPLICIT NONE
  PRIVATE
  Public :: mlf_MaxGaussian, mlf_MaxGaussian_c, mlf_EvalGaussian, mlf_EvalGaussian_C
Contains
  ! Infert a multivariate normal distribution using array of weighted vector
  real(c_double) function mlf_MaxGaussian(X, Proba, mu, C) result(sumPi)
    real(c_double), intent(in) :: X(:,:), Proba(:) ! Vector and associated weight
    real(c_double), intent(out) :: mu(:), C(:,:)
    real(c_double), allocatable :: W(:)
    integer :: NX, ND
    NX = size(X,2); ND = size(X,1)
    sumPi = sum(Proba)
    ALLOCATE(W(size(Proba)))
    W(:) = Proba(:)/sumPi
    call GECovMatrix(C, X, Mu, W, 1d-3/real(NX))
  end function mlf_MaxGaussian
  ! C function for using MaxGaussian
  real(c_double) function mlf_MaxGaussian_c(ND, NX, X, Proba, mu, C) result(sumPi)&
      bind(C, name="mlf_maxGaussian")
    integer(c_int), value :: ND, NX
    real(c_double), intent(in) :: X(ND,NX), Proba(NX)
    real(c_double), intent(out) :: mu(ND), C(ND,ND)
    sumPi = mlf_MaxGaussian(X, Proba, mu, C)
  end function mlf_MaxGaussian_c

  ! Evaluate PDF of vectors given a peculiar normal law Norm(mu, C)
  ! with a conditional probability lambda
  integer function mlf_EvalGaussian(X, P, mu, C, lambda, sumLL) result(info)
    real(c_double), intent(in) :: X(:,:), mu(:), C(:,:), lambda
    real(c_double), intent(out) :: P(:), sumLL
    real(c_double), allocatable :: invC(:,:)
    real(c_double) :: valProd, lnDet
    integer :: ND, NX, i, N
    ND = size(X,1); NX = size(X,2); N = size(C,1)
    ALLOCATE(invC(N,N))
    info = InverseSymMatrix(C, invC, lnDet, .TRUE.)
    if(info /= 0) RETURN
    sumLL = NX*(log(lambda)-lnDet)
    do i=1,NX
      valProd = -0.5d0*dot_product(X(:,i)-mu,matmul(invC,X(:,i)-mu))
      sumLL = sumLL + valProd
      P(i) = lambda*exp(valProd-lnDet)
    end do
  end function mlf_EvalGaussian
  ! C function for using EvalGaussian
  integer(c_int) function mlf_EvalGaussian_C(ND, NX, X, P, mu, C, lambda, sumLL) &
      result(info) bind(C, name="mlf_evalGaussian")
    integer(c_int), value :: ND, NX
    real(c_double), intent(in) :: X(ND,NX), mu(ND), C(ND,ND)
    real(c_double), value :: lambda
    real(c_double), intent(out) :: P(NX), sumLL
    info = mlf_EvalGaussian(X, P, mu, C, lambda, sumLL)
  end function mlf_EvalGaussian_C
 
End Module mlf_gaussian
