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

Module mlf_normal_dist
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_utils
  Use mlf_cfuns
  Use mlf_distribution
  IMPLICIT NONE
  PRIVATE

  Type, Public, Extends(mlf_distributionWithQuantile_type) :: mlf_normal_distribution
    real(c_double) :: mu, sigma
  Contains
    procedure :: fitWithData => Normal_fitWithData
    procedure :: quantile => Normal_quantile
    procedure :: computeCDF => Normal_computeCDF
    procedure :: computePDF => Normal_computePDF
  End Type mlf_normal_distribution

  Type, Public, Extends(mlf_normal_distribution) :: mlf_logNormal_distribution
  Contains
    procedure :: fitWithData => logNormal_fitWithData
    procedure :: quantile => logNormal_quantile
    procedure :: computeCDF => logNormal_computeCDF
    procedure :: computePDF => logNormal_computePDF
  End Type mlf_logNormal_distribution
Contains
  Integer Function Normal_fitWithData(this, Points, W) Result(info)
    class(mlf_normal_distribution), intent(inout) :: this
    real(c_double), intent(in) :: Points(:,:)
    real(c_double), intent(in), optional :: W(:)
    If(PRESENT(W)) Then
      this%mu = SUM(W*Points(1,:))/SUM(W)
      this%sigma = SUM(W*(Points(1,:)-this%mu)**2)/(1/SUM(W)-SUM(W*W)*SUM(W))
    Else
      this%mu = SUM(Points(1,:))/SIZE(Points,2)
      this%sigma = SUM((Points(1,:)-this%mu)**2/(SIZE(Points,2)-1d0))
    Endif
    info = 0
  End Function Normal_fitWithData

  Real(c_double) Function Normal_computeCDF(this, x) Result(y)
    class(mlf_normal_distribution), intent(in) :: this
    real(c_double), intent(in) :: x
    y = 0.5d0*(1+c_erf((x-this%mu)/(this%sigma*SQRT(2d0))))
  End Function Normal_computeCDF

  Real(c_double) Function Normal_computePDF(this, x) Result(y)
    class(mlf_normal_distribution), intent(in) :: this
    real(c_double), intent(in) :: x
    y = 1d0/SQRT(2d0*mlf_PI*this%sigma)*EXP(-(x-this%mu)**2/(2*this%sigma**2))
  End Function Normal_computePDF

  Real(c_double) Function Normal_quantile(this, y) Result(x)
    class(mlf_normal_distribution), intent(in) :: this
    real(c_double), intent(in) :: y
    x = this%mu - this%sigma*SQRT(2d0)*c_erf_inv(2*(1-y))
  End Function Normal_quantile

  Real(c_double) Function logNormal_quantile(this, y) Result(x)
    class(mlf_logNormal_distribution), intent(in) :: this
    real(c_double), intent(in) :: y
    x = EXP(this%mu - this%sigma*SQRT(2d0)*c_erf_inv(2*(1-y)))
  End Function logNormal_quantile

  Integer Function logNormal_fitWithData(this, Points, W) Result(info)
    class(mlf_logNormal_distribution), intent(inout) :: this
    real(c_double), intent(in) :: Points(:,:)
    real(c_double), intent(in), optional :: W(:)
    info = Normal_fitWithData(this, LOG(Points), W)
  End Function logNormal_fitWithData

  Real(c_double) Function logNormal_computeCDF(this, x) Result(y)
    class(mlf_logNormal_distribution), intent(in) :: this
    real(c_double), intent(in) :: x
    y = 0.5d0*(1+c_erf((LOG(x)-this%mu)/(this%sigma*SQRT(2d0))))
  End Function logNormal_computeCDF

  Real(c_double) Function logNormal_computePDF(this, x) Result(y)
    class(mlf_logNormal_distribution), intent(in) :: this
    real(c_double), intent(in) :: x
    y = 1d0/(x*this%sigma*SQRT(2d0*mlf_PI))*EXP(-(LOG(x)-this%mu)**2/(2*this%sigma**2))
  End Function logNormal_computePDF
End Module mlf_normal_dist

