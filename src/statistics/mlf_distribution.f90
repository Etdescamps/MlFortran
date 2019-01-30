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

Module mlf_distribution
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_utils
  Use mlf_cfuns
  IMPLICIT NONE
  PRIVATE

  Integer, Public, Parameter :: mlf_dist_mode = 1, mlf_dist_median = 2, mlf_dist_mean = 3
  Integer, Public, Parameter :: mlf_dist_geom_mean = 4

  Type, Public, Abstract :: mlf_distribution_abstract
  End Type mlf_distribution_abstract

  Type, Public, Abstract, Extends(mlf_distribution_abstract) :: mlf_distribution_type
  Contains
    procedure(mlf_distribution_fitWithData), deferred :: fitWithData
    procedure(mlf_distribution_computePDF), deferred :: computePDF
    procedure(mlf_distribution_computePDF), deferred :: computeLogPDF
  End Type mlf_distribution_type

  Type, Public, Abstract, Extends(mlf_distribution_type) :: mlf_distribution_univariate
  Contains
    procedure(mlf_distribution_getStats), deferred :: getStats
  End Type mlf_distribution_univariate

  Type, Public, Abstract, Extends(mlf_distribution_univariate) :: mlf_distributionWithCDF_type
  Contains
    procedure(mlf_distribution_computeCDF), deferred :: computeCDF
    procedure :: integrateIncomplete => mlf_distribution_integrateIncomplete
  End Type mlf_distributionWithCDF_type

  Type, Public, Abstract, Extends(mlf_distributionWithCDF_type) :: mlf_distributionWithQuantile_type
  Contains
    procedure(mlf_distribution_quantile), deferred :: quantile
    procedure :: quantileTable => mlf_distribution_quantileTable
  End Type mlf_distributionWithQuantile_type

  Type, Public, Abstract, Extends(mlf_distributionWithQuantile_type) :: mlf_distributionWithPrior_type
  Contains
    procedure(mlf_distribution_fitWithDataWithPrior), deferred :: fitWithDataWithPrior
  End Type mlf_distributionWithPrior_type

  Abstract Interface
    Integer Function mlf_distribution_fitWithData(this, Points, W)
      Use iso_c_binding
      import :: mlf_distribution_type
      class(mlf_distribution_type), intent(inout) :: this
      real(c_double), intent(in) :: Points(:,:)
      real(c_double), intent(in), optional :: W(:)
    End Function mlf_distribution_fitWithData

    Function mlf_distribution_computePDF(this, X) Result(y)
      Use iso_c_binding
      import :: mlf_distribution_type
      class(mlf_distribution_type), intent(in) :: this
      real(c_double), intent(in) :: X(:)
      real(c_double) :: y
    End Function mlf_distribution_computePDF


    Function mlf_distribution_getStats(this, statType) Result(y)
      Use iso_c_binding
      import :: mlf_distribution_univariate
      class(mlf_distribution_univariate), intent(in) :: this
      integer, intent(in) :: statType
      real(c_double) :: y
    End Function mlf_distribution_getStats

    Integer Function mlf_distribution_fitWithDataWithPrior(this, Points, prior)
      Use iso_c_binding
      import :: mlf_distributionWithPrior_type, mlf_distribution_abstract
      class(mlf_distributionWithPrior_type), intent(inout) :: this
      real(c_double), intent(in) :: Points(:,:)
      class(mlf_distribution_abstract), intent(in) :: prior
    End Function mlf_distribution_fitWithDataWithPrior

    Function mlf_distribution_quantile(this, y) Result(x)
      Use iso_c_binding
      import :: mlf_distributionWithQuantile_type
      class(mlf_distributionWithQuantile_type), intent(in) :: this
      real(c_double), intent(in) :: y
      real(c_double) :: x
    End Function mlf_distribution_quantile

    Function mlf_distribution_computeCDF(this, x) Result(y)
      Use iso_c_binding
      import :: mlf_distributionWithCDF_type
      class(mlf_distributionWithCDF_type), intent(in) :: this
      real(c_double), intent(in) :: x
      real(c_double) :: y
    End Function mlf_distribution_computeCDF

  End Interface
Contains
  Subroutine mlf_distribution_quantileTable(this, X)
    class(mlf_distributionWithQuantile_type), intent(in) :: this
    real(c_double), intent(out) :: X(:)
    integer :: i
    Do i= 1,SIZE(X)
      X(i) = this%quantile(REAL(i-1,8)/REAL(SIZE(X)-1,8))
    End Do
  End Subroutine mlf_distribution_quantileTable

  Real(c_double) Function mlf_distribution_integrateIncomplete(this, X, z, dydz) Result(y)
    class(mlf_distributionWithCDF_type), intent(in) :: this
    real(c_double), intent(in) :: X(:), z
    real(c_double), intent(out), optional :: dydz
    real(c_double), parameter :: icoeff(3) =  [3d0/8d0, 7d0/6d0, 23d0/24d0]
    integer :: i, N
    N = SIZE(X)
    y = DOT_PRODUCT(icoeff, [IC(X(1)), IC(X(2)), IC(X(3))])
    Do i = 4, N-3
      y = y + IC(X(i))
    End Do
    y = y + DOT_PRODUCT(icoeff, [IC(X(N)), IC(X(N-1)), IC(X(N-2))])
    If(PRESENT(dydz)) Then
      dydz = DOT_PRODUCT(icoeff, [ID(X(1)), ID(X(2)), ID(X(3))])
      Do i = 4, N-3
        dydz = dydz + ID(X(i))
      End Do
      dydz = dydz + DOT_PRODUCT(icoeff, [ID(X(N)), ID(X(N-1)), ID(X(N-2))])
    Endif
  Contains
    Real(c_double) Function IC(t) Result(r)
      real(c_double), intent(in) :: t
      real(c_double) :: delta
      delta = ABS(t-z)
      If(delta == 0d0) Then
        r = 0
      Else
        r = this%computeCDF(t+delta) - this%computeCDF(t-delta)
      Endif
    End Function IC
    Real(c_double) Function ID(t) Result(r)
      real(c_double), intent(in) :: t
      real(c_double) :: delta
      delta = ABS(t-z)
      If(delta == 0d0) Then
        r = 0
      Else
        r = SIGN(this%computePDF([t+delta]) - this%computePDF([t-delta]), t-z)
      Endif
    End Function ID
  End Function mlf_distribution_integrateIncomplete

End Module mlf_distribution

