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

  Type, Public, Abstract :: mlf_distribution_abstract
  End Type mlf_distribution_abstract

  Type, Public, Abstract, Extends(mlf_distribution_abstract) :: mlf_distribution_type
  Contains
    procedure(mlf_distribution_fitWithData), deferred :: fitWithData
    procedure(mlf_distribution_inverseInterval), deferred :: inverseInterval
    procedure(mlf_distribution_integrateIncomplete), deferred :: integrateIncomplete
  End Type mlf_distribution_type

  Type, Public, Abstract, Extends(mlf_distribution_type) :: mlf_distributionWithPrior_type
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

    Integer Function mlf_distribution_fitWithDataWithPrior(this, Points, prior)
      Use iso_c_binding
      import :: mlf_distributionWithPrior_type, mlf_distribution_abstract
      class(mlf_distributionWithPrior_type), intent(inout) :: this
      real(c_double), intent(in) :: Points(:,:)
      class(mlf_distribution_abstract), intent(in) :: prior
    End Function mlf_distribution_fitWithDataWithPrior

    Subroutine mlf_distribution_inverseInterval(this, X)
      Use iso_c_binding
      import :: mlf_distribution_type
      class(mlf_distribution_type), intent(in) :: this
      real(c_double), intent(out) :: X(:)
    End Subroutine mlf_distribution_inverseInterval

    Function mlf_distribution_integrateIncomplete(this, X, z) Result(y)
      Use iso_c_binding
      import :: mlf_distribution_type
      class(mlf_distribution_type), intent(in) :: this
      real(c_double), intent(in) :: X(:), z
      real(c_double) :: y
    End Function mlf_distribution_integrateIncomplete
  End Interface
End Module mlf_distribution

