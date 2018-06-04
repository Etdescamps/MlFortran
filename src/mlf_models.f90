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

Module mlf_models
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_step_algo
  IMPLICIT NONE
  PRIVATE
  
  Type, public, abstract, extends(mlf_step_obj) :: mlf_class_model
  Contains
    procedure (mlf_model_getClass), deferred :: getClass
  End Type mlf_class_model
   
  Type, public, abstract, extends(mlf_class_model) :: mlf_class_proba_model
  Contains
    procedure (mlf_model_getProba), deferred :: getProba
    procedure (mlf_model_getNumClasses), deferred :: getNumClasses
    procedure :: getClass => mlf_model_getClass_proba
  End Type mlf_class_proba_model

  Abstract Interface
    integer Function mlf_model_getClass(this, X, Cl)
      use iso_c_binding
      import :: mlf_class_model
      class(mlf_class_model), intent(in), target :: this
      real(c_double), intent(in) :: X(:,:)
      integer(c_int), intent(out) :: Cl(:)
    End Function mlf_model_getClass

    integer Function mlf_model_getNumClasses(this)
      use iso_c_binding
      import :: mlf_class_proba_model
      class(mlf_class_proba_model), intent(in), target :: this
    End Function mlf_model_getNumClasses

    integer Function mlf_model_getProba(this, X, Proba)
      use iso_c_binding
      import :: mlf_class_proba_model
      class(mlf_class_proba_model), intent(in), target :: this
      real(c_double), intent(in) :: X(:,:)
      real(c_double), intent(out) :: Proba(:,:)
    End Function mlf_model_getProba
  End Interface
Contains
  integer Function mlf_model_getClass_proba(this, X, Cl) result(info)
    class(mlf_class_proba_model), intent(in), target :: this
    real(c_double), intent(in) :: X(:,:)
    integer(c_int), intent(out) :: Cl(:)
    real(c_double), allocatable :: Proba(:,:)
    integer :: nX, nC
    nX = size(X,2); nC = this%getNumClasses()
    allocate(Proba(nX, nC))
    info = this%getProba(X, Proba)
    cl = maxloc(Proba, dim=2)
  End Function mlf_model_getClass_proba

End Module mlf_models

