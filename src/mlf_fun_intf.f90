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

Module mlf_fun_intf
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  IMPLICIT NONE
  PRIVATE
 
  Public :: c_objfunction, c_basisfunction
  Type, bind(C) :: mlf_objfuninfo
    integer(c_int) :: nDimIn, nDimCstr, nDimOut
  End Type mlf_objfuninfo

  ! Function handler ((multi/mono) objective fitness function)
  Type, Public, abstract, extends(mlf_obj) :: mlf_objective_fun
    integer :: nD, nC, nY
  Contains
    procedure (mlf_obj_eval), deferred :: eval
    procedure (mlf_obj_eval), deferred :: constraints
  End Type mlf_objective_fun

  ! C function wrapper
  Type, Public, extends(mlf_objective_fun) :: mlf_objective_fun_c
    procedure (mlf_obj_eval_c), nopass, pointer :: evalC => NULL()
    procedure (mlf_obj_eval_c), nopass, pointer :: constraintsC => NULL()
    type(c_ptr) :: ptr
  Contains
    procedure :: eval => mlf_obj_c_eval, constraints => mlf_obj_c_constraints
  End Type mlf_objective_fun_c

  ! Function handler for basis functions
  Type, Public, abstract, extends(mlf_obj) :: mlf_basis_fun
  Contains
    procedure (mlf_basis_eval), deferred :: eval
  End Type mlf_basis_fun

  ! C function wrapper
  Type, Public, extends(mlf_basis_fun) :: mlf_basis_fun_c
    procedure (mlf_basis_eval_c), nopass, pointer :: evalC => NULL()
    type(c_ptr) :: ptr
  Contains
    procedure :: eval => mlf_basis_c_eval
  End Type mlf_basis_fun_c

  ! Function handler for weight functions
  Type, Public, abstract, extends(mlf_obj) :: mlf_weight_fun
  Contains
    procedure (mlf_weight_eval), deferred :: eval
  End Type mlf_weight_fun

  ! Simple weight functions handler
  Type, Public, abstract, extends(mlf_weight_fun) :: mlf_weight_simple_fun
    procedure (mlf_weight_simple_eval), nopass, pointer :: evalS
  Contains
    procedure :: eval => mlf_weight_simple_ev
  End Type mlf_weight_simple_fun

  ! ODE function evaluator
  Type, Public, abstract, extends(mlf_obj) :: mlf_ode_fun
    integer :: idConst = -1
  Contains
    procedure (mlf_ode_eval), deferred :: eval
  End Type mlf_ode_fun

  Abstract Interface

    ! Abstract ODE function type
    Subroutine mlf_ode_eval(this, t, X, F)
      Use iso_c_binding
      import :: mlf_ode_fun
      class(mlf_ode_fun), intent(in), target :: this
      real(c_double), intent(in) :: t
      real(c_double), intent(in), target :: X(:)
      real(c_double), intent(out), target :: F(:)
    End Subroutine mlf_ode_eval

    ! Abstract objective function type
    integer Function mlf_obj_eval(this, X, Y)
      Use iso_c_binding
      import :: mlf_objective_fun
      class(mlf_objective_fun), intent(in), target :: this
      real(c_double), intent(in), target :: X(:,:)
      real(c_double), intent(inout), target:: Y(:,:)
    End Function mlf_obj_eval

    ! C interface for objective functions
    Function mlf_obj_eval_c(X, Y, ND, NY, lambda, ptr) bind(C)
      Use iso_c_binding
      integer(c_int), intent(in), value :: ND, NY, lambda
      integer(c_int) :: mlf_obj_eval_c
      type(c_ptr), value :: ptr, X, Y
    End Function mlf_obj_eval_c

    ! Abstract basis function type (for dimension reduction)
    integer Function mlf_basis_eval(this, X, rpar, Y)
      use iso_c_binding
      import :: mlf_basis_fun
      class(mlf_basis_fun), intent(in), target :: this
      real(c_double), intent(in), target :: X(:), rpar(:,:)
      real(c_double), intent(out), target :: Y(:,:)
    End Function mlf_basis_eval

    ! C interface for basis functions
    Function mlf_basis_eval_c(X, rpar, Y, nX, nPar, sPar, ptr) bind(C)
      use iso_c_binding
      integer(c_int), intent(in), value :: nX, nPar, sPar
      type(c_ptr), value :: ptr, X, Y, rpar
      integer(c_int) :: mlf_basis_eval_c
    End Function mlf_basis_eval_c

    ! Abstract weight function type
    Subroutine mlf_weight_eval(this, Y, lambda)
      use iso_c_binding
      import :: mlf_weight_fun
      class(mlf_weight_fun), intent(in), target :: this
      real(c_double), intent(out) :: Y(:)
      integer :: lambda
    End Subroutine mlf_weight_eval

    ! Simple weight function type
    Subroutine mlf_weight_simple_eval(Y, lambda)
      use iso_c_binding
      real(c_double), intent(out) :: Y(:)
      integer :: lambda
    End Subroutine mlf_weight_simple_eval
  End Interface
Contains
  integer Function mlf_basis_c_eval(this, X, rpar, Y) result(info)
    class(mlf_basis_fun_c), intent(in), target :: this
    real(c_double), intent(in), target :: X(:), rpar(:,:)
    real(c_double), intent(out), target :: Y(:,:)
    integer(c_int) :: nX, sPar, nPar
    if(.NOT. associated(this%evalC)) then
      info = mlf_UNINIT
      RETURN
    endif
    nX = size(X,1)
    sPar = size(rpar,1)
    nPar = size(rpar,2)
    info = this%evalC(c_loc(X), c_loc(rpar), c_loc(Y), nX, nPar, sPar, this%ptr)
  End Function mlf_basis_c_eval

  integer Function mlf_obj_c_eval(this, X, Y) result(info)
    class(mlf_objective_fun_c), intent(in), target :: this
    real(c_double), intent(in), target :: X(:,:)
    real(c_double), intent(inout), target :: Y(:,:)
    integer(c_int) :: ND, NY, lambda
    if(.NOT. associated(this%evalC)) then
      info = mlf_UNINIT
      RETURN
    endif
    lambda = size(X,2)
    ND = size(X,1)
    NY = size(Y,1)
    info = this%evalC(c_loc(X), c_loc(Y), ND, nY, lambda, this%ptr)
  End Function mlf_obj_c_eval

  integer Function mlf_obj_c_constraints(this, X, Y) result(info)
    class(mlf_objective_fun_c), intent(in), target :: this
    real(c_double), intent(in), target :: X(:,:)
    real(c_double), intent(inout), target :: Y(:,:)
    integer(c_int) :: ND, NY, lambda
    if(associated(this%constraintsC)) then
      lambda = size(X,2)
      ND = size(X,1)
      NY = size(Y,1)
      info = this%constraintsC(c_loc(X), c_loc(Y), ND, nY, lambda, this%ptr)
    else
      info = mlf_OK
      Y=0
    endif
  End Function mlf_obj_c_constraints

  type(c_ptr) Function c_objfunction(cfun, cptr, ccst, funinfo) bind(C, name="mlf_objfunction")
    type(c_ptr), value :: cptr
    type(mlf_objfuninfo) :: funinfo
    type(c_funptr), value :: cfun, ccst
    type(mlf_objective_fun_c), pointer :: x
    class (*), pointer :: obj
    ALLOCATE(x)
    x%ptr = cptr
    x%nD = funinfo%nDimIn; x%nY = funinfo%nDimOut; x%nC = funinfo%nDimCstr
    if(C_ASSOCIATED(cfun)) call C_F_PROCPOINTER(cfun, x%evalC)
    if(C_ASSOCIATED(ccst)) call C_F_PROCPOINTER(ccst, x%constraintsC)
    obj => x
    c_objfunction = c_allocate(obj)
  End Function c_objfunction

  type(c_ptr) Function c_basisfunction(cfun, cptr) bind(C, name="mlf_basisfunction")
    type(c_ptr), value :: cptr
    type(c_funptr), value :: cfun
    type(mlf_basis_fun_c), pointer :: x
    class (*), pointer :: obj
    ALLOCATE(x)
    x%ptr = cptr
    if(C_ASSOCIATED(cfun)) call C_F_PROCPOINTER(cfun, x%evalC)
    obj => x
    c_basisfunction = c_allocate(obj)
  End Function c_basisfunction

  ! Abstract weight function type
  Subroutine mlf_weight_simple_ev(this, Y, lambda)
    class(mlf_weight_simple_fun), intent(in), target :: this
    real(c_double), intent(out) :: Y(:)
    integer :: lambda
    call this%evalS(Y, lambda)
  End Subroutine mlf_weight_simple_ev

End Module mlf_fun_intf
