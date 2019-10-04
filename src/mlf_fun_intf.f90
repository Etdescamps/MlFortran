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

Module mlf_fun_intf
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_errors
  Use mlf_rand
  IMPLICIT NONE
  PRIVATE
 
  Public :: c_objfunction, c_basisfunction
  Type, bind(C) :: mlf_objfuninfo
    integer(c_int) :: nDimIn
    integer(c_int) :: nDimOut
    integer(c_int) :: nDimCstr
  End Type mlf_objfuninfo

  ! Function handler ((multi/mono) objective fitness function)
  Type, Public, Abstract, Extends(mlf_obj) :: mlf_objective_fun
    ! Input nD parameter and output nC constraints and nY objectives
    ! Optimisation methods (mono-objectives) tends to use the lexicographic order
    ! So the first elements will be prioretary to the other elements in the objective vector
    integer :: nD, nY, nC
  Contains
    procedure (mlf_obj_eval), deferred :: eval
  End Type mlf_objective_fun

  ! C function wrapper
  Type, Public, Extends(mlf_objective_fun) :: mlf_objective_fun_c
    procedure (mlf_obj_eval_c), nopass, pointer :: evalC => NULL()
    type(c_ptr) :: ptr
  Contains
    procedure :: eval => mlf_obj_c_eval
  End Type mlf_objective_fun_c

  ! Function handler for basis functions
  Type, Public, Abstract, Extends(mlf_obj) :: mlf_basis_fun
  Contains
    procedure (mlf_basis_eval), deferred :: eval
  End Type mlf_basis_fun

  Type, Public, Abstract, Extends(mlf_basis_fun) :: mlf_basis_fun_inv
  Contains
    procedure (mlf_basis_integral), deferred :: integral
    procedure (mlf_basis_inv_integ), deferred :: inv_integ
  End Type mlf_basis_fun_inv

  ! Function handler for basis functions
  Type, Public, Abstract, Extends(mlf_obj) :: mlf_basis_funder
  Contains
    procedure (mlf_basis_evalder), deferred :: eval
  End Type mlf_basis_funder

  ! C function wrapper
  Type, Public, Extends(mlf_basis_fun) :: mlf_basis_fun_c
    procedure (mlf_basis_eval_c), nopass, pointer :: evalC => NULL()
    type(c_ptr) :: ptr
  Contains
    procedure :: eval => mlf_basis_c_eval
  End Type mlf_basis_fun_c

  ! Function handler for weight functions
  Type, Public, Abstract, Extends(mlf_obj) :: mlf_weight_fun
  Contains
    procedure (mlf_weight_eval), deferred :: eval
  End Type mlf_weight_fun

  ! Simple weight functions handler
  Type, Public, Abstract, Extends(mlf_weight_fun) :: mlf_weight_simple_fun
    procedure (mlf_weight_simple_eval), nopass, pointer :: evalS
  Contains
    procedure :: eval => mlf_weight_simple_ev
  End Type mlf_weight_simple_fun

  ! Bijection transformation object
  Type, Public, Abstract, Extends(mlf_obj) :: mlf_bijection_fun
  Contains
    procedure(mlf_bijection_proj), deferred :: proj
    procedure(mlf_bijection_invert), deferred :: invert
  End Type mlf_bijection_fun

  ! ODE class has two type of constraint:
  !   - a "hard" constraint that prevent the evaluation of derivatives:
  !     in this case the derivation function cannot be evaluated
  !     -> the function eval return a positive value instead of 0 in this case
  !     -> a negative return value will imply the stop of the algorithm
  !       -> handled by eval function of mlf_ode_fun
  !   - a "soft" constraint that permit the evaluation of the derivative
  !     -> in this case, it uses dense evaluation to find where the limit is reached

  ! ODE function evaluator
  Type, Public, abstract, extends(mlf_obj) :: mlf_ode_fun
  Contains
    procedure (mlf_ode_eval), deferred :: eval
  End Type mlf_ode_fun

  Type, Public, abstract, extends(mlf_ode_fun) :: mlf_ode_funCstr
    integer :: NCstr
  Contains
    procedure (mlf_ode_funCstr_update), deferred :: updateCstr
    procedure (mlf_ode_funCstr_getHMax), deferred :: getHMax
    procedure (mlf_ode_funCstr_reach), deferred :: reachCstr
    procedure :: checkCstr => mlf_ode_checkCstr
    procedure (mlf_ode_funCstr_getDerivatives), deferred :: getDerivatives
  End Type mlf_ode_funCstr

  Type, Public, abstract, extends(mlf_ode_funCstr) :: mlf_ode_funStop
  Contains
    procedure :: funStop => mlf_ode_dummy_funStop
  End Type mlf_ode_funStop

  Type, Public, abstract, extends(mlf_ode_funStop) :: mlf_ode_funCstrIds
    real(c_double), allocatable :: cstrValRef(:)
    real(c_double), allocatable :: cstrAlpha(:)
    integer, allocatable :: cstrSelIds(:)
    real(c_double) :: cstrT
    integer :: cstrId
  Contains
    procedure :: allocateCstr => mlf_ode_allocateCstrIds
    procedure :: getHMax => mlf_ode_getHMaxIds
    procedure :: updateCstr => mlf_ode_updateCstrIds
    procedure :: reachCstr => mlf_ode_reachCstrIds
    procedure :: getDerivatives => mlf_ode_getDerivativesIds
  End Type mlf_ode_funCstrIds


  Type, Public, abstract, extends(mlf_ode_funStop) :: mlf_ode_funCstrVect
    real(c_double), allocatable :: cstrVect(:,:)
    real(c_double), allocatable :: cstrValRef(:)
    real(c_double), allocatable :: cstrLastVal(:)
    real(c_double), allocatable :: cstrLastDer(:)
    real(c_double), allocatable :: cstrTmp(:,:)
    real(c_double), allocatable :: cstrAlpha(:)
    real(c_double) :: cstrT, epsilonT
    integer :: cstrId
  Contains
    procedure :: mlf_ode_allocateCstr, mlf_ode_allocateCstrVect
    generic :: allocateCstr => mlf_ode_allocateCstr, mlf_ode_allocateCstrVect
    procedure :: getHMax => mlf_ode_getHMax
    procedure :: updateCstr => mlf_ode_updateCstr
    procedure :: reachCstr => mlf_ode_reachCstr
    procedure :: getDerivatives => mlf_ode_getDerivatives
  End Type mlf_ode_funCstrVect

  Integer, Parameter, Public :: mlf_ODE_Continue = 0, mlf_ODE_SoftCstr = 1
  Integer, Parameter, Public :: mlf_ODE_StopTime = 2, mlf_ODE_HardCstr = 3
  Integer, Parameter, Public :: mlf_ODE_ReevaluateDer = 4, mlf_ODE_StopODE = 5
  Integer, Parameter, Public :: mlf_ODE_FunError = -1

  Abstract Interface
    Integer Function mlf_ode_funCstr_reach(this, t, tMin, tMax, ids, X, F)
      Use iso_c_binding
      import :: mlf_ode_funCstr
      class(mlf_ode_funCstr), intent(inout), target :: this
      real(c_double), intent(inout) :: t
      real(c_double), intent(in) :: tMin, tMax
      real(c_double), intent(inout), target :: X(:), F(:)
      integer, intent(in) :: ids(:)
    End Function mlf_ode_funCstr_reach

    Subroutine mlf_ode_funCstr_getDerivatives(this, ids, X0, X, K, C0, C, Q)
      Use iso_c_binding
      import :: mlf_ode_funCstr
      class(mlf_ode_funCstr), intent(inout), target :: this
      real(c_double), intent(in), target :: K(:,:), X0(:), X(:)
      real(c_double), intent(out), target :: C0(:), C(:), Q(:,:)
      integer, intent(in), target :: ids(:)
    End Subroutine mlf_ode_funCstr_getDerivatives

    Integer Function mlf_ode_funCstr_update(this, t, X0, X, F0, F, ids, hMax) Result(N)
      Use iso_c_binding
      import :: mlf_ode_funCstr
      class(mlf_ode_funCstr), intent(inout), target :: this
      real(c_double), intent(in) :: t
      real(c_double), intent(in), target :: X0(:), X(:), F0(:), F(:)
      integer, intent(out), target :: ids(:)
      real(c_double), intent(inout) :: hMax
    End Function mlf_ode_funCstr_update

    Function mlf_ode_funCstr_getHMax(this, t, X, F)
      Use iso_c_binding
      import :: mlf_ode_funCstr
      class(mlf_ode_funCstr), intent(inout), target :: this
      real(c_double), intent(in) :: t
      real(c_double), intent(in), target :: X(:), F(:)
      real(c_double) :: mlf_ode_funCstr_getHMax
    End Function mlf_ode_funCstr_getHMax

    ! Abstract ODE function type
    Integer Function mlf_ode_eval(this, t, X, F)
      Use iso_c_binding
      import :: mlf_ode_fun
      class(mlf_ode_fun), intent(inout), target :: this
      real(c_double), intent(in) :: t
      real(c_double), intent(in), target :: X(:)
      real(c_double), intent(out), target :: F(:)
    End Function mlf_ode_eval

    ! Abstract objective function type
    Integer Function mlf_obj_eval(this, X, Y)
      Use iso_c_binding
      import :: mlf_objective_fun
      class(mlf_objective_fun), intent(inout), target :: this
      real(c_double), intent(in), target :: X(:,:)
      real(c_double), intent(inout), target :: Y(:,:)
    End Function mlf_obj_eval

    ! C interface for objective functions
    Function mlf_obj_eval_c(X, Y, ND, NY, lambda, ptr) bind(C)
      Use iso_c_binding
      integer(c_int), intent(in), value :: ND, NY, lambda
      real(c_double) :: mlf_obj_eval_c
      type(c_ptr), value :: ptr, X, Y
    End Function mlf_obj_eval_c

    ! Abstract basis function type (for dimension reduction)
    Integer Function mlf_basis_eval(this, X, rpar, Y)
      use iso_c_binding
      import :: mlf_basis_fun
      class(mlf_basis_fun), intent(in), target :: this
      real(c_double), intent(in), target :: X(:), rpar(:,:)
      real(c_double), intent(out), target :: Y(:,:)
    End Function mlf_basis_eval

    ! Abstract basis function type (for dimension reduction)
    Integer Function mlf_basis_integral(this, x0, xEnd, rpar, Z)
      use iso_c_binding
      import :: mlf_basis_fun_inv
      class(mlf_basis_fun_inv), intent(in), target :: this
      real(c_double), intent(in) :: x0, xEnd
      real(c_double), intent(in), target :: rpar(:,:)
      real(c_double), intent(out), target :: Z(:)
    End Function mlf_basis_integral

    ! Abstract basis function type (for dimension reduction)
    Integer Function mlf_basis_inv_integ(this, x0, Y, rpar, X)
      use iso_c_binding
      import :: mlf_basis_fun_inv
      class(mlf_basis_fun_inv), intent(in), target :: this
      real(c_double), intent(in) :: x0
      real(c_double), intent(in), target :: Y(:), rpar(:)
      real(c_double), intent(out), target :: X(:)
    End Function mlf_basis_inv_integ

    ! Abstract basis function type (for dimension reduction)
    Integer Function mlf_basis_evalder(this, X, rpar, Y)
      use iso_c_binding
      import :: mlf_basis_funder
      class(mlf_basis_funder), intent(in), target :: this
      real(c_double), intent(in), target :: X(:), rpar(:,:)
      real(c_double), intent(out), target :: Y(:,:,:)
    End Function mlf_basis_evalder

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

    ! Functions definition for the bijection object
    Pure Subroutine mlf_bijection_proj(this, X, Y)
      use iso_c_binding
      import :: mlf_bijection_fun
      class(mlf_bijection_fun), intent(in) :: this
      real(c_double), intent(in) :: X(:)
      real(c_double), intent(out) :: Y(:)
    End Subroutine mlf_bijection_proj

    ! Functions definition for the bijection object
    Pure Subroutine mlf_bijection_invert(this, Y, X)
      use iso_c_binding
      import :: mlf_bijection_fun
      class(mlf_bijection_fun), intent(in) :: this
      real(c_double), intent(in) :: Y(:)
      real(c_double), intent(out) :: X(:)
    End Subroutine mlf_bijection_invert

   End Interface
Contains
  Integer Function mlf_ode_checkCstr(this, t, X0, Xd, K0, K1, ids, K) Result(N)
    class(mlf_ode_funCstr), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X0(:), Xd(:), K0(:), K1(:)
    integer, intent(inout), target :: ids(:)
    integer, intent(in) :: K
    ! Default behaviour: keep all the selected ids
    N = K
  End Function mlf_ode_checkCstr

  Subroutine mlf_ode_allocateCstr(this, N, M)
    class(mlf_ode_funCstrVect), intent(inout), target :: this
    integer, intent(in) :: N, M
    ALLOCATE(this%cstrVect(N,M), this%cstrValRef(M), this%cstrLastVal(M), &
      this%cstrLastDer(M), this%cstrTmp(M, 2), this%cstrAlpha(M))
    this%NCstr = M
    this%cstrValRef = 0
    this%cstrLastVal = 0
    this%cstrLastDer = 0
    this%cstrAlpha = 1.5d0
    this%cstrId = -1
    ! Set the round error of the dot-product
    this%epsilonT = 8d0*epsilon(1d0)*sqrt(real(N))
  End Subroutine mlf_ode_allocateCstr

  Subroutine mlf_ode_allocateCstrVect(this, Vect, ValRef)
    class(mlf_ode_funCstrVect), intent(inout), target :: this
    real(c_double) :: Vect(:,:)
    real(c_double), optional :: ValRef(:)
    CALL mlf_ode_allocateCstr(this, size(Vect,1), size(Vect,2))
    this%cstrVect = Vect
    if(PRESENT(ValRef)) this%cstrValRef = ValRef
  End Subroutine mlf_ode_allocateCstrVect

  Subroutine mlf_ode_allocateCstrIds(this, Ids, ValRef)
    class(mlf_ode_funCstrIds), intent(inout), target :: this
    integer, intent(in) :: Ids(:)
    real(c_double), optional :: ValRef(:)
    integer :: N
    N = size(Ids)
    ALLOCATE(this%cstrValRef(N), this%cstrAlpha(N))
    this%NCstr = N
    this%cstrValRef = 0
    this%cstrAlpha = 1.5d0
    this%cstrId = -1
    ALLOCATE(this%cstrSelIds, SOURCE = Ids)
    if(PRESENT(ValRef)) this%cstrValRef = ValRef
  End Subroutine mlf_ode_allocateCstrIds

  Real(c_double) Function GetHMaxFromCU(C, U, Alpha, id, T) result(hMax)
    real(c_double), intent(in) :: C(:), U(:), Alpha(:)
    real(c_double), intent(inout) :: T
    integer, intent(out) :: id
    integer :: i
    real(c_double) :: h
    id = -1; hMax = HUGE(hMax)
    Do i=1, size(C)
      If(U(i)*C(i)<0) Then
        h = -Alpha(i)*C(i)/U(i)
        If(h >= hMax) CYCLE
        hMax = h
        id = i
      Endif
    End Do
    T = T + hMax
  End Function GetHMaxFromCU

  ! Dummy function for function event
  Integer Function mlf_ode_dummy_funStop(this, t, ids, X, F) result(info)
    class(mlf_ode_funStop), intent(inout), target :: this
    real(c_double), intent(inout), target :: t, X(:), F(:)
    integer, intent(in) :: ids(:)
    info = mlf_ODE_SoftCstr
  End Function mlf_ode_dummy_funStop

  ! Default function that update the alpha if t is reached after cstrT
  ! Reevaluate dF/dt and move slightly t and X for reaching the condition
  ! The root should be near the value t (< tol)
  Integer Function mlf_ode_reachCstr(this, t, tMin, tMax, ids, X, F) result(info)
    class(mlf_ode_funCstrVect), intent(inout), target :: this
    real(c_double), intent(inout) :: t
    real(c_double), intent(in) :: tMin, tMax
    real(c_double), intent(inout), target :: X(:), F(:)
    integer, intent(in) :: ids(:)
    real(c_double) :: h, Uid
    integer :: id
    If(ANY(this%cstrId == ids) .AND. t > this%cstrT) Then
      this%cstrAlpha(this%cstrId) = 1.5d0*this%cstrAlpha(this%cstrId)
    Endif
    id = RandFromIArr(ids)
    this%cstrId = id
    ! Compute the value of the constraints
    Uid = DOT_PRODUCT(X, this%cstrVect(:,id))-this%cstrValRef(id)

    ! Compute h such as <X(t+h),cstrVect(:,id)> = 0
    ! As h is very small, X(t+h)=X(t)+h*F(t)+O(h²)
    ! So with h = -<X(t),cstrVect(:,id)>/<F(t),cstrVect(:,id)>
    ! we will have <X(t)+h*F(t),cstrVect(:,id)> = 0
    h = -Uid/DOT_PRODUCT(F, this%cstrVect(:,id))
    ! Add epsilon, so t will be slightly after collision
    h = h + abs(h)*this%epsilonT
    h = MIN(MAX(h, tMin-t), tMax-t)
    t = t + h
    X = X + h*F
    info = this%funStop(t, ids, X, F)
    ! Check if there is an error or a hard constraints
    If(info < 0 .OR. info == mlf_ODE_HardCstr .OR. info == mlf_ODE_StopTime) RETURN
    info = this%eval(t, X, F)
  End Function mlf_ode_reachCstr
  
  Integer Function mlf_ode_reachCstrIds(this, t, tMin, tMax, ids, X, F) result(info)
    class(mlf_ode_funCstrIds), intent(inout), target :: this
    real(c_double), intent(inout) :: t
    real(c_double), intent(in) :: tMin, tMax
    real(c_double), intent(inout), target :: X(:), F(:)
    integer, intent(in) :: ids(:)
    real(c_double) :: h, Uid
    integer :: k, id
    If(ANY(this%cstrId == ids) .AND. t > this%cstrT) Then
      this%cstrAlpha(this%cstrId) = 1.5d0*this%cstrAlpha(this%cstrId)
    Endif
    id = RandFromIArr(ids)
    k = this%cstrSelIds(id)
    this%cstrId = id
    ! Compute the value of the constraints
    Uid = X(k)-this%cstrValRef(id)
    h = MIN(MAX(-Uid/F(k), tMin-t), tMax-t)
    t = t + h
    X = X + h*F
    X(ids) = 0
    info = this%funStop(t, ids, X, F)
    ! Check if there is an error or a hard constraints
    If(info < 0 .OR. info == mlf_ODE_HardCstr .OR. info == mlf_ODE_StopTime) RETURN
    info = this%eval(t, X, F)
  End Function mlf_ode_reachCstrIds

  Integer Function SelectIdsCrossing(ids, X0, F0, X) result(n)
    real(c_double), intent(in) :: X0(:), F0(:), X(:)
    integer, intent(inout) :: ids(:)
    real(c_double) :: U
    integer :: i
    n = 0
    Do i = 1,size(X0)
      If(X(i) /= 0) then
        U = X0(i)*X(i)
        If(U > 0) CYCLE
        If(U == 0) then ! Initial point has reached X0 == 0
          U = F0(i)*X(i) ! Use the derivative instead to look if X == 0 has been crossed
          If(U >= 0) CYCLE
        Endif
      Endif ! else X == 0 has been reached
      n = n + 1
      ids(n) = i
    End Do
  End Function SelectIdsCrossing

  Subroutine mlf_ode_getDerivatives(this, ids, X0, X, K, C0, C, Q)
    class(mlf_ode_funCstrVect), intent(inout), target :: this
    real(c_double), intent(in), target :: K(:,:), X0(:), X(:)
    real(c_double), intent(out), target :: C0(:), C(:), Q(:,:)
    integer, intent(in), target :: ids(:)
    integer :: N, i
    N = size(ids)
    C0(1:N) = this%cstrLastVal(ids)
    C(1:N) = this%cstrTmp(ids,1)
    Do i = 1,N
      Q(i,:) = MATMUL(this%cstrVect(:,ids(i)), K)
    End Do
    ! Q(1:N,:) = MATMUL(TRANSPOSE(this%cstrVect(:,ids)), K)
  End Subroutine mlf_ode_getDerivatives

  Subroutine mlf_ode_getDerivativesIds(this, ids, X0, X, K, C0, C, Q)
    class(mlf_ode_funCstrIds), intent(inout), target :: this
    real(c_double), intent(in), target :: K(:,:), X0(:), X(:)
    real(c_double), intent(out), target :: C0(:), C(:), Q(:,:)
    integer, intent(in), target :: ids(:)
    integer :: N
    N = size(ids)
    C0(1:N) = X0(ids)
    C(1:N) = X(ids)
    Q = K(this%cstrSelIds(ids),:)
  End Subroutine mlf_ode_getDerivativesIds

  Real(c_double) Function mlf_ode_getHMaxIds(this, t, X, F) Result(hMax)
    class(mlf_ode_funCstrIds), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:), F(:)
    integer :: i, j
    real(c_double) :: Z
    hMax = HUGE(hMax)
    Do i = 1, size(this%cstrSelIds)
      j = this%cstrSelIds(i)
      Z = X(j) - this%cstrValRef(i)
      If(F(j)*Z < 0) hMax = MIN(hMax, -this%cstrAlpha(i)*Z/F(j))
    End Do
    this%cstrT = t
  End Function mlf_ode_getHMaxIds

  Real(c_double) Function mlf_ode_getHMax(this, t, X, F) Result(hMax)
    class(mlf_ode_funCstrVect), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:), F(:)
    this%cstrLastVal = MATMUL(X, this%cstrVect) - this%cstrValRef
    this%cstrLastDer = MATMUL(F, this%cstrVect)
    this%cstrT = t
    hMax = GetHMaxFromCU(this%cstrLastVal, this%cstrLastDer, this%cstrAlpha, &
      this%cstrId, this%cstrT)
  End Function mlf_ode_getHMax

  ! Update function for single value constraints
  Integer Function mlf_ode_updateCstrIds(this, t, X0, X, F0, F, ids, hMax) &
      result(N)
    class(mlf_ode_funCstrIds), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X0(:), X(:), F0(:), F(:)
    integer, intent(out), target :: ids(:)
    real(c_double), intent(inout) :: hMax
    integer :: i, j
    real(c_double) :: U, V, Z
    this%cstrId = -1
    N = 0
    Do i = 1, size(this%cstrSelIds)
      j = this%cstrSelIds(i)
      V = this%cstrValRef(i)
      Z = X(j) - V
      If(Z /= 0) Then
        U = (X0(j)-V)*Z
        If(U > 0 .OR. (U == 0 .AND. F0(j)*Z >= 0)) Then
          If(F(j) /= 0) hMax = MIN(hMax, -this%cstrAlpha(i)*Z/F(j))
          CYCLE
        Endif
      Endif
      N = N + 1
      ids(N) = j
    End Do
    If(N>0) RETURN
    this%cstrT = t
  End Function mlf_ode_updateCstrIds

  ! Update function for vector constraints
  Integer Function mlf_ode_updateCstr(this, t, X0, X, F0, F, ids, hMax) &
      result(N)
    class(mlf_ode_funCstrVect), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X0(:), X(:), F0(:), F(:)
    integer, intent(out), target :: ids(:)
    real(c_double), intent(inout) :: hMax
    this%cstrId = -1
    this%cstrTmp(:,1) = MATMUL(X, this%cstrVect) - this%cstrValRef
    this%cstrTmp(:,2) = MATMUL(F, this%cstrVect)
    N = SelectIdsCrossing(ids, this%cstrLastVal, this%cstrLastDer, &
                          this%cstrTmp(:,1))
    If(N>0) RETURN
    this%cstrLastVal = this%cstrTmp(:,1)
    this%cstrLastDer = this%cstrTmp(:,2)
    this%cstrT = t
    hMax = GetHMaxFromCU(this%cstrLastVal, this%cstrLastDer, this%cstrAlpha, &
      this%cstrId, this%cstrT)
  End Function mlf_ode_updateCstr

  Integer Function mlf_basis_c_eval(this, X, rpar, Y) result(info)
    class(mlf_basis_fun_c), intent(in), target :: this
    real(c_double), intent(in), target :: X(:), rpar(:,:)
    real(c_double), intent(out), target :: Y(:,:)
    integer(c_int) :: nX, sPar, nPar
    If(.NOT. ASSOCIATED(this%evalC)) Then
      info = mlf_UNINIT
      RETURN
    Endif
    nX = size(X,1)
    sPar = size(rpar,1)
    nPar = size(rpar,2)
    info = this%evalC(C_LOC(X), C_LOC(rpar), C_LOC(Y), nX, nPar, sPar, this%ptr)
  End Function mlf_basis_c_eval

  Integer Function mlf_obj_c_eval(this, X, Y) result(info)
    class(mlf_objective_fun_c), intent(inout), target :: this
    real(c_double), intent(in), target :: X(:,:)
    real(c_double), intent(inout), target :: Y(:,:)
    integer(c_int) :: ND, NY, lambda
    If(.NOT. associated(this%evalC)) Then
      info = mlf_UNINIT
      RETURN
    Endif
    lambda = size(X,2)
    ND = size(X,1)
    NY = size(Y,1)
    info = this%evalC(C_LOC(X), C_LOC(Y), ND, nY, lambda, this%ptr)
  End Function mlf_obj_c_eval

  Type(c_ptr) Function c_objfunction(cfun, cptr, funinfo) bind(C, name="mlf_objfunction")
    type(c_ptr), value :: cptr
    type(mlf_objfuninfo) :: funinfo
    type(c_funptr), value :: cfun
    type(mlf_objective_fun_c), pointer :: x
    class (*), pointer :: obj
    ALLOCATE(x)
    x%ptr = cptr
    x%nD = funinfo%nDimIn; x%nY = funinfo%nDimOut; x%nC = funinfo%nDimCstr
    If(C_ASSOCIATED(cfun)) call C_F_PROCPOINTER(cfun, x%evalC)
    obj => x
    c_objfunction = c_allocate(obj)
  End Function c_objfunction

  Type(c_ptr) Function c_basisfunction(cfun, cptr) bind(C, name="mlf_basisfunction")
    type(c_ptr), value :: cptr
    type(c_funptr), value :: cfun
    type(mlf_basis_fun_c), pointer :: x
    class (*), pointer :: obj
    ALLOCATE(x)
    x%ptr = cptr
    If(C_ASSOCIATED(cfun)) call C_F_PROCPOINTER(cfun, x%evalC)
    obj => x
    c_basisfunction = c_allocate(obj)
  End Function c_basisfunction

  ! Abstract weight function type
  Subroutine mlf_weight_simple_ev(this, Y, lambda)
    class(mlf_weight_simple_fun), intent(in), target :: this
    real(c_double), intent(out) :: Y(:)
    integer :: lambda
    CALL this%evalS(Y, lambda)
  End Subroutine mlf_weight_simple_ev

End Module mlf_fun_intf
