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

  ! Function handler for basis functions
  Type, Public, abstract, extends(mlf_obj) :: mlf_basis_funder
  Contains
    procedure (mlf_basis_evalder), deferred :: eval
  End Type mlf_basis_funder


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
  Contains
    procedure (mlf_ode_funCstr_update), deferred :: updateCstr
    procedure (mlf_ode_funCstr_reach), deferred :: reachCstr
    procedure (mlf_ode_funCstr_getDerivatives), deferred :: getDerivatives
  End Type mlf_ode_funCstr

  Type, Public, abstract, extends(mlf_ode_funCstr) :: mlf_ode_funStop
  Contains
    procedure :: funStop => mlf_ode_dummy_funStop
  End Type mlf_ode_funStop

  Type, Public, abstract, extends(mlf_ode_funStop) :: mlf_ode_funCstrIds
    real(c_double), allocatable :: cstrValRef(:)
    real(c_double), allocatable :: cstrLastVal(:)
    real(c_double), allocatable :: cstrLastDer(:)
    real(c_double), allocatable :: cstrAlpha(:)
    real(c_double), allocatable :: cstrTmp(:)
    integer, allocatable :: cstrSelIds(:)
    integer, allocatable :: cstrIds(:)
    real(c_double) :: cstrT
    integer :: cstrId
  Contains
    procedure :: allocateCstr => mlf_ode_allocateCstrIds
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
    integer, allocatable :: cstrIds(:)
    real(c_double) :: cstrT, epsilonT
    integer :: cstrId
  Contains
    procedure :: mlf_ode_allocateCstr, mlf_ode_allocateCstrVect
    generic :: allocateCstr => mlf_ode_allocateCstr, mlf_ode_allocateCstrVect
    procedure :: updateCstr => mlf_ode_updateCstr
    procedure :: reachCstr => mlf_ode_reachCstr
    procedure :: getDerivatives => mlf_ode_getDerivatives
  End Type mlf_ode_funCstrVect

  Integer, Parameter, Public :: mlf_ODE_StopT = 2, mlf_ODE_SoftCstr = 3, mlf_ODE_HardCstr = 4

  Abstract Interface

    Integer Function mlf_ode_funCstr_reach(this, t, id, X, F, hMax)
      Use iso_c_binding
      import :: mlf_ode_funCstr
      class(mlf_ode_funCstr), intent(inout), target :: this
      real(c_double), intent(inout), target :: t, X(:), F(:), hMax
      integer, intent(in) :: id
    End Function mlf_ode_funCstr_reach

    Subroutine mlf_ode_funCstr_getDerivatives(this, ids, K, C0, C, Q)
      Use iso_c_binding
      import :: mlf_ode_funCstr
      class(mlf_ode_funCstr), intent(inout), target :: this
      real(c_double), intent(in), target :: K(:,:)
      real(c_double), intent(out), target :: C0(:), C(:), Q(:,:)
      integer, intent(inout), target :: ids(:)
    End Subroutine mlf_ode_funCstr_getDerivatives

    Function mlf_ode_funCstr_update(this, t, X, F, ids)
      Use iso_c_binding
      import :: mlf_ode_funCstr
      class(mlf_ode_funCstr), intent(inout), target :: this
      real(c_double), intent(in) :: t
      real(c_double), intent(in), target :: X(:), F(:)
      integer, intent(out), optional, pointer :: ids(:)
      real(c_double) :: mlf_ode_funCstr_update
    End Function mlf_ode_funCstr_update

    ! Abstract ODE function type
    Integer Function mlf_ode_eval(this, t, X, F)
      Use iso_c_binding
      import :: mlf_ode_fun
      class(mlf_ode_fun), intent(in), target :: this
      real(c_double), intent(in) :: t
      real(c_double), intent(in), target :: X(:)
      real(c_double), intent(out), target :: F(:)
    End Function mlf_ode_eval

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

    ! Abstract basis function type (for dimension reduction)
    integer Function mlf_basis_evalder(this, X, rpar, Y)
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
  End Interface
Contains

  Subroutine mlf_ode_allocateCstr(this, N, M)
    class(mlf_ode_funCstrVect), intent(inout), target :: this
    integer, intent(in) :: N, M
    ALLOCATE(this%cstrVect(N,M), this%cstrValRef(M), this%cstrLastVal(M), &
      this%cstrLastDer(M), this%cstrTmp(M, 2), this%cstrAlpha(M), this%cstrIds(M))
    this%cstrValRef = 0
    this%cstrLastVal = 0
    this%cstrLastDer = 0
    this%cstrAlpha = 1.5d0
    this%cstrId = -1
    this%epsilonT = 8d0*epsilon(1d0)*sqrt(real(N)) ! ~ Round error of the dot-product
  End Subroutine mlf_ode_allocateCstr

  Subroutine mlf_ode_allocateCstrVect(this, Vect, ValRef)
    class(mlf_ode_funCstrVect), intent(inout), target :: this
    real(c_double) :: Vect(:,:)
    real(c_double), optional :: ValRef(:)
    call mlf_ode_allocateCstr(this, size(Vect,1), size(Vect,2))
    this%cstrVect = Vect
    if(present(ValRef)) this%cstrValRef = ValRef
  End Subroutine mlf_ode_allocateCstrVect

  Subroutine mlf_ode_allocateCstrIds(this, Ids, ValRef)
    class(mlf_ode_funCstrIds), intent(inout), target :: this
    integer, intent(in) :: Ids(:)
    real(c_double), optional :: ValRef(:)
    integer :: N
    N = size(Ids)
    ALLOCATE(this%cstrValRef(N), this%cstrLastVal(N), this%cstrLastDer(N), &
      this%cstrAlpha(N), this%cstrSelIds(N), this%cstrIds(N), this%cstrTmp(N))
    this%cstrSelIds = Ids
    this%cstrValRef = 0
    this%cstrLastVal = 0
    this%cstrLastDer = 0
    this%cstrAlpha = 1.5d0
    this%cstrId = -1
    if(present(ValRef)) this%cstrValRef = ValRef
  End Subroutine mlf_ode_allocateCstrIds

  Real(c_double) Function GetHMaxFromCU(C, U, Alpha, id, T) result(hMax)
    real(c_double), intent(in) :: C(:), U(:), Alpha(:)
    real(c_double), intent(inout) :: T
    integer, intent(out) :: id
    integer :: i
    real(c_double) :: h
    id = -1; hMax = HUGE(hMax)
    Do i=1, size(C)
      if(U(i)*C(i)<0) then
        h = -Alpha(i)*C(i)/U(i)
        if(h >= hMax) CYCLE
        hMax = h
        id = i
      endif
    End Do
    T = T + hMax
  End Function GetHMaxFromCU

  ! Dummy function for function event
  Integer Function mlf_ode_dummy_funStop(this, t, id, X, F) result(info)
    class(mlf_ode_funStop), intent(inout), target :: this
    real(c_double), intent(inout), target :: t, X(:), F(:)
    integer, intent(in) :: id
    info = mlf_ODE_SoftCstr
  End Function mlf_ode_dummy_funStop

  ! Default function that update the alpha if t is reached after cstrT
  ! Reevaluate dF/dt and move slightly t and X for reaching the condition
  ! The root should be near the value t (< tol)
  Integer Function mlf_ode_reachCstr(this, t, id, X, F, hMax) result(info)
    class(mlf_ode_funCstrVect), intent(inout), target :: this
    real(c_double), intent(inout), target :: t, X(:), F(:), hMax
    integer, intent(in) :: id
    real(c_double) :: h, Uid
    if(this%cstrId == id .AND. t > this%cstrT) then
      this%cstrAlpha(id) = 1.5d0*this%cstrAlpha(id)
    endif
    this%cstrId = id
    ! Compute the value of the constraints
    Uid = DOT_PRODUCT(X, this%cstrVect(:,id))-this%cstrValRef(id)
    info = this%eval(t, X, F)
    if(info /= 0) RETURN ! If there is an error or a hard constraints
    ! Compute h such as <X(t+h),cstrVect(:,id)> = 0
    ! As h is very small, X(t+h)=X(t)+h*F(t)+O(h²)
    ! So with h = -<X(t),cstrVect(:,id)>/<F(t),cstrVect(:,id)>
    ! we will have <X(t)+h*F(t),cstrVect(:,id)> = 0
    h = -Uid/DOT_PRODUCT(F, this%cstrVect(:,id))
    h = h + abs(h)*this%epsilonT ! Add epsilon, so t will be slightly after collison
    t = t + h
    X = X + h*F
    info = this%funStop(t, id, X, F)
    this%cstrLastVal = MATMUL(X, this%cstrVect) - this%cstrValRef
    this%cstrLastDer = MATMUL(F, this%cstrVect)
    this%cstrT = t
    hMax = GetHMaxFromCU(this%cstrLastVal, this%cstrLastDer, this%cstrAlpha, this%cstrId, this%cstrT)
  End Function mlf_ode_reachCstr
  
  Integer Function mlf_ode_reachCstrIds(this, t, id, X, F, hMax) result(info)
    class(mlf_ode_funCstrIds), intent(inout), target :: this
    real(c_double), intent(inout), target :: t, X(:), F(:), hMax
    integer, intent(in) :: id
    real(c_double) :: h, Uid
    integer :: k
    if(this%cstrId == id .AND. t > this%cstrT) then
      this%cstrAlpha(id) = 1.5d0*this%cstrAlpha(id)
    endif
    k = this%cstrSelIds(id)
    this%cstrId = id
    ! Compute the value of the constraints
    Uid = X(k)-this%cstrValRef(id)
    info = this%eval(t, X, F)
    if(info /= 0) RETURN ! If there is an error or a hard constraints
    h = -Uid/F(k)
    t = t + h
    X = X + h*F
    X(k) = 0
    info = this%funStop(t, id, X, F)
    this%cstrLastVal = X(this%cstrSelIds) - this%cstrValRef
    this%cstrLastDer = F(this%cstrSelIds)
    this%cstrT = t
    hMax = GetHMaxFromCU(this%cstrLastVal, this%cstrLastDer, this%cstrAlpha, this%cstrId, this%cstrT)
  End Function mlf_ode_reachCstrIds

  Integer Function SelectIdsCrossing(ids, X0, F0, X) result(n)
    real(c_double), intent(in) :: X0(:), F0(:), X(:)
    integer, intent(inout) :: ids(:)
    real(c_double) :: U
    integer :: i
    n = 0
    Do i = 1,size(X0)
      if(X(i) /= 0) then
        U = X0(i)*X(i)
        if(U > 0) CYCLE
        if(U == 0) then ! Initial point has reached X0 == 0
          U = F0(i)*X(i) ! Use the derivative instead to look if X == 0 has been crossed
          if(U >= 0) CYCLE
        endif
      endif ! else X == 0 has been reached
      n = n + 1
      ids(n) = i
    End Do
  End Function SelectIdsCrossing

  Subroutine mlf_ode_getDerivatives(this, ids, K, C0, C, Q)
    class(mlf_ode_funCstrVect), intent(inout), target :: this
    real(c_double), intent(in), target :: K(:,:)
    real(c_double), intent(out), target :: C0(:), C(:), Q(:,:)
    integer, intent(inout), target :: ids(:)
    integer :: N
    N = size(ids)
    C0(1:N) = this%cstrLastVal(ids)
    C(1:N) = this%cstrTmp(ids,1)
    Q = MATMUL(TRANSPOSE(this%cstrVect(:,ids)), K)
  End Subroutine mlf_ode_getDerivatives

  Subroutine mlf_ode_getDerivativesIds(this, ids, K, C0, C, Q)
    class(mlf_ode_funCstrIds), intent(inout), target :: this
    real(c_double), intent(in), target :: K(:,:)
    real(c_double), intent(out), target :: C0(:), C(:), Q(:,:)
    integer, intent(inout), target :: ids(:)
    integer :: N
    N = size(ids)
    C0(1:N) = this%cstrLastVal(ids)
    C(1:N) = this%cstrTmp(ids)
    Q = K(this%cstrSelIds(ids),:)
  End Subroutine mlf_ode_getDerivativesIds

  ! Update function for single value constraints
  Real(c_double) Function mlf_ode_updateCstrIds(this, t, X, F, ids) result(hMax)
    class(mlf_ode_funCstrIds), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:), F(:)
    integer, intent(out), optional, pointer :: ids(:)
    integer :: i
    hmax = HUGE(hMax)
    this%cstrId = -1
    if(PRESENT(ids)) then
      ids => NULL()
      this%cstrTmp = X(this%cstrSelIds) - this%cstrValRef
      i = SelectIdsCrossing(this%cstrIds, this%cstrLastVal, this%cstrLastDer, &
                            this%cstrTmp)
      if(i > 0) then
        ids => this%cstrIds(1:i)
        RETURN
      endif
    endif
    this%cstrLastVal = X(this%cstrSelIds) - this%cstrValRef
    this%cstrLastDer = F(this%cstrSelIds)
    this%cstrT = t
    hMax = GetHMaxFromCU(this%cstrLastVal, this%cstrLastDer, this%cstrAlpha, this%cstrId, this%cstrT)
  End Function mlf_ode_updateCstrIds

  ! Update function for vector constraints
  Real(c_double) Function mlf_ode_updateCstr(this, t, X, F, ids) result(hMax)
    class(mlf_ode_funCstrVect), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:), F(:)
    integer, intent(out), optional, pointer :: ids(:)
    integer :: i
    hmax = HUGE(hMax)
    this%cstrId = -1
    if(PRESENT(ids)) then
      ids => NULL()
      this%cstrTmp(:,1) = MATMUL(X, this%cstrVect) - this%cstrValRef
      this%cstrTmp(:,2) = MATMUL(F, this%cstrVect)
      i = SelectIdsCrossing(this%cstrIds, this%cstrLastVal, this%cstrLastDer, &
                            this%cstrTmp(:,1))
      if(i > 0) then
        ids => this%cstrIds(1:i)
        RETURN
      endif
      this%cstrLastVal = this%cstrTmp(:,1)
      this%cstrLastDer = this%cstrTmp(:,2)
    else
      this%cstrLastVal = MATMUL(X, this%cstrVect) - this%cstrValRef
      this%cstrLastDer = MATMUL(F, this%cstrVect)
    endif
    this%cstrT = t
    hMax = GetHMaxFromCU(this%cstrLastVal, this%cstrLastDer, this%cstrAlpha, this%cstrId, this%cstrT)
  End Function mlf_ode_updateCstr

  Integer Function mlf_basis_c_eval(this, X, rpar, Y) result(info)
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

  Integer Function mlf_obj_c_eval(this, X, Y) result(info)
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

  Integer Function mlf_obj_c_constraints(this, X, Y) result(info)
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
