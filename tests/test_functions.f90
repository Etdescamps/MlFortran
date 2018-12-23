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

Module test_functions
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_utils
  Use mlf_fun_intf
  IMPLICIT NONE
  
  ! Test function wrapper for fun_basis
  Type, extends(mlf_basis_fun) :: mlf_basis_test
    procedure (mlf_basis_test_fun), nopass, pointer :: evalT => NULL()
  Contains
    procedure :: eval => mlf_basis_test_eval
  End Type mlf_basis_test

  ! Test function wrapper for optim_algo
  Type, extends(mlf_objective_fun) :: mlf_objective_test
    procedure (mlf_obj_test_fun), nopass, pointer :: evalT => NULL()
    procedure (mlf_obj_test_fun), nopass, pointer :: constraintsT => NULL()
  Contains
    procedure :: eval => mlf_obj_test_eval, constraints => mlf_obj_test_constraints
    procedure :: init => mlf_obj_test_init
  End Type mlf_objective_test

  ! Test function ODE
  Type, extends(mlf_ode_fun) :: mlf_odeTest
    procedure (odeTest_fun), nopass, pointer :: evalF => NULL()
  Contains
    procedure :: eval => odeTest_eval
    procedure :: init => odeTest_init
  End Type mlf_odeTest

  ! Test function ODE
  Type, extends(mlf_ode_funCstrVect) :: mlf_odeTestCstr
    procedure (odeTest_fun), nopass, pointer :: evalF => NULL()
  Contains
    procedure :: eval => odeTestCstr_eval
    procedure :: init => odeTestCstr_init
  End Type mlf_odeTestCstr

  ! Test function ODE
  Type, extends(mlf_ode_funCstrIds) :: mlf_odeTestCstrIds
    procedure (odeTest_fun), nopass, pointer :: evalF => NULL()
  Contains
    procedure :: eval => odeTestCstrIds_eval
    procedure :: init => odeTestCstrIds_init
  End Type mlf_odeTestCstrIds


  real(c_double) , Parameter :: X0Arenstorf(4) = [0.994d0, 0d0, 0d0, -2.00158510637908252240537862224d0]
  real(c_double) , Parameter :: TEndArenstorf = 17.0652165601579625588917206249d0

  Abstract Interface
    ! Test basis function type (for dimension reduction)
    Subroutine mlf_basis_test_fun(X, rpar, Y)
      use iso_c_binding
      real(c_double), intent(in) :: X(:), rpar(:,:)
      real(c_double), intent(out) :: Y(:,:)
    End Subroutine mlf_basis_test_fun
    ! Test objective function type (for optimisation algorithms)
    Pure Function mlf_obj_test_fun(X) result(Y)
      use iso_c_binding
      real(c_double), intent(in) :: X(:)
      real(c_double) :: Y
    End Function mlf_obj_test_fun

    ! Test function for ODE
    Subroutine odeTest_fun(t, X, F)
      use iso_c_binding
      real(c_double), intent(in) :: X(:), t
      real(c_double), intent(out) :: F(:)
    End Subroutine odeTest_fun
  End Interface

Contains
  Integer Function odeTest_eval(this, t, X, F) result(info)
    class(mlf_odeTest), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:)
    real(c_double), intent(out), target :: F(:)
    call this%evalF(t, X, F)
    info = 0
  End Function odeTest_eval

  Subroutine odeTest_init(this, fun)
    class(mlf_odeTest), intent(inout), target :: this
    procedure (odeTest_fun) :: fun
    this%evalF => fun
  End Subroutine odeTest_init

  Integer Function odeTestCstrIds_eval(this, t, X, F) result(info)
    class(mlf_odeTestCstrIds), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:)
    real(c_double), intent(out), target :: F(:)
    call this%evalF(t, X, F)
    info = 0
  End Function odeTestCstrIds_eval

  Subroutine odeTestCstrIds_init(this, fun, Ids, ValRef)
    class(mlf_odeTestCstrIds), intent(inout), target :: this
    procedure (odeTest_fun) :: fun
    integer, intent(in) :: Ids(:)
    real(c_double), optional :: ValRef(:)
    call this%allocateCstr(Ids, ValRef)
    this%evalF => fun
  End Subroutine odeTestCstrIds_init

  Integer Function odeTestCstr_eval(this, t, X, F) result(info)
    class(mlf_odeTestCstr), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(in), target :: X(:)
    real(c_double), intent(out), target :: F(:)
    call this%evalF(t, X, F)
    info = 0
  End Function odeTestCstr_eval

  Subroutine odeTestCstr_init(this, fun, vC, ValRef)
    class(mlf_odeTestCstr), intent(inout), target :: this
    procedure (odeTest_fun) :: fun
    real(c_double), intent(in) :: vC(:,:)
    real(c_double), optional :: ValRef(:)
    call this%allocateCstr(vC, ValRef)
    this%evalF => fun
  End Subroutine odeTestCstr_init

  Subroutine mlf_obj_test_init(this, evalF, nD, cstrF, nC)
    class(mlf_objective_test), intent(inout), target :: this
    procedure(mlf_obj_test_fun) :: evalF
    procedure(mlf_obj_test_fun), optional :: cstrF
    integer :: nD
    integer, optional :: nC
    this%nD = nD; this%nC = 0; this%nY = 1
    this%evalT => evalF
    if(present(nC)) this%nC = nC
    if(present(cstrF)) this%constraintsT => cstrF
  End Subroutine mlf_obj_test_init

  integer Function mlf_basis_test_eval(this, X, rpar, Y) result(info)
    class(mlf_basis_test), intent(in), target :: this
    real(c_double), intent(in), target :: X(:), rpar(:,:)
    real(c_double), intent(out), target :: Y(:,:)
    info = 0
    call this%evalT(X, rpar, Y)
  End Function mlf_basis_test_eval
  
  integer Function mlf_obj_test_eval(this, X, Y) result(info)
    class(mlf_objective_test), intent(in), target :: this
    real(c_double), intent(in), target :: X(:,:)
    real(c_double), intent(inout), target :: Y(:,:)
    integer :: i
    info = 0
    Do i=1,size(X,2)
      Y(1,i) = this%evalT(X(:,i))
    End Do
  End Function mlf_obj_test_eval

  integer Function mlf_obj_test_constraints(this, X, Y) result(info)
    class(mlf_objective_test), intent(in), target :: this
    real(c_double), intent(in), target :: X(:,:)
    real(c_double), intent(inout), target :: Y(:,:)
    integer :: i
    info = 0
    if(associated(this%constraintsT)) then
      Do i=1,size(X,2)
        Y(1,i) = this%constraintsT(X(:,i))
      End Do
    else
      Y=0
    endif
  End Function mlf_obj_test_constraints

  Subroutine FArenstorf(t, X, F)
    real(c_double), intent(in) :: X(:), t
    real(c_double), intent(out) :: F(:)
    real(c_double), parameter :: mu = 0.012277471d0
    real(c_double) :: D1, D2, mu1
    mu1 = 1d0-mu
    ASSOCIATE(y1 => X(1), y1d => X(2), y2 => X(3), y2d => X(4), &
        dy1 => F(1), dy1d => F(2), dy2 => F(3), dy2d => F(4))
      dy1 = y1d; dy2 = y2d
      D1 = ((y1+mu)**2+y2**2)**(3d0/2d0)
      D2 = ((y1-mu1)**2+y2**2)**(3d0/2d0)
      dy1d = y1+2d0*y2d-mu1*(y1+mu)/D1-mu*(y1-mu1)/D2
      dy2d = y2-2d0*y1d-mu1*y2/D1-mu*y2/D2
    END ASSOCIATE
  End Subroutine FArenstorf

  ! F_a,b (x) = 1-(1/(1+a*exp(-b*x))-1)**2
  Subroutine FExpInv(X, rpar, Y)
    real(c_double), intent(in) :: X(:), rpar(:,:)
    real(c_double), intent(out) :: Y(:,:)
    integer :: N, M, i
    real(c_double) :: a, b
    N = size(rpar,2)
    M = size(X,1)
    Do i=1,N
      a = rpar(1,i); b = rpar(2,i)
      Y(:,i) = 1d0-(1d0/(1d0+a*exp(-b*X))-1d0)**2
    End Do
  End Subroutine FExpInv

  pure real(c_double) function OFsphere(X)
    real(c_double), intent(in) :: X(:)
    OFsphere = sum(X*X)
  end function OFsphere
  pure real(c_double) function OFcigar(X)
    real(c_double), intent(in) :: X(:)
    OFcigar = X(1)*X(1) + 1.0d6 * sum(X(2:)*X(2:))
  end function OFcigar
  pure real(c_double) function OFtablet(X)
    real(c_double), intent(in) :: X(:)
    OFtablet = 1.0d6 * X(1)*X(1) + sum(X(2:)*X(2:))
  end function OFtablet
  pure real(c_double) function OFellipsoid(X)
    real(c_double), intent(in) :: X(:)
    integer :: i
    OFellipsoid = sum( (/ (10**((i-1.0)/(size(X, 1)-1.0)), i=1,size(X, 1)) /)*X*X )
  end function OFellipsoid
  pure real(c_double) function OFparabolicRidge(X)
    real(c_double), intent(in) :: X(:)
    OFparabolicRidge = -X(1) +100.0*sum(X(2:)*X(2:))
  end function OFparabolicRidge
  pure real(c_double) function OFsharpRidge(X)
    real(c_double), intent(in) :: X(:)
    OFsharpRidge = -X(1) +100.0*norm2(X(2:))
  end function OFsharpRidge
  pure real(c_double) function OFrosenbrock(X)
    real(c_double), intent(in) :: X(:)
    real(c_double) :: U(size(X,1)-1)
    U(:) = 10d0*(X(2:)-X(:(size(X,1)-1))*X(:(size(X,1)-1)))
    OFrosenbrock = sum(U*U) 
    U(:) = X(:size(X,1)-1)-1d0
    OFrosenbrock = OFrosenbrock + sum(U*U)
  end function OFrosenbrock
End Module test_functions

