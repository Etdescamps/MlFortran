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

Module mlf_poly
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_cfuns
  Use mlf_utils
  IMPLICIT NONE
  PRIVATE

  Public :: mlf_solve4DPoly, mlf_rootsPoly, mlf_polyFromRoots, mlf_polyVal
  Public :: FindRoot2DInterp, FindRoot2PDerivative

Contains

  Elemental Real(c_double) Function FindRoot2DInterp(x0, x1, dx1) Result(y)
    ! Evaluate the root of the 2nd order polynomial f with:
    ! f(0) = x0; f(1) = x1; f'(1) = dx1
    real(c_double), intent(in) :: x0, x1, dx1
    real(c_double) :: a, b, z
    a = dx1 - x1 + x0
    If(ABS(a) < 1d-4*ABS(dx1)) Then
      ! 1D interpolation (ignore derivate at 1)
      y = x0/(x0-x1)
      RETURN
    Endif
    b = 1d0 - 0.5d0*dx1/a
    z = SQRT(b**2-x0/a)
    y = b + z ! Choose the root closer to one (if < 1)
    If(y < 1d0) RETURN
    y = b - z
  End Function FindRoot2DInterp

  Integer Function ZRootsQuadratic(R, b, c, eps) Result(nr)
    complex(c_double), intent(in) :: b, c
    complex(c_double), intent(out) :: R(:)
    real(c_double), intent(in) :: eps
    complex(c_double) :: delta
    nr = 0
    delta = 0.25d0*b**2-c
    If(ABS(delta) < eps*ABS(c)) Then
      R(1) = -0.5d0*b
      nr = 1
      RETURN
    Endif
    delta = SQRT(delta)
    R(1) = -0.5d0*b-delta
    R(2) = -0.5d0*b+delta
    nr = 2
  End Function ZRootsQuadratic

  ! Find the solutions of the equation z^3+p*z+q=0
  Integer Function ZRootsCardano(R, p, q, eps) Result(nr)
    complex(c_double), intent(in) :: p, q
    complex(c_double), intent(out) :: R(:)
    real(c_double), intent(in) :: eps
    complex(c_double) :: delta
    complex(c_double) :: a
    integer :: i
    complex(c_double), parameter :: j = (-0.5d0, 8.660254037844386467869d-1)
    nr = 1 ! At least one solution
    If(ABS(q)<eps) Then
      R(1) = (0d0, 0d0)
      If(ABS(p)<eps) RETURN
      nr = 3
      a = SQRT(-p)
      R(2:3) = [a, -a]
      RETURN
    Endif
    If(abs(p)<eps) Then
      nr = 3
      a = p**(1d0/3d0)
      R(1:3) = [a, a*j, a*j*j]
    Endif
    delta = 1d0/27d0*p**31d0/4d0*q**2
    delta = SQRT(delta)
    a = (delta-0.5d0*q)**(1d0/3d0)
    Do i=1,3
      R(i) = a-p/(3d0*a)
      a = a*j
    End Do
    nr = 3
    ! Remove double solutions
    If(ABS(R(1)-R(3))<eps .OR. ABS(R(2)-R(3))<eps) Then
      nr = 2
    Else If(ABS(R(1)-R(2))<eps) Then
      R(2) = R(3)
      nr = 2
    Endif
  End Function ZRootsCardano

  Integer Function DRootsQuadratic(R, b, c, eps) Result(nr)
    real(c_double), intent(in) :: b, c, eps
    real(c_double), intent(out) :: R(:)
    real(c_double) :: delta
    nr = 0
    delta = 0.25d0*b**2-c
    If(ABS(delta) < eps*ABS(c)) Then
      R(1) = -0.5d0*b
      nr = 1
      RETURN
    Endif
    If(delta<0d0) RETURN
    delta = SQRT(delta)
    R(1) = -0.5d0*b-delta
    R(2) = -0.5d0*b+delta
    nr = 2
  End Function DRootsQuadratic

  ! Find the real solutions of the equation x^3+p*x+q=0
  Integer Function DRootsCardano(R, p, q, eps) Result(nr)
    real(c_double), intent(in) :: p, q, eps
    real(c_double), intent(out) :: R(:)
    real(c_double) :: delta
    complex(c_double) :: a
    complex(c_double), parameter :: j = (-0.5d0, 8.660254037844386467869d-1)
    nr = 1 ! At least one solution
    delta = -4d0*p**3-27d0*q**2
    If(ABS(delta)<eps*ABS(p+q)) Then ! discriminant = 0
      If(ABS(p) < eps) Then ! p == q == 0, so only the solution z=0
        R(1) = 0
      Else
        R(1) = 3d0*q/p
        R(2) = (-3d0*q)/(2d0*p) ! The double root of the equation
        ! So in this case: x^3+p*x+q=(x-R(1)*(x-R(2))^2
        nr = 2
      Endif
      RETURN
    Endif
    If(delta<0) Then ! One real solution (+ 2 complex solutions)
      delta = SQRT(-delta/27d0)
      R(1) = c_cbrt(0.5d0*(-q+delta))+c_cbrt(0.5d0*(-q-delta))
    Else ! Three real solutions
      nr = 3
      a = CMPLX(-0.5d0*q, 0.5*SQRT(delta/27d0), KIND=c_double)**(1d0/3d0)
      ! Get one complex root of (-q+i*sqrt(delta/27))/2
      ! the other complex roots are a*j and a*j^2
      ! The solutions have the form a*j^k+conjg(a*j^k) = 2*real(a*j^k)
      R(1) = 2d0*REAL(a, kind= c_double)
      a = a*j
      R(2) = 2d0*REAL(a, kind= c_double)
      a = a*j
      R(3) = 2d0*REAL(a, kind= c_double)
    Endif
  End Function DRootsCardano

  Real(c_double) Function FindRoot2PTaylorInterpolation(y0, d0, y1, d1) Result(x)
    real(c_double), intent(in) :: y0, d0, y1, d1
    real(c_double) :: a, b, c, d, R(3), eps
    integer :: i, nR
    eps = 1d-30*(MAX(ABS(y0), ABS(y1))/MAX(ABS(d0), ABS(d1)))
    a = d1 - 2*(y1-y0)
    b = 3*(y1-y0) - d1 - d0
    If(ABS(a) < eps) Then
      If(ABS(b) < 1d2*ABS(a)) Then
        x = y0/(y0-y1)
        RETURN
      Endif
      nR = DRootsQuadratic(R, d0/b, y0/b, eps)
    Else
      b = b/a; c = d0/a; d = y0/a
      nR = DrootsCardano(R, c-(b**2)/3d0, b/27d0*(2d0*b**2-9d0*c)+d, eps)
      R(1:nR) = R(1:nR) - b/3d0
    Endif
    Do i = 1,nR
      x = R(i)
      If(x > 0 .AND. x < 1) RETURN
    End Do
  End Function FindRoot2PTaylorInterpolation

  Real(c_double) Function FindRoot2PDerivative(y0, d0, y1, d1) Result(x)
    ! Approximate root of a function where f(0)=y0; f'(0)=d0; f(1)=y1; f'(1)=d1
    real(c_double), intent(in) :: y0, d0, y1, d1
    real(c_double) :: dy, m0, m1
    real(c_double) :: a, b, t
    dy = y1-y0
    If(y0 == 0) Then
      x = 0
    Else If(y1 == 1) Then
      x = 1
    Else If(SIGN(y1,y0) == y1) Then
      ! This case is not supported by the algorithm
      x = IEEE_VALUE(x, IEEE_QUIET_NAN)
    Else If(d0 == 0 .AND. d1 == 0) Then
      ! No way to do a proper interpolation so use the secant method instead
      x = -y0/dy
    Else
      ! Main part of the algorithm
      If(ABS(d0) == HUGE(1d0) .OR. ABS(d1) == HUGE(1d0)) Then
        m0 = HUGE(1d0)
      Else
        m0 = MAX(ABS(d0/dy), ABS(d1/dy))
      Endif
      If(d0 == 0d0 .OR. d1 == 0d0) Then
        m1 = HUGE(1d0)
      Else
        m1 = MAX(ABS(dy/d0), ABS(dy/d1))
      Endif
      If(MIN(m0, m1) > 1d8) Then
        ! Use quadrant approximation
        If(ABS(d0/dy) > 1d0) Then
          x = 1-SQRT(1-(y0/dy)**2)
        Else
          x = SQRT(1-(y1/dy)**2)
        Endif
      Else If(m0 < m1) Then
        ! Interpolate p(x) = y with p(0)=y0; p'(0)=d0; p(1)=y1; p'(1)=d1
        x = FindRoot2PTaylorInterpolation(y0, d0, y1, d1)
      Else
        ! Interpolate p(y) = x with p(0)=0; p'(0)=dy/d0; p(1)=1; p'(1)=dy/d1
        ! Rotate, translate and distord the function in order to fit in the square [0;1]²
        m0 = dy/d0
        m1 = dy/d1
        a = m1 - 2
        b = 3 - m1 - m0
        t = -y0/dy
        ! x is the evaluation at p(y0/dy)
        x = t*(m0+t*(b+t*a))
      Endif
    Endif
  End Function FindRoot2PDerivative

  ! Roots of a depressed quartic using Ferrari's method
  ! TODO: Debug!!
  ! Method too unstable for actual computation! Use DRootsCompanion instead
  Integer Function DRootsFerrari(A, p, q, r, eps) result(nr)
    real(c_double), intent(in) :: p, q, r, eps
    real(c_double), intent(out) :: A(:)
    real(c_double) :: u, v, mv, R2(3), m
    integer :: i, j
    nr = 0
    mv = max(abs(p), abs(q), abs(r))
    if(abs(q) < eps*mv) then
      ! Solve x**4+p*x**2+r=0
      j = DRootsQuadratic(R2, p, r, eps)
      Do i=1,j
        if(R2(i) < eps*mv) CYCLE
        if(R2(i) < eps*mv) then
          nr = nr+1
          A(nr) = 0d0
        else
          nr = nr+2
          u = sqrt(R2(i))
          A(nr-1:nr) = [-u, u]
        endif
      End Do
      RETURN
    endif
    ASSOCIATE(b => p, c => (0.25d0*p**2-r), d => (1d0/8d0*q**2))
      ! Equation P(x)=x^3+b*x^2+c*x+d is transformed in depressed form:
      ! x is replaced by z = x+b/3, so the equation has the form z^3+p*z+q
      j = DRootsCardano(R2, c-(b**2)/3d0, b/27d0*(2d0*b**2-9d0*c)+d, eps)
      R2(1:j) = R2(1:j) - b/3d0 ! x = z-b/3
    END ASSOCIATE
    m = Maxval(R2(1:j))
    if(m<-1d-12*mv) RETURN ! No solutions
    if(m<1d-12*mv) then ! Shall not happen (q /= 0)
      nr = DRootsQuadratic(A, 0d0, 0.5d0*p, eps)
      RETURN
    endif
    u = sqrt(2d0*m)
    v = 0.5d0*q/u
    nr = DRootsQuadratic(A, u, 0.5d0*p+m-v, eps)
    nr = nr + DRootsQuadratic(A(nr+1:), -u, 0.5d0*p+m+v, eps)
  End Function DRootsFerrari

  ! Classical method for solving the root of a polynomial using a companion Matrix
  Integer Function DRootsCompanion(A, R, eps) result(nr)
    real(c_double), intent(in) :: A(:), eps
    real(c_double), intent(out) :: R(:)
    real(c_double), allocatable :: M(:,:), Work(:), WR(:), WI(:)
    integer :: N, LWORK, INFO, i
    nr = 0
    N = size(A)
    LWORK = N*max(N,11)
    allocate(M(N,N), Work(LWORK), WR(N), WI(N))
    M = 0
    forall(i=1:N-1) M(i+1,i) = 1d0
    M(:,N) = -A
    !call DGEEV('N', 'N', N, M, N, WR, WI, C_NULL_PTR, N, C_NULL_PTR, N, Work, LWORK, INFO)
    call DHSEQR('E', 'N', N, 1, N, M, N, WR, WI, C_NULL_PTR, N, WORK, LWORK, INFO)
    if(info < 0) RETURN
    Do i=N,info+1,-1
      if(abs(WI(i))<eps) then
        nr = nr +1
        R(nr) = WR(i)
      endif
    End Do
  End Function DRootsCompanion

  ! Find real roots of a polynomial
  Integer Function DRootsPoly(A, R, eps) result(nr)
    real(c_double), intent(in) :: A(:), eps
    real(c_double), intent(out) :: R(:)
    real(c_double) :: delta
    integer :: N
    N = size(A)
    nr = 0
    select case (N)
      case(0)
        RETURN
      case(1)
        R(1) = -A(1)
        nr = 1
      case(2) ! Solve analytically the second order equation
        nr = DRootsQuadratic(R, A(2), A(1), eps)
      case(3) ! Use Cardano's method for 3rd order equations
        ASSOCIATE(b => A(3), c => A(2), d => A(1))
          ! Equation P(x)=x^3+b*x^2+c*x+d is transformed in depressed form:
          ! x is replaced by z = x+b/3, so the equation has the form z^3+p*z+q
          nr = DrootsCardano(R, c-(b**2)/3d0, b/27d0*(2d0*b**2-9d0*c)+d, eps)
          R(1:nr) = R(1:nr) - b/3d0 ! x = z-b/3
        END ASSOCIATE
      ! TODO : debug Ferrari's method
      !case(4) ! Use Ferrari's method
      !  ASSOCIATE(b => A(4), c => A(3), d => A(2), e => A(1))
      !    nr = rootsFerrari(R, c-3d0/8d0*b**2, 1d0/8d0*b**3-0.5*b*c+d, &
      !      e-3d0/256d0*b**4-0.25*b*d+1d0/16d0*b**2*c)
      !    R(1:nr) = R(1:nr) - b/4d0
      !  END ASSOCIATE
      case default
      ! Use companion matrix to determine polynomial
      nr = DRootsCompanion(A, R, eps)
    end select
  End Function DRootsPoly

  ! Interface for rootsPoly, check if the top coefficent of A can be negleted
  Integer Function mlf_rootsPoly(A, R, eps) result(nr)
    real(c_double), intent(in) :: A(:), eps
    real(c_double), intent(out) :: R(:)
    real(c_double) :: m
    integer :: i, j, N
    m = maxval(abs(A))
    N = size(A)
    nr = 0
    ! Remove head and tail null elements
    Do i=N,1,-1
      if(abs(A(i)) > m*eps) EXIT
    End Do
    Do j=1,i ! Simplify by X if a_j = 0 (from X^j)
      if(abs(A(j)) > abs(A(i))*eps) EXIT
    End Do
    if(i>1) nr = DRootsPoly(A(j:i-1)/A(i), R, eps)
    if(j>1) then ! If there is null element at tail -> X=0 is a solution
      nr = nr +1
      R(nr) = 0d0
    endif
  End Function mlf_rootsPoly

  real(c_double) Function mlf_solve4DPoly(t, a0, a1, a2, a3, a4) result(x)
    real(c_double), intent(in) :: t, a0, a1, a2, a3, a4
    real(c_double) :: R(4)
    integer :: nr, i, j
    nr = mlf_rootsPoly([a0, a1, a2, a3, a4], R, 1d-12)
    Do i=1,nr
      if(R(i)>t) EXIT
    End Do
    Do j=i+1,nr
      if(R(j)>t .AND. R(j)<R(i)) i = j
    End Do
    if(i>0) then
      x = R(i)
    else
      x = IEEE_VALUE(x, IEEE_QUIET_NAN)
    endif
  End Function mlf_solve4DPoly

  Pure real(c_double) Function mlf_polyVal(P, x) result(Y)
    real(c_double), intent(in) :: P(:), x
    integer :: N, i
    N =size(P)
    Y = P(N)
    Do i=N-1,1,-1
      Y = x*Y+P(i)
    End Do
  End Function mlf_polyVal

  Subroutine mlf_polyFromRoots(R, P)
    real(c_double), intent(in) :: R(:)
    real(c_double), intent(out) :: P(:)
    real(c_double) :: a
    integer :: k, N
    N = size(R)
    if(N == 1) then
      P(1) = -R(1)
      P(2) = 1d0
      RETURN
    endif
    P(1) = R(1)*R(2)
    P(2) = -R(1)-R(2)
    Do k=3,N
      a = R(k)
      P(k) = P(k-1)-a
      P(2:k-1) = P(1:k-2)-a*P(2:k-1)
      P(1) = -a*P(1)
    End Do
    P(N+1) = 1d0
  End Subroutine mlf_polyFromRoots


End Module mlf_poly

