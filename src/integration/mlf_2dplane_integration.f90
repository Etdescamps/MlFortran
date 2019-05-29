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

Module mlf_2dplane_integration
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_models
  Use mlf_errors
  IMPLICIT NONE
  PRIVATE

  Type, Public, Abstract :: mlf_2dplane_h_fun
  Contains
    procedure(mlf_2dplane_h_getValDer), deferred :: getValDer
    procedure :: integrateOnPlane => mlf_2dplane_h_integrateOnPlane
  End Type mlf_2dplane_h_fun

  Type, Public :: mlf_2d_h_val
    real(c_double) :: val
    real(c_double) :: der(2)
    real(c_double) :: hes(3)
  End Type mlf_2d_h_val

  Abstract Interface
    Subroutine mlf_2dplane_h_getValDer(this, x, y, D)
      Use iso_c_binding
      import :: mlf_2dplane_h_fun, mlf_2d_h_val
      class(mlf_2dplane_h_fun), intent(inout) :: this
      real(c_double), intent(in) :: x, y
      type(mlf_2d_h_val), intent(out) :: D
    End Subroutine mlf_2dplane_h_getValDer
  End Interface
Contains
  Real(c_double) Function mlf_2dplane_h_integrateOnPlane(this, eps, hmax, x0, y0, xMin, yMin, xMax, yMax) Result(res)
    class(mlf_2dplane_h_fun), intent(inout) :: this
    real(c_double), intent(in) :: eps
    integer, intent(in) :: hmax
    real(c_double), intent(in), optional :: x0, y0, xMin, yMin, xMax, yMax
    real(c_double) :: xA0, yA0, F0, dx0, dy0, d2x0, d2y0, R1, R2, R3, R4
    type(mlf_2d_h_val) :: E(3,3)
    xA0 = 0; yA0 = 0
    If(PRESENT(x0)) xA0 = x0
    If(PRESENT(y0)) yA0 = y0
    CALL this%getValDer(xA0, yA0, E(2,2))
    F0 = E(2,2)%val; dx0 = E(2,2)%der(1); dy0 = E(2,2)%der(2)
    d2x0 = E(2,2)%hes(1); d2y0 = E(2,2)%hes(2)
    If(PRESENT(yMin)) Then
      CALL this%getValDer(xA0, yMin, E(2,1))
    Else
      E(2,1)%val = yA0 - SqrtDev(F0, -dy0, d2y0)
    Endif
    If(PRESENT(yMax)) Then
      CALL this%getValDer(xA0, yMax, E(2,3))
    Else
      E(2,3)%val = yA0 + SqrtDev(F0, dy0, d2y0)
    Endif
    If(PRESENT(xMin)) Then
      CALL this%getValDer(xMin, yA0, E(1,2))
      If(PRESENT(yMin)) CALL this%getValDer(xMin, yMin, E(1,1))
      If(PRESENT(yMax)) CALL this%getValDer(xMin, yMax, E(1,3))
    Else
      E(1,2)%val = xA0 - SqrtDev(F0, -dx0, d2x0)
    Endif
    If(PRESENT(xMax)) Then
      CALL this%getValDer(xMax, yA0, E(3,2))
      If(PRESENT(yMin)) CALL this%getValDer(xMax, yMin, E(3,1))
      If(PRESENT(yMax)) CALL this%getValDer(xMax, yMax, E(3,3))
    Else
      E(3,2)%val = xA0 + SqrtDev(F0, dx0, d2x0)
    Endif
    !!$OMP PARALLEL SECTIONS
    !!$OMP SECTION
      R1 = Integr(0.25d0*eps, hmax-2, E(1,1), E(1,2), E(2,1), E(2,2), xMin, yMin, xA0, yA0)
    !!$OMP SECTION
      R2 = Integr(0.25d0*eps, hmax-2, E(1,2), E(1,3), E(2,2), E(2,3), xMin, y0, xA0, yMax)
    !!$OMP SECTION
      R3 = Integr(0.25d0*eps, hmax-2, E(2,1), E(2,2), E(3,1), E(3,2), xA0, yMin, xMax, yA0)
    !!$OMP SECTION
      R4 = Integr(0.25d0*eps, hmax-2, E(2,2), E(2,3), E(3,2), E(3,3), xA0, yA0, xMax, yMax)
    !!$OMP END PARALLEL SECTIONS
    res = R1 + R2 + R3 + R4
  Contains
    Real(c_double) Function SqrtDev(F, D1, D2) Result(y)
      real(c_double), intent(in) :: F, D1, D2
      real(c_double) :: u, v, z
      u = F/D2; v = D1/D2
      z = v**2-u
      If(z<0) Then
        y = ABS(v)
      Else
        z = SQRT(z)
        y = -z-v
        If(y < 0) y = ABS(z-v)
      Endif
    End Function SqrtDev

    Real(c_double) Function SquareI(err, D11, D13, D31, D33, D22, dx, dy) Result(S)
      real(c_double), intent(in) :: dx, dy
      real(c_double), intent(out) :: err
      type(mlf_2d_h_val), intent(in) :: D11, D13, D31, D33, D22
      real(c_double) :: u(0:1, 0:1), ux(0:1, 0:1), uy(0:1, 0:1)
      real(c_double) :: uxx(0:1, 0:1), uxy(0:1, 0:1), uyy(0:1, 0:1)
      real(c_double) :: vxx(0:1, 0:1), vxy(0:1, 0:1), vyy(0:1, 0:1)
      real(c_double) :: Ig, Ih, K
      real(c_double) :: w, wxx, wyy, sxy, IGr
      u   =       RESHAPE([D11%val, D31%val, D13%val, D33%val], [2,2])
      ux  =    dx*RESHAPE([D11%der(1), D31%der(1), D13%der(1), D33%der(1)], [2,2])
      uy  =    dy*RESHAPE([D11%der(2), D31%der(2), D13%der(2), D33%der(2)], [2,2])
      uxx = dx**2*RESHAPE([D11%hes(1), D31%hes(1), D13%hes(1), D33%hes(1)], [2,2])
      uxy = dx*dy*RESHAPE([D11%hes(3), D31%hes(3), D13%hes(3), D33%hes(3)], [2,2])
      uyy = dy**2*RESHAPE([D11%hes(2), D31%hes(2), D13%hes(2), D33%hes(2)], [2,2])
      K = u(1,1)-u(1,0)-u(0,1)+u(0,0)

      Ig = dx*dy*(1d0/4d0 *(u(0,0)+u(0,1)+u(1,0)+u(1,1)) &
                 + 1d0/24d0*(ux(0,0)-ux(1,0)+ux(0,1)-ux(1,1)) &
                 + 1d0/24d0*(uy(0,0)-uy(0,1)+uy(1,0)-uy(1,1)))

      vxx = uxx + 2*RESHAPE([ &
            -ux(1,0) - 2*ux(0,0) + 3*(u(1,0) - u(0,0)), & ! gxx(0,0)
            2*ux(1,0) + ux(0,0) - 3*(u(1,0) - u(0,0)), & ! gxx(1,0)
            -ux(1,1) - 2*ux(0,1) + 3*(u(1,1) - u(1,0)), & ! gxx(0,1)
            2*ux(1,1) + ux(0,1) - 3*(u(1,1) - u(0,1)) & ! gxx(1,1)
          ], [2,2])

      vyy = uyy + 2*RESHAPE([ &
            -uy(0,1) - 2*uy(0,0) + 3*(u(0,1) - u(0,0)), & ! gxx(0,0)
            -uy(1,1) - 2*uy(1,0) + 3*(u(1,1) - u(1,0)), & ! gxx(1,0)
            2*uy(0,1) + uy(0,0) - 3*(u(0,1) - u(0,0)), & ! gxx(0,1)
            2*uy(1,1) + uy(1,0) - 3*(u(1,1) - u(1,0)) & ! gxx(1,1)
          ], [2,2])

      vxy = uxy + K - RESHAPE([ &
            uy(0,1) - uy(0,0) + ux(1,0) - ux(0,0), & ! gxy(0,0)
            uy(0,1) - uy(0,0) + ux(1,1) - ux(1,0), & ! gxy(1,0)
            uy(1,1) - uy(0,1) + ux(1,0) - ux(0,0), & ! gxy(0,1)
            uy(1,1) - uy(1,0) + ux(1,1) - ux(1,0) & ! gxy(1,1)
          ], [2,2])

      Ih = dx*dy*(3d0/720d0*SUM(vxx+vyy) &
                 + 5d0/720d0*(vxy(1,1)-vxy(1,0)-vxy(0,1)+vxy(0,0))) 
      sxy = vxy(1,1) - vxy(1,0) - vxy(0,1) + vxy(0,0)
      w = D22%val - 1d0/4d0*SUM(u) - 1d0/128*SUM(vxx+vyy) - 1d0/64d0*sxy &
        + 1d0/16d0*(ux(1,1)+uy(1,1)-ux(1,0)-uy(1,0)-ux(0,1)-uy(0,1)+ux(0,0)+uy(0,0))
        
      wxx = dx**2*D22%hes(1) - 0.5d0*(ux(1,1)+ux(1,0)-ux(0,1)-ux(0,0)) &
          + 1d0/8d0*(SUM(vxx) + sxy)
      wyy = dy**2*D22%hes(2) - 0.5d0*(uy(1,1)-uy(1,0)+uy(0,1)-uy(0,0)) &
          + 1d0/8d0*(SUM(vyy) + sxy)
      IGr = dx*dy/1575d0*(8*wyy+8*wxx+704*w)
      S = Ih + Ig + IGr
      err = ABS(IGr)
    End Function SquareI

    Real(c_double) Recursive Function I4(eps, hmax, D11, D13, D31, D33, D22, x1, y1, x2, y2, eval) Result(S)
      real(c_double), intent(in) :: eps, x1, y1, x2, y2
      integer, intent(in) :: hmax
      type(mlf_2d_h_val), intent(in) :: D11, D13, D31, D33, D22
      real(c_double), optional, intent(in) :: eval
      type(mlf_2d_h_val) :: D12, D21, D32, D23, U11, U12, U21, U22
      real(c_double) :: dx, dy, x0, y0, xu1, xu2, yu1, yu2, eps0
      real(c_double) :: V(0:1, 0:1), error(0:1, 0:1)
      x0 = 0.5d0*(x1+x2); y0 = 0.5d0*(y1+y2)
      dx = x2-x1; dy = y2-y1
      xu1 = 0.5d0*(x1+x0); xu2 = 0.5d0*(x0+x2)
      yu1 = 0.5d0*(y1+y0); yu2 = 0.5d0*(y0+y2)
      ! Divide the block into four parts
      CALL this%getValDer(x1, y0, D12)
      CALL this%getValDer(x2, y0, D32)
      CALL this%getValDer(x0, y1, D21)
      CALL this%getValDer(x0, y2, D23)

      CALL this%getValDer(xu1, yu1, U11)
      CALL this%getValDer(xu2, yu1, U21)
      CALL this%getValDer(xu1, yu2, U12)
      CALL this%getValDer(xu2, yu2, U22)

      V(0,0) = SquareI(error(0,0), D11, D12, D21, D22, U11, 0.5d0*dx, 0.5d0*dy)
      V(1,0) = SquareI(error(1,0), D21, D22, D31, D32, U21, 0.5d0*dx, 0.5d0*dy)
      V(0,1) = SquareI(error(0,1), D12, D13, D22, D23, U12, 0.5d0*dx, 0.5d0*dy)
      V(1,1) = SquareI(error(1,1), D22, D23, D32, D33, U22, 0.5d0*dx, 0.5d0*dy)

      S = SUM(V)
      If(hmax <= 1) RETURN
      If(PRESENT(eval)) Then
        If(ABS(eval-S) < 2*eps) RETURN
      Endif

      If(SUM(error) <= 2.1d0*eps) RETURN

      eps0 = (2d0*eps - SUM(error, error <= eps*0.5d0))/COUNT(error > eps*0.5d0)

      If(error(0,0) <= eps*0.5d0) Then
        S =     V(0,0)
      Else
        S =     I4(eps0, hmax-2, D11, D12, D21, D22, U11, x1, y1, x0, y0, V(0,0))
      Endif
      If(error(1,0) <= eps*0.5d0) Then
        S = S + V(1,0)
      Else
        S = S + I4(eps0, hmax-2, D21, D22, D31, D32, U21, x0, y1, x2, y0, V(1,0))
      Endif
      If(error(0,1) <= eps*0.5d0) Then
        S = S + V(0,1)
      Else
        S = S + I4(eps0, hmax-2, D12, D13, D22, D23, U12, x1, y0, x0, y2, V(0,1))
      Endif
      If(error(1,1) <= eps*0.5d0) Then
        S = S + V(1,1)
      Else
        S = S + I4(eps0, hmax-2, D22, D23, D32, D33, U22, x0, y0, x2, y2, V(1,1))
      Endif
    End Function I4
    !Real(c_double) Recursive Function I3() Result(S)
    !End Function I3
    Real(c_double) Recursive Function Integr(eps, hmax, D11, D13, D31, D33, x1, y1, x2, y2) Result(S)
      real(c_double), intent(in) :: eps
      integer, intent(in) :: hmax
      type(mlf_2d_h_val), intent(in) :: D11, D13, D31, D33
      real(c_double), intent(in), optional :: x1, y1, x2, y2
      type(mlf_2d_h_val) :: D22
      logical :: pr(4)
      pr = [PRESENT(x1), PRESENT(x2), PRESENT(y1), PRESENT(y2)]
      S = 0
      Select Case(COUNT(pr))
      Case(4)
        CALL this%getValDer(0.5d0*(x1+x2), 0.5d0*(y1+y2), D22)
        S = I4(eps, hmax, D11, D13, D31, D33, D22, x1, y1, x2, y2)
      Case(3)
        If(     .NOT. PRESENT(x1)) Then
        Else If(.NOT. PRESENT(x2)) Then
        Else If(.NOT. PRESENT(y1)) Then
        Else
        Endif
      Case(2)
      End Select
    End Function Integr
  End Function mlf_2dplane_h_integrateOnPlane
End Module mlf_2dplane_integration

