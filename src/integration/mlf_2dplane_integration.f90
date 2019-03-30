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
      class(mlf_2dplane_h_fun), intent(in) :: this
      real(c_double), intent(in) :: x, y
      type(mlf_2d_h_val), intent(out) :: D
    End Subroutine mlf_2dplane_h_getValDer
  End Interface
Contains
  Real(c_double) Function mlf_2dplane_h_integrateOnPlane(this, eps, x0, y0, xMin, yMin, xMax, yMax, nsteps0) Result(res)
    class(mlf_2dplane_h_fun), intent(in) :: this
    real(c_double), intent(in) :: eps
    integer, intent(out), optional :: nsteps0
    real(c_double), intent(in), optional :: x0, y0, xMin, yMin, xMax, yMax
    real(c_double) :: xA0, yA0, F0, dx0, dy0, d2x0, d2y0, eps0, R1, R2, R3, R4
    integer :: nsteps
    type(mlf_2d_h_val) :: E(3,3)
    nsteps = 0
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
    eps0 = MAXVAL(ABS(E%val))*eps
    !$OMP PARALLEL SECTIONS
    !$OMP SECTION
      R1 = Integr(E(1,1), E(1,2), E(2,1), E(2,2), xMin, yMin, xA0, yA0)
    !$OMP SECTION
      R2 = Integr(E(1,2), E(1,3), E(2,2), E(2,3), xMin, y0, xA0, yMax)
    !$OMP SECTION
      R3 = Integr(E(2,1), E(2,2), E(3,1), E(3,2), xA0, yMin, xMax, yA0)
    !$OMP SECTION
      R4 = Integr(E(2,2), E(2,3), E(3,2), E(3,3), xA0, yA0, xMax, yMax)
    !$OMP END PARALLEL SECTIONS
    res = R1 + R2 + R3 + R4
    If(PRESENT(nsteps0)) nsteps0 = nsteps
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
    Real(c_double) Recursive Function I4(D11, D13, D31, D33, x1, y1, x2, y2) Result(S)
      real(c_double), intent(in) :: x1, y1, x2, y2
      type(mlf_2d_h_val), intent(in) :: D11, D13, D31, D33
      real(c_double) :: V(2,2), V1(2,2), vm, dx, dy, dx2, dy2, x0, y0, ax, ay
      dx = 0.5d0*(x2-x1); dy = 0.5d0*(y2-y1)
      dx2 = 0.5d0*dx**2; dy2 = 0.5d0*dy**2
      V(1,1) = D11%val + dx*D11%der(1) + dy*D11%der(2) &
             + dx2*D11%hes(1) + dy2*D11%hes(2) + dx*dy*D11%hes(3)

      V(1,2) = D13%val + dx*D13%der(1) - dy*D13%der(2) &
             + dx2*D13%hes(1) + dy2*D13%hes(2) - dx*dy*D13%hes(3)

      V(2,1) = D31%val - dx*D31%der(1) + dy*D31%der(2) &
             + dx2*D31%hes(1) + dy2*D31%hes(2) - dx*dy*D31%hes(3)

      V(2,2) = D33%val - dx*D33%der(1) - dy*D33%der(2) &
             + dx2*D33%hes(1) + dy2*D33%hes(2) + dx*dy*D33%hes(3)
      vm = 0.25d0*SUM(V)
      If(MAXVAL(ABS(V-vm)) < eps0) Then
        ! Use the second order bivariate Taylor approximation
        S = dx*dy*( D11%val + D13%val + D31%val + D33%val &
                  + 0.5d0*dx*(    D11%der(1) - D13%der(1) + D31%der(1) - D33%der(1)) &
                  + 0.5d0*dy*(    D11%der(2) + D13%der(2) - D31%der(2) - D33%der(2)) &
                  + 1d0/3d0*dx2*( D11%hes(1) + D13%hes(1) + D31%hes(1) + D33%hes(1)) &
                  + 1d0/3d0*dy2*( D11%hes(2) + D13%hes(2) + D31%hes(2) + D33%hes(2)) &
                  + 0.25d0*dx*dy*(D11%hes(3) - D13%hes(3) - D31%hes(3) + D33%hes(3))  )
        nsteps = nsteps+1
      Else
        x0 = 0.5d0*(x1+x2); y0 = 0.5d0*(y1+y2)
        V1(1,1) = (D11%val + dx*D11%der(1) + dx2*D11%hes(1)) &
               - (D31%val - dx*D31%der(1) + dx2*D31%hes(1))

        V1(1,2) = (D11%val + dy*D11%der(2) + dy2*D11%hes(2)) &
               - (D13%val - dy*D13%der(2) + dy2*D13%hes(2))

        V1(2,1) = (D33%val - dy*D33%der(2) + dy2*D33%hes(2)) &
               - (D31%val + dy*D31%der(2) + dy2*D31%hes(2))

        V1(2,2) = (D33%val - dx*D33%der(1) + dx2*D33%hes(1)) &
               - (D13%val + dx*D13%der(1) + dx2*D13%hes(1))
        ax = ABS(V1(1,1))+ABS(V1(2,2))
        ay = ABS(V1(2,1))+ABS(V1(1,2))
        If(ax < 0.2d0*MIN(eps0, ay)) Then
          ! Divide horizontally the block
          BLOCK
            type(mlf_2d_h_val) :: D12, D32
            CALL this%getValDer(x1, y0, D12)
            CALL this%getValDer(x2, y0, D32)
            S =     I4(D11, D12, D31, D32, x1, y1, x2, y0)
            S = S + I4(D12, D13, D32, D33, x1, y0, x2, y2)
          END BLOCK
        Else If(ay < 0.2d0*MIN(eps0, ax)) Then
          ! Divide vertically the block
          BLOCK
            type(mlf_2d_h_val) :: D21, D23
            CALL this%getValDer(x0, y1, D21)
            CALL this%getValDer(x0, y2, D23)
            S =     I4(D11, D13, D21, D23, x1, y1, x0, y2)
            S = S + I4(D21, D23, D31, D33, x0, y1, x2, y2)
          END BLOCK
        Else
          ! Divide the block into four parts
          BLOCK
            type(mlf_2d_h_val) :: D12, D21, D22, D23, D32
            CALL this%getValDer(x0, y0, D22)
            If(ABS(D22%val-vm) < eps0) Then
              S = dx*dy*(8d0/36d0*(D22%val-vm) + D11%val + D13%val + D31%val + D33%val &
                  + 0.5d0*dx*(    D11%der(1) - D13%der(1) + D31%der(1) - D33%der(1)) &
                  + 0.5d0*dy*(    D11%der(2) + D13%der(2) - D31%der(2) - D33%der(2)) &
                  + 1d0/3d0*dx2*( D11%hes(1) + D13%hes(1) + D31%hes(1) + D33%hes(1)) &
                  + 1d0/3d0*dy2*( D11%hes(2) + D13%hes(2) + D31%hes(2) + D33%hes(2)) &
                  + 0.25d0*dx*dy*(D11%hes(3) - D13%hes(3) - D31%hes(3) + D33%hes(3))  )
              nsteps = nsteps+1
            Else
              CALL this%getValDer(x1, y0, D12)
              CALL this%getValDer(x2, y0, D32)
              CALL this%getValDer(x0, y1, D21)
              CALL this%getValDer(x0, y2, D23)
              S =     I4(D11, D12, D21, D22, x1, y1, x0, y0)
              S = S + I4(D12, D13, D22, D23, x1, y0, x0, y2)
              S = S + I4(D21, D22, D31, D32, x0, y1, x2, y0)
              S = S + I4(D22, D23, D32, D33, x0, y0, x2, y2)
            Endif
          END BLOCK
        Endif
      Endif
    End Function I4
    !Real(c_double) Recursive Function I3() Result(S)
    !End Function I3
    Real(c_double) Recursive Function Integr(D11, D13, D31, D33, x1, y1, x2, y2) Result(S)
      real(c_double), intent(in), optional :: x1, y1, x2, y2
      type(mlf_2d_h_val), intent(in) :: D11, D13, D31, D33
      logical :: pr(4)
      pr = [PRESENT(x1), PRESENT(x2), PRESENT(y1), PRESENT(y2)]
      S = 0
      Select Case(COUNT(pr))
      Case(4)
        S = I4(D11, D13, D31, D33, x1, y1, x2, y2)
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

