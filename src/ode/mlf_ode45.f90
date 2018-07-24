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

Module mlf_ode45
  Use ieee_arithmetic
  Use iso_c_binding
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_models
  Use mlf_step_algo
  Use mlf_fun_intf
  Use mlf_utils
  IMPLICIT NONE
  PRIVATE

  ! Contains a Runge-Kutta ODE solver based on Dormand and Prince method
  ! This code is based on formulæ from the book Solving Ordinary Differential Equation I
  ! by Hairer, Nørsett and Wanner (DOPRI5) (2009 Springer, ISBN 978-3-642-05163-0)

  ! TODO: add stiff detection
  real(c_double), Parameter :: A2 = 0.2d0
  real(c_double), Parameter :: A3(2) = [3d0/40d0, 9d0/40d0]
  real(c_double), Parameter :: A4(3) = [44d0/45d0, -56d0/15d0, 32d0/9d0]
  real(c_double), Parameter :: A5(4) = [9372d0/6561d0, -25360d0/2187d0, 64448d0/6561d0, -212d0/729d0]
  real(c_double), Parameter :: A6(5) = [9017d0/3168d0, -355d0/33d0, 46732d0/5247d0, 49d0/176d0, -5103d0/18656d0]
  real(c_double), Parameter :: A7(6) = [35d0/384d0, 0d0, 500d0/1113d0, 125d0/192d0, -2187d0/6784d0, 11d0/84d0]
  real(c_double), Parameter :: C(7) = [0d0, 0.2d0, 0.3d0, 0.8d0, 8d0/9d0, 1d0, 1d0]
  real(c_double), Parameter :: EC(7) = [71d0/57600d0, 0d0, -71d0/16695d0, &
    71d0/1920d0, -17253d0/339200d0, 22d0/525d0, -0.025d0]
  real(c_double), Parameter :: DC(7) = [-12715105075d0/11282082432d0, 0d0, 87487479700d0/32700410799d0, &
    -10690763975d0/1880347072d0, 701980252875d0/199316789632d0, -1453857185d0/822651844d0, 69997945d0/29380423d0]

  Type, Public, Abstract, extends(mlf_step_obj) :: mlf_ode45_obj
    real(c_double), pointer :: atoli, rtoli, facMin, facMax, hMax
    real(c_double), pointer :: lastErr, lastTheta, fac, alpha
    real(c_double), pointer :: K(:,:), X0(:), X(:), t0, t, tMax
    real(c_double), allocatable :: Cont(:,:)
    real(c_double) :: lastT
    integer(c_int64_t), pointer :: nAccept, nReject, nStiff, nFun
    class(mlf_ode_fun), pointer :: fun
  Contains
    procedure :: reinitT => mlf_ode45_reinitT
    procedure :: reinit => mlf_ode45_reinit
    procedure :: denseEvaluation => mlf_ode45_denseEvaluation
    procedure :: updateDense => mlf_ode45_updateDense
    procedure :: errorFun => mlf_ode45_errorFun
    procedure :: deltaFun => mlf_ode45_deltaFun
    procedure :: stepF => mlf_ode45_stepFun
    procedure :: findRoot => mlf_ode45_findRoot
    procedure :: init =>mlf_ode45_init
  End Type mlf_ode45_obj


Contains
  integer Function mlf_ode45_reinit(this) result(info)
    class(mlf_ode45_obj), intent(inout), target :: this
    info = mlf_step_obj_reinit(this)
    this%nAccept = 0; this%nReject = 0; this%nFun = 0
    this%Cont = 0; this%lastT = -HUGE(this%lastT)
    this%lastErr = -1d0; this%lastTheta = -1d0
    this%facMin = 0.2d0; this%facMax = 10d0
    this%hMax = HUGE(this%hMax); this%fac = 0.8
    this%X0 = 0; this%X = 0; this%t0 = 0; this%t = 0
    this%tMax = HUGE(this%tMax); this%alpha = 1.5d0
  End Function mlf_ode45_reinit

  integer Function mlf_ode45_reinitT(this, X0, t0, tMax, atoli, rtoli, fac, facMin, facMax, hMax) result(info)
    class(mlf_ode45_obj), intent(inout), target :: this
    real(c_double), intent(in), optional :: X0(:), t0, atoli, rtoli, fac, facMin, facMax, hMax, tMax
    info = this%reinit()
    if(present(X0)) then
      this%X0 = X0
      this%X = X0
    endif
    if(present(t0)) then
      this%t0 = t0
      this%t = t0
    endif
    if(present(atoli)) this%atoli = atoli
    if(present(rtoli)) this%rtoli = rtoli
    if(present(fac)) this%fac = fac
    if(present(facMin)) this%facMin = facMin
    if(present(facMax)) this%facMax = facMax
    if(present(hMax)) this%hMax = hMax
    if(present(tMax)) this%tMax = tMax
  End Function mlf_ode45_reinitT

  integer Function mlf_ode45_init(this, fun, X0, t0, tMax, atoli, rtoli, fac, facMin, facMax, hMax, data_handler) result(info)
    class(mlf_ode45_obj), intent(inout), target :: this
    class(mlf_ode_fun), intent(inout), target :: fun
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_double), intent(in), optional :: X0(:), t0, tMax, atoli, rtoli
    real(c_double), intent(in), optional :: fac, facMin, facMax, hMax
    type(mlf_step_numFields) :: numFields
    integer(c_int64_t) :: N=-1, nK(2)
    this%fun => fun
    call numFields%initFields(nIPar = 1, nRPar = 7, nRsc = 3, nIVar = 3, nRVar = 5)
    info = mlf_step_obj_init(this, numFields, data_handler)
    if(present(X0)) N = size(X0)
    info = this%add_rarray(numFields, N, this%X0, C_CHAR_"X0", data_handler = data_handler)
    info = this%add_rarray(numFields, N, this%X, C_CHAR_"X", data_handler = data_handler)
    nK = [N, 7_8]
    info = this%add_rmatrix(numFields, nK, this%K, C_CHAR_"K", data_handler = data_handler)
    ! Integer parameters
    call this%addIPar(numFields, this%nStiff, "nStiff")
    ! Integer variables
    call this%addIVar(numFields, this%nAccept, "nAccept")
    call this%addIVar(numFields, this%nReject, "nReject")
    call this%addIVar(numFields, this%nFun, "nFun")
    ! Real parameters
    call this%addRPar(numFields, this%atoli, "atoli")
    call this%addRPar(numFields, this%rtoli, "rtoli")
    call this%addRPar(numFields, this%fac, "fac")
    call this%addRPar(numFields, this%facMin, "facMin")
    call this%addRPar(numFields, this%facMax, "facMax")
    call this%addRPar(numFields, this%hMax, "hMax")
    call this%addRPar(numFields, this%tMax, "tMax")
    ! Real variables
    call this%addRVar(numFields, this%t, "t")
    call this%addRVar(numFields, this%t0, "t0")
    call this%addRVar(numFields, this%lastErr, "lastErr")
    call this%addRVar(numFields, this%lastTheta, "lastTheta")
    call this%addRVar(numFields, this%alpha, "alpha")
    if(.NOT. ALLOCATED(this%Cont)) ALLOCATE(this%Cont(N,3))
    if(.NOT. present(data_handler)) then
      info = this%reinitT(X0, t0, tMax, atoli, rtoli, fac, facMin, facMax, hMax)
    endif
  End Function mlf_ode45_init

  real(c_double) Function mlf_ode45_errorFun(this, E, X0, X, theta) result(err)
    class(mlf_ode45_obj), intent(inout) :: this
    real(c_double), intent(in) :: E(:), X0(:), X(:), theta
    real(c_double) :: U(size(E))
    ASSOCIATE( at => this%atoli, rt => this%rtoli)
      U = U/(at+rt*MAX(ABS(X0), ABS(X)))
    END ASSOCIATE
    err = sqrt(sum(U*U)/size(E))
    this%lastTheta = theta
    this%lastErr = err
    if(err <= 1.0) then
      this%nAccept = this%nAccept + 1
    else
      this%nReject = this%nReject + 1
    endif
  End Function mlf_ode45_errorFun

  Subroutine mlf_ode45_updateDense(this)
    class(mlf_ode45_obj), intent(inout) :: this
    real(c_double) :: dt, dX(size(this%X0))
    dt = this%t-this%t0
    dX = this%X-this%X0
    ASSOCIATE(A => this%Cont, K => this%K)
      A(:,1) = dt*K(:,1)-dX
      A(:,2) = -dt*K(:,7)+dX-A(1,:)
      A(:,3) = dt*matmul(K,DC)
    END ASSOCIATE
    this%lastT = this%t
  End Subroutine mlf_ode45_updateDense

  real(c_double) Function mlf_ode45_findRoot(this, id) result(th)
    class(mlf_ode45_obj), intent(inout) :: this
    integer, intent(in) :: id
    integer :: i
    real(c_double) :: Q0, Q1, Q2, Q3, Q4, Q5, Q6
    real(c_double) :: th1, Y, F, U
    if(this%lastT < this%t) call this%updateDense()
    ASSOCIATE(X0 => this%X0, X => this%X, A => this%Cont)
      Q0 = X0(id); Q1 = X(id)-X0(id)
      Q2 = A(id,1); Q3 = A(id,2); Q4 = A(id,3)
      Q5 = -4d0*Q4-2d0*Q3; Q6 = Q4+Q3-Q2
      ! Do a bissection step and then a secant step
      ! Evaluate Y(0.5)
      Y = Q0+0.5d0*(Q1+0.5d0*(Q2+0.5d0*(Q3+0.5d0*Q4)))
      if(Y*Q0>0d0) then ! The root is between th= 0.5 and 1
        U = 1d0/X(id)
        th = 0.5d0+0.5d0*abs(U*(1d0/Y-U))
      else ! The root is between th= 0 and 0.5
        U = 1d0/Y
        th = 0.5d0*abs(U*(1d0/X0(id)-U))
      endif
      ! Use Newton-Ralphson to polish the root th
      Do i=1,4
        th1 = 1d0-th
        U = Q1+th1*(Q2+th*(Q3+th1*Q4))
        F = U+th*(th*(3d0*Q4*th+Q5)+Q6)
        Y = Q0+th*U
        th = th-Y/F
      End Do
    END ASSOCIATE
  End Function mlf_ode45_findRoot

  Subroutine mlf_ode45_denseEvaluation(this, t, Y)
    class(mlf_ode45_obj), intent(inout) :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(out) :: Y(:)
    real(c_double) :: th, th1 
    if(this%lastT < this%t) call this%updateDense()
    th = (T-this%t0)/(this%T-this%t0)
    th1 = 1d0 - th
    ASSOCIATE(X0 => this%X0, X => this%X, A => this%Cont)
      Y = X0+th*(X-X0+th1*(A(:,1)+th*(A(:,2)+th1*A(:,3))))
    END ASSOCIATE
  End Subroutine mlf_ode45_denseEvaluation

  real(c_double) Function mlf_ode45_deltaFun(this, hMax) result(hopt)
    class(mlf_ode45_obj), intent(inout) :: this
    real(c_double), intent(in) :: hMax
    real(c_double) :: d0, d1, d2, dM, h0, h1
    if(this%lastErr <= 0) then
      ! Starting Step Size as described in II.4 from Hairer and Wanner
      ! Evaluate the norm of |f(x_0,y_0)| and |y_0| using sc_i
      ASSOCIATE(at => this%atoli, rt => this%rtoli, SC => this%K(:,2))
        SC = at+rt*ABS(this%X)
        d0 = norm2(this%X(:)/SC)
        d1 = norm2(this%K(:,1)/SC)
      END ASSOCIATE
      ! Do a first guess on hopt
      if(d0 <= 1d-5 .OR. d1 <= 1d-5) then
        h0 = 1d-6
      else
        h0 = d0/d1*0.01d0
      endif
      h0 = MIN(h0, hMax)
      ! Do an explicit euler step for evaluating the second derivative
      call this%fun%eval(this%t+h0, this%X+h0*this%K(:,1), this%K(:,2))
      d2 = norm2(this%K(:,2)-this%K(:,1))/h0
      dM = max(d1,d2)
      if(dM<1d-15) then
        h1 = max(1d-6, h0*1e-3)
      else
        h1 = (0.01d0/dM)**(1d0/5d0)
      endif
      hopt = min(100d0*h0, h1)
    else
      hopt = this%lastTheta*min(this%facMax, max(this%facMin, &
        this%fac*(1d0/this%lastErr)**(1d0/5d0)))
    endif
    hopt = min(hopt, hMax)
    this%lastTheta = hopt
  End Function mlf_ode45_deltaFun

  Integer Function mlf_ode45_stepFun(this, niter) result(info)
    class(mlf_ode45_obj), intent(inout), target :: this
    integer(kind=8), intent(inout), optional :: niter
    integer(kind=8) :: i
    integer :: idc
    real(c_double) :: h, hMax, err, Y(size(this%X)), Ysti(size(this%X))
    real(c_double) :: th
    logical :: exitLoop = .FALSE.
    info = -1; idC = this%fun%idConst
    ASSOCIATE(t => this%t, K => this%K, X0 =>this%X0, X => this%X)
      hMax = this%hMax
      call this%fun%eval(t, X, K(:,1))
      this%nFun = this%nFun +1
      do i=1,niter
        this%t0 = t
        X0 = X
        if(idC > 0) then
          hMax = min(this%hMax, this%alpha*X0(idC)/K(idC,1))
        endif
        hMax = min(hMax, this%tMax-t)
        h = this%deltaFun(hMax)
        if(h<0) RETURN
        call this%fun%eval(t+C(2)*h, X0+h*A2*K(:,1), K(:,2))
        call this%fun%eval(t+C(3)*h, X0+h*MATMUL(K(:,1:2), A3), K(:,3))
        call this%fun%eval(t+C(4)*h, X0+h*MATMUL(K(:,1:3), A4), K(:,4))
        call this%fun%eval(t+C(5)*h, X0+h*MATMUL(K(:,1:4), A5), K(:,5))
        ! Ysti is used by DOPRI5 for stiffness detection
        Ysti = X0+h*MATMUL(K(:,1:5), A6)
        call this%fun%eval(t+C(6)*h, Ysti, K(:,6))
        ! Y Contains the value of X(t+h)
        Y = X0+h*MATMUL(K(:,1:6), A7)
        call this%fun%eval(t+C(7)*h, Y, K(:,7))
        this%nFun = this%nFun + 6
        err = this%errorFun(MATMUL(K,EC), X0, Y, h)
        if(err > 1d0) CYCLE
        ! Update X and t
        X = Y
        t = t + h
        this%t = t
        if(t > this%tMax) then
          t = this%tMax
          call this%denseEvaluation(t, Y)
          exitLoop = .TRUE.
        endif
        ! Check if the constraint is negative
        if(idC > 0) then
          if(Y(idC) < 0) then ! Check the value at t=min(this%tMax, this%t)
            th = this%findRoot(idC)
            t = this%t0+th*h
            call this%denseEvaluation(t, Y)
            exitLoop = .TRUE.
          endif
        endif
        if(exitLoop) then
            this%t = t
            X = Y
            info = 0
            niter = i
            RETURN
        endif
        K(:,1) = K(:,7)
      end do
    END ASSOCIATE
    info = 0
  End Function mlf_ode45_stepFun
End Module mlf_ode45

