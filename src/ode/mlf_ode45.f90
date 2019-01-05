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
  Use iso_fortran_env
  Use mlf_intf
  Use mlf_rsc_array
  Use mlf_models
  Use mlf_step_algo
  Use mlf_fun_intf
  Use mlf_utils
  Use mlf_ode_class
  IMPLICIT NONE
  PRIVATE

  ! Contains a Runge-Kutta ODE solver based on Dormand and Prince method
  ! This code is based on formulæ from the book Solving Ordinary Differential
  ! Equation I by Hairer, Nørsett and Wanner (DOPRI5)
  ! (2009 Springer, ISBN 978-3-642-05163-0)

  Real(c_double), Parameter :: DOPRI5_A2 = 0.2d0
  Real(c_double), Parameter :: DOPRI5_A3(2) = [3d0/40d0, 9d0/40d0]
  Real(c_double), Parameter :: DOPRI5_A4(3) = [44d0/45d0, -56d0/15d0, 32d0/9d0]
  Real(c_double), Parameter :: DOPRI5_A5(4) = [19372d0/6561d0, -25360d0/2187d0,&
    64448d0/6561d0, -212d0/729d0]
  Real(c_double), Parameter :: DOPRI5_A6(5) = [9017d0/3168d0, -355d0/33d0, &
    46732d0/5247d0, 49d0/176d0, -5103d0/18656d0]
  Real(c_double), Parameter :: DOPRI5_A7(6) = [35d0/384d0, 0d0, 500d0/1113d0, &
    125d0/192d0, -2187d0/6784d0, 11d0/84d0]
  Real(c_double), Parameter :: DOPRI5_C(7) = [0d0, 0.2d0, 0.3d0, 0.8d0, &
    8d0/9d0, 1d0, 1d0]
  Real(c_double), Parameter :: DOPRI5_EC(7) = [71d0/57600d0, 0d0, &
    -71d0/16695d0, 71d0/1920d0, -17253d0/339200d0, 22d0/525d0, -0.025d0]
  Real(c_double), Parameter :: DOPRI5_DC(7) = [-12715105075d0/11282082432d0, &
    0d0, 87487479700d0/32700410799d0, -10690763975d0/1880347072d0, &
    701980252875d0/199316789632d0, -1453857185d0/822651844d0, &
    69997945d0/29380423d0]

  Type, Public, Extends(mlf_ode_algo) :: mlf_ode45_obj
    real(c_double), pointer :: atoli, rtoli, facMin, facMax
    real(c_double), pointer :: lastErr, lastTheta, fac
    integer(c_int64_t), pointer :: nStiff
    integer(c_int64_t) :: nAccept, nReject
    integer :: nonStiff, iStiff
  Contains
    procedure :: reinitT => mlf_ode45_reinitT
    procedure :: reinit => mlf_ode45_reinit
    procedure :: denseEvaluation => mlf_ode45_denseEvaluation
    procedure :: stepF => mlf_ode45_stepFun
    procedure :: init => mlf_ode45_init
    procedure :: initHandler => mlf_ode45_initHandler
  End Type mlf_ode45_obj

  Integer, Parameter, Public :: mlf_ODE_FunError = -1, mlf_ODE_Stiff = -2 

Contains
  Integer Function mlf_ode45_reinit(this) result(info)
    class(mlf_ode45_obj), intent(inout), target :: this
    info = mlf_ode_algo_reinit(this)
    this%nAccept = 0; this%nReject = 0
    this%lastErr = -1d0; this%lastTheta = -1d0
    this%facMin = 0.2d0; this%facMax = 10d0
    this%fac = 0.9
    this%nStiff = 1000_8; this%nonStiff = 0; this%iStiff = 0
  End Function mlf_ode45_reinit

  Integer Function mlf_ode45_reinitT(this, atoli, rtoli, fac, facMin, facMax, &
      nStiff) Result(info)
    class(mlf_ode45_obj), intent(inout), target :: this
    real(c_double), intent(in), optional :: atoli, rtoli, fac, facMin, facMax
    integer(c_int64_t), intent(in), optional :: nStiff
    info = this%reinit()
    If(PRESENT(atoli)) this%atoli = atoli
    If(PRESENT(rtoli)) this%rtoli = rtoli
    If(PRESENT(fac)) this%fac = fac
    If(PRESENT(facMin)) this%facMin = facMin
    If(PRESENT(facMax)) this%facMax = facMax
    If(PRESENT(nStiff)) this%nStiff = nStiff
  End Function mlf_ode45_reinitT

  Integer Function mlf_ode45_initHandler(this, data_handler) result(info)
    class(mlf_ode45_obj), intent(inout), target :: this
    class(mlf_data_handler), intent(inout) :: data_handler
    info = mlf_ode45_init(this, data_handler = data_handler)
  End Function mlf_ode45_initHandler
 
  Integer Function mlf_ode45_init(this, nX, atoli, rtoli, fac, &
      facMin, facMax, nStiff, data_handler) result(info)
    class(mlf_ode45_obj), intent(inout), target :: this
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_double), intent(in), optional :: atoli, rtoli, fac, facMin, facMax
    integer(c_int64_t), intent(in), optional :: nX, nStiff
    type(mlf_step_numFields) :: numFields
    CALL numFields%initFields(nIPar = 1, nRPar = 5, nRVar = 2)
    info = mlf_ode_algo_init(this, numFields, 7, nX, data_handler)
    If(info < 0) RETURN
    ! Integer parameters
    CALL this%addIPar(numFields, this%nStiff, "nStiff")

    ! Real parameters
    CALL this%addRPar(numFields, this%atoli, "atoli")
    CALL this%addRPar(numFields, this%rtoli, "rtoli")
    CALL this%addRPar(numFields, this%fac, "fac")
    CALL this%addRPar(numFields, this%facMin, "facMin")
    CALL this%addRPar(numFields, this%facMax, "facMax")

    ! Real variables
    CALL this%addRVar(numFields, this%lastErr, "lastErr")
    CALL this%addRVar(numFields, this%lastTheta, "lastTheta")
    If(.NOT. PRESENT(data_handler)) Then
      info = this%reinitT(atoli, rtoli, fac, facMin, facMax, nStiff)
    Endif
  End Function mlf_ode45_init

  ! Stiff detection algorithm used by DOPRI5 
  Logical Function ODE45StiffDetect(this, h, X, Xsti, K7, K6) result(is_stiff)
    class(mlf_ode45_obj), intent(inout) :: this
    real(c_double), intent(in) :: X(:), Xsti(:), K7(:), K6(:), h
    real(c_double) :: Xdist
    is_stiff = .FALSE.
    Xdist = norm2(X-Xsti)
    If(h*norm2(K7-K6) > Xdist*3.25d0) Then
      this%nonStiff = 0
      this%iStiff = this%iStiff + 1
      If(this%iStiff == 15) is_stiff = .TRUE.
    Else If(this%nonStiff < 6) Then
      this%nonStiff = this%nonStiff + 1
      if(this%nonStiff == 6) this%iStiff = 0
    Endif
  End Function ODE45StiffDetect

  real(c_double) Function ODE45ErrorFun(this, E, X0, X, theta) result(err)
    class(mlf_ode45_obj), intent(inout) :: this
    real(c_double), intent(in) :: E(:), X0(:), X(:), theta
    real(c_double) :: U(size(E))
    ASSOCIATE( at => this%atoli, rt => this%rtoli)
      U = theta*E/(at+rt*MAX(ABS(X0), ABS(X)))
    END ASSOCIATE
    err = sqrt(sum(U*U)/size(E))
    this%lastTheta = theta
    this%lastErr = err
    If(err <= 1.0) Then
      this%nAccept = this%nAccept + 1
    Else
      this%nReject = this%nReject + 1
    Endif
  End Function ODE45ErrorFun

  Integer Function ODE45FindRoot(this, fun, t, X, hMax) result(info)
    class(mlf_ode45_obj), intent(inout) :: this
    class(mlf_ode_funCstr), intent(inout) :: fun
    real(c_double), intent(inout) :: hMax, t, X(:)
    real(c_double) :: dt, h, hMax0
    integer :: i, N, ids(fun%NCstr), K
    hMax0 = hMax
    info = mlf_ODE_Continue
    N = fun%updateCstr(t, this%X0, X, this%K(:,1), this%K(:,7), ids, hMax)
    If(N == 0) RETURN ! No constraints reached
    If(this%lastT < this%t) CALL ODE45UpdateDense(this)
    dt = this%t - this%t0
    BLOCK
      real(c_double) :: C0(N), C(N), Q(N,7)
      integer :: nIds(N)
      CALL fun%getDerivatives(ids(:N), this%X0, X, this%K, C0, C, Q)
      Q = dt * Q
      K = ODE45FindCstrRoot(this%rtoli, this%atoli, Q, C0, C, nIds, h)
      h = dt * h
      ids(:K) = ids(nIds(:K))
    END BLOCK
    If(K < 0) GOTO 10
    t = this%t0 + h
    CALL this%denseEvaluation(t, X, this%K(:,7))
    info = fun%reachCstr(t, this%t0, this%t, ids(1:K), X, this%K(:,7))
    If(info < 0) GOTO 11
    If(info == mlf_ODE_HardCstr .OR. info == mlf_ODE_StopTime) RETURN
    If(t >= this%tMax) Then
      info = mlf_ODE_StopTime
      RETURN
    Endif
    If(t < this%t0 .OR. t > this%t) GOTO 10
    ! reachCstr makes a supplementary evaluation of the function
    this%nFun = this%nFun + 1
    hMax = MIN(hMax0, fun%getHMax(t, X, this%K(:,7)))
    RETURN
 11 info = -1
    WRITE (error_unit, *) "ODE45 findRoot error with reachCstr: ", info
    WRITE (error_unit, *) "t0: ", this%t0, this%X0
    WRITE (error_unit, *) "t: ", this%t, this%X
    RETURN
 10 info = -1
    WRITE (error_unit, *) "ODE45 findRoot error: h/dt=", (t-this%t0)/dt, " t0=", this%t0, &
      " ids:", ids(:N), " nIds:", ids(:K)
    WRITE (error_unit, *) "t0: ", this%t0, this%X0
    WRITE (error_unit, *) "t: ", this%t, this%X
  End Function ODE45FindRoot

  Real(c_double) Function FindNewtonRalphson(t0, A0, A, V) &
      Result(th)
    real(c_double), intent(in) :: t0, A0, A(6), V
    real(c_double) :: th1, W, F, Y
    integer :: i
    th = t0
    Do i=1,16
      th1 = 1d0-th
      W = A(1)+th1*(A(2)+th*(A(3)+th1*A(4)))
      F = W+th*(th*(3d0*A(4)*th+A(5))+A(6))
      Y = A0+th*W
      th = th-Y/F
      If(abs(Y)<V) RETURN
    End Do
  End Function FindNewtonRalphson

  ! Find root of the constraints using dense output
  ! CORNER CASE: the case where C0=C(t0)=0 and dC/dt(t0)=Q(:,1)=0 shall be avoided
  ! TODO: handle case when C0 and C have the same sign
  Integer Function ODE45FindCstrRoot(rtol, atol, Q, C0, C, ids, th) Result(N)
    real(c_double), intent(in) :: Q(:, :), C0(:), C(:), rtol, atol
    real(c_double), intent(inout) :: th
    integer, intent(out) :: ids(:)
    integer :: i, id
    real(c_double) :: A(SIZE(C),6), Y(SIZE(C)), F(SIZE(C)), thI(SIZE(C)), W(SIZE(C))
    real(c_double) :: V, th1, thMin, thMax
    id = -1
    A(:,1) = C-C0
    A(:,2) = Q(:,1)-A(:,1)
    A(:,3) = -Q(:,7)+A(:,1)-A(:,2)
    A(:,4) = MATMUL(Q,DOPRI5_DC)
    A(:,5) = -4d0*A(:,4)-2d0*A(:,3)
    A(:,6) = A(:,4)+A(:,3)-A(:,2)
    ! Do a bissection step and then a secant step
    ! Evaluate Y(0.5)
    Y = C0+0.5d0*(A(:,1)+0.5d0*(A(:,2)+0.5d0*(A(:,3)+0.5d0*A(:,4))))
    ! Determine if the crossing is negative or positive using C(t0) or dC/dt(t0)
    ! THE CASE C0 == 0 IS NOT CONSIDERED AS A CROSSING!!!
    F = Q(:,1)
    Where(C0 /= 0) F = C0
    Where(Y*F > 0)
      thI = 0.5d0+0.5d0*Y/(Y-C)
    ElseWhere
      thI = 0.5d0*C0/(C0-Y)
    EndWhere
    If(SIZE(C) == 1) Then ! Most common case
      ! Use Newton-Ralphson to polish the root th
      ids(1) = 1
      V = 1d-3*(atol + rtol*ABS(C0(1)-C(1)))
      th = FindNewtonRalphson(thI(1), C0(1), A(1,:), V)
      N = 1
      RETURN
    Endif
    ! Use more intensively the Newton-Ralphson in order to get one unique root
    Do i = 1, 4
      W = A(:,1)+(1d0-thI)*(A(:,2)+thI*(A(:,3)+(1d0-thI)*A(:,4)))
      F = W+thI*(thI*(3d0*A(:,4)*thI+A(:,5))+A(:,6))
      Y = C0+thI*W
      thI = thI-Y/F
    End Do
    thMin = 0; thMax = 1
    id = MINLOC(thI, DIM = 1)
    th = thI(id)
    Do
      th1 = 1d0-th
      W = A(:,1)+th1*(A(:,2)+th*(A(:,3)+th1*A(:,4)))
      F = W+th*(th*(3d0*A(:,4)*th+A(:,5))+A(:,6))
      Y = C0+th*W
      thI = -Y/F
      id = MINLOC(thI, DIM = 1)
      N = COUNT(thI <= 0)
      If(thMax - thMin < 1d-2*(atol+rtol)) Then
        ! When multiple event occurs within a short time span, put multiple event in the list
        ! As the dichotomous search is expensive, this case has to be rare
        N = 0
        Do i = 1, SIZE(C)
          If(thI(i) <= 0) Then
            N = N + 1
            ids(N) = i
          Endif
        End Do
        RETURN
      Endif
      Select Case(N)
      Case(0)
        thMin = th
        th = MIN(th + 2*thI(id), thMax)
      Case(2:)
        thMax = th
        th = 0.5d0*(thMin+thMax)
      Case(1)
        thMax = th
        V = 1d-3*MIN(atol, rtol*ABS(C0(id)-C(id)), ABS(C0(id)))
        th = MAX(MIN(th-Y(id)/F(id), thMax), thMin)
        th = FindNewtonRalphson(th, C0(id), A(id,:), V)
        th = MIN(th, thMax)
        ids(1) = id
        RETURN
      End Select
    End Do
  End Function ODE45FindCstrRoot

  Subroutine ODE45UpdateDense(this)
    class(mlf_ode45_obj), intent(inout) :: this
    real(c_double) :: dt
    If(.NOT. ALLOCATED(this%Cont)) ALLOCATE(this%Cont(SIZE(this%X),4))
    dt = this%t-this%t0
    ASSOCIATE(A => this%Cont, K => this%K)
      A(:,1) = this%X-this%X0
      A(:,2) = dt*K(:,1)-A(:,1)
      A(:,3) = -dt*K(:,7)+A(:,1)-A(:,2)
      A(:,4) = dt*matmul(K,DOPRI5_DC)
    END ASSOCIATE
    this%lastT = this%t
  End Subroutine ODE45UpdateDense

  Subroutine mlf_ode45_denseEvaluation(this, t, Y, K)
    class(mlf_ode45_obj), intent(inout), target :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(out), target :: Y(:)
    real(c_double), intent(out), optional, target :: K(:)
    real(c_double) :: th, th1
    If(this%lastT < this%t) CALL ODE45UpdateDense(this)
    th = (t-this%t0)/(this%lastT-this%t0)
    th1 = (1d0 - th)*th
    ASSOCIATE(X0 => this%X0, A => this%Cont)
      Y = X0+th*A(:,1)+th1*A(:,2)+th*th1*A(:,3)+th1*th1*A(:,4)
      If(PRESENT(K)) Then
        K = A(:,1)+(1-2*th)*A(:,2)+th*(2-3*th)*A(:,3) &
          + 2*th1*(1-2*th)*A(:,4)
      Endif
    END ASSOCIATE
  End Subroutine mlf_ode45_denseEvaluation

  Integer Function ODE45DeltaFun(this, hMax, hOpt) Result(info)
    class(mlf_ode45_obj), intent(inout) :: this
    real(c_double), intent(in) :: hMax
    real(c_double), intent(out) :: hOpt
    real(c_double) :: d0, d1, d2, dM, h0, h1, hCstr
    integer :: i
    info = -1
    hCstr = hMax
    If(this%lastErr <= 0) Then
      ! Starting Step Size as described in II.4 from Hairer and Wanner
      ! Evaluate the norm of |f(x_0,y_0)| and |y_0| using sc_i
      ASSOCIATE(at => this%atoli, rt => this%rtoli, SC => this%K(:,2))
        SC = at+rt*ABS(this%X)
        d0 = NORM2(this%X(:)/SC)
        d1 = NORM2(this%K(:,1)/SC)
      END ASSOCIATE
      ! Do a first guess on hopt
      If(d0 <= 1d-5 .OR. d1 <= 1d-5) Then
        h0 = 1d-6
      Else
        h0 = d0/d1*0.01d0
      Endif
      h0 = MIN(h0, hCstr)
      ! Do an explicit euler step for evaluating the second derivative
      Do i=1,16
        info = this%fun%eval(this%t+h0, this%X+h0*this%K(:,1), this%K(:,2))
        if(info < 0) RETURN
        if(info <= 1) EXIT ! The point (t+h0,X+h0*F0) is valid
        h0 = 0.5d0*h0 ! The point (t+h0,X+h0*F0) is outside the constraints space
        hCstr = h0
      End Do
      If(info /= 0) RETURN
      d2 = NORM2(this%K(:,2)-this%K(:,1))/h0
      dM = MAX(d1,d2)
      If(dM<1d-15) then
        h1 = MAX(1d-6, h0*1e-3)
      Else
        h1 = (0.01d0/dM)**(1d0/5d0)
      Endif
      hopt = MIN(100d0*h0, h1)
    Else
      info = 0
      hopt = this%lastTheta*MIN(this%facMax, MAX(this%facMin, &
        this%fac*(1d0/this%lastErr)**(1d0/5d0)))
    Endif
    hOpt = MIN(hopt, hCstr)
    this%lastTheta = hOpt
  End Function ODE45DeltaFun

  ! Dichotomous search of the point where the hard constraints is triggered
  Real(c_double) Function ODE45SearchHardCstr(this, K, t, hMax) Result(h)
    class(mlf_ode45_obj), intent(inout), target :: this
    real(c_double), intent(inout) :: K(:,:)
    real(c_double), intent(in) :: t, hMax
    real(c_double) :: h0, h1
    integer :: info
    h0 = 0; h1 = hMax
    Do While(h1-h0 > MAX(this%atoli, hMax*this%rtoli))
      h = 0.5d0*(h0+h1)
      info = this%fun%eval(t+h, this%X0+h*DOPRI5_A2*K(:,1), K(:,2))
      If(info < 0) Then
        h = -1; RETURN
      Endif
      If(info <= 1) Then
        h0 = h
      Else
        h1 = h
      Endif
      this%nFun = this%nFun + 1
    End Do
  End Function ODE45SearchHardCstr

  Integer Function mlf_ode45_stepFun(this, niter) result(info)
    class(mlf_ode45_obj), intent(inout), target :: this
    integer(kind=8), intent(inout), optional :: niter
    integer(kind=8) :: i, niter0, nHard
    real(c_double) :: h, hMax, err, Xsti(size(this%X))
    real(c_double) :: t
    logical :: lastHard, stopTime
    i=0; niter0=1; info = 0; nHard = 0
    lastHard = .FALSE.; stopTime = .FALSE.
    If(PRESENT(niter)) niter0 = niter
    info = mlf_ODE_FunError
    If(this%t == this%tMax) Then
      info = mlf_ODE_StopTime
      niter = 0
      RETURN
    Endif
    t = this%t
    this%lastT = -HUGE(this%lastT)
    ASSOCIATE(K => this%K, X0 =>this%X0, X => this%X, fun => this%fun)
      hMax = MIN(this%hMax, this%tMax-this%t)
      X0 = X
      info = fun%eval(t, X, K(:,1))
      If(info /= 0) GOTO 20
      this%nFun = this%nFun + 1
      Select Type(fun)
      Class is (mlf_ode_funCstr)
        hMax = MIN(hMax, fun%getHMax(t, X, K(:,1)))
      End Select
      Do While(i < niter0)
        this%t0 = t
        info = ODE45DeltaFun(this, hMax, h)
        If(t+h >= this%tMax) h = this%tMax-t
        If(h<0) GOTO 12
        info = fun%eval(t+DOPRI5_C(2)*h, X0+h*DOPRI5_A2*K(:,1), K(:,2))
        If(info<0) GOTO 11; If(info>1) Then
          ! Help to find a path that does not trigger a hard constraint
          this%nFun = this%nFun + 1; nHard = nHard + 1
          hMax = ODE45SearchHardCstr(this, K(:,1:2), t, DOPRI5_C(2)*h);
          lastHard = .TRUE.
          If(hMax < 0) GOTO 10
          If(nHard < 5) CYCLE ! Continue the evaluation
          ! Stop the evaluation: cannot avoid hard constaint
          info = mlf_ODE_HardCstr; RETURN
        Endif
        info = fun%eval(t+DOPRI5_C(3)*h, X0+h*MATMUL(K(:,1:2), DOPRI5_A3), K(:,3))
        If(info<0) GOTO 11; If(info>1) Then
          this%nFun = this%nFun + 2; hMax = DOPRI5_C(2)*h; CYCLE
        Endif
        info = fun%eval(t+DOPRI5_C(4)*h, X0+h*MATMUL(K(:,1:3), DOPRI5_A4), K(:,4))
        If(info<0) GOTO 11; If(info>1) Then
          this%nFun = this%nFun + 3; hMax = DOPRI5_C(3)*h; CYCLE
        Endif
        info = fun%eval(t+DOPRI5_C(5)*h, X0+h*MATMUL(K(:,1:4), DOPRI5_A5), K(:,5))
        If(info<0) GOTO 11; If(info>1) Then
          this%nFun = this%nFun + 4; hMax = DOPRI5_C(4)*h; CYCLE
        Endif

        ! Xsti is used by DOPRI5 for stiffness detection
        Xsti = X0+h*MATMUL(K(:,1:5), DOPRI5_A6)
        info = fun%eval(t+DOPRI5_C(6)*h, Xsti, K(:,6))
        If(info<0) GOTO 11; If(info>1) Then
          this%nFun = this%nFun + 5; hMax = DOPRI5_C(5)*h; CYCLE
        Endif

        ! X <- fun(t+h)
        X = X0+h*MATMUL(K(:,1:6), DOPRI5_A7)
        info = fun%eval(t+DOPRI5_C(7)*h, X, K(:,7))
        this%nFun = this%nFun + 6
        If(info<0) GOTO 11; If(info>1) Then;
          hMax = 0.5*SUM(DOPRI5_C(5:6))*h; CYCLE
        Endif
        err = ODE45ErrorFun(this, MATMUL(K,DOPRI5_EC), X0, X, h)
        If(err > 1d0) CYCLE
        ! Update X and t
        t = t + h
        this%t = t
        If(info == 1) stopTime = .TRUE.
        ! If an iteration does not trigger a hard constraint -> reinit nHard counter
        If(.NOT. lastHard) nHard = 0
        If(MOD(this%nAccept, this%nStiff) == 0 .OR. this%iStiff > 0) Then
          If(ODE45StiffDetect(this, h, X, Xsti, K(:,7), K(:,6))) Then
            info = mlf_ODE_Stiff
            EXIT
          Endif
        Endif
        i = i+1
        hMax = MIN(this%hMax, this%tMax-this%t)
        ! Check if the constraint is present
        Select Type(fun)
        Class is (mlf_ode_funCstr)
          info = ODE45FindRoot(this, fun, t, X, hMax)
          If(info < 0) RETURN
          this%t = t
          If(info == mlf_ODE_StopTime .OR. info == mlf_ODE_HardCstr) EXIT
        End Select
        If(this%tMax <= t .OR. stopTime) Then
          info = mlf_ODE_StopTime
          EXIT
        Endif
        If(i == niter0) Then
          info = mlf_ODE_Continue
          EXIT
        Endif
        K(:,1) = K(:,7)
        X0 = X
        lastHard = .FALSE.
      End Do
    END ASSOCIATE
    If(PRESENT(niter)) niter = i
    If(info < 0) info = 0
    RETURN
    ! Handle the case when 
 20 Select Case(info)
    Case(1)
      info = mlf_ODE_StopTime
    Case(2:)
      info = mlf_ODE_HardCstr
    Case(:-1)
      WRITE (error_unit, *) "Error evaluating at start t:", this%t0
      WRITE (error_unit, *) this%X0
    End Select
    RETURN
 12 WRITE (error_unit, *) "h negative ", i
    WRITE (error_unit, *) this%t0, this%X0
    RETURN
 11 WRITE (error_unit, *) "OdeEval error ", i
    WRITE (error_unit, *) this%t0, this%X0
    RETURN
 10 info = -1
    WRITE (error_unit, *) "Error with hardconstraint: ", i
    WRITE (error_unit, *) this%t0, this%X0
  End Function mlf_ode45_stepFun
End Module mlf_ode45

