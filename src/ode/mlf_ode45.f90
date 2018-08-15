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

  real(c_double), Parameter :: DOPRI5_A2 = 0.2d0
  real(c_double), Parameter :: DOPRI5_A3(2) = [3d0/40d0, 9d0/40d0]
  real(c_double), Parameter :: DOPRI5_A4(3) = [44d0/45d0, -56d0/15d0, 32d0/9d0]
  real(c_double), Parameter :: DOPRI5_A5(4) = [19372d0/6561d0, -25360d0/2187d0, 64448d0/6561d0, -212d0/729d0]
  real(c_double), Parameter :: DOPRI5_A6(5) = [9017d0/3168d0, -355d0/33d0, 46732d0/5247d0, &
    49d0/176d0, -5103d0/18656d0]
  real(c_double), Parameter :: DOPRI5_A7(6) = [35d0/384d0, 0d0, 500d0/1113d0, 125d0/192d0, &
    -2187d0/6784d0, 11d0/84d0]
  real(c_double), Parameter :: DOPRI5_C(7) = [0d0, 0.2d0, 0.3d0, 0.8d0, 8d0/9d0, 1d0, 1d0]
  real(c_double), Parameter :: DOPRI5_EC(7) = [71d0/57600d0, 0d0, -71d0/16695d0, &
    71d0/1920d0, -17253d0/339200d0, 22d0/525d0, -0.025d0]
  real(c_double), Parameter :: DOPRI5_DC(7) = [-12715105075d0/11282082432d0, 0d0, 87487479700d0/32700410799d0, &
    -10690763975d0/1880347072d0, 701980252875d0/199316789632d0, -1453857185d0/822651844d0, 69997945d0/29380423d0]

  Type, Public, extends(mlf_step_obj) :: mlf_ode45_obj
    real(c_double), pointer :: atoli, rtoli, facMin, facMax, hMax
    real(c_double), pointer :: lastErr, lastTheta, fac, alpha
    real(c_double), pointer :: K(:,:), X0(:), X(:), t0, t, tMax
    real(c_double), allocatable :: Cont(:,:)
    real(c_double) :: lastT
    integer(c_int64_t), pointer :: nAccept, nReject, nStiff, nFun
    integer :: nonStiff, iStiff
    class(mlf_ode_fun), pointer :: fun
  Contains
    procedure :: reinitT => mlf_ode45_reinitT
    procedure :: cancelStep => mlf_ode45_cancelStep
    procedure :: reinit => mlf_ode45_reinit
    procedure :: denseEvaluation => mlf_ode45_denseEvaluation
    procedure :: updateDense => mlf_ode45_updateDense
    procedure :: errorFun => mlf_ode45_errorFun
    procedure :: deltaFun => mlf_ode45_deltaFun
    procedure :: stepF => mlf_ode45_stepFun
    procedure :: findRoot => mlf_ode45_findRoot
    procedure :: init => mlf_ode45_init
    procedure :: stiffDetect => mlf_ode45_stiffDetect
  End Type mlf_ode45_obj

  Integer, Parameter, Public :: mlf_ODE_FunError = -1, mlf_ODE_Stiff = -2 

Contains
  integer Function mlf_ode45_reinit(this) result(info)
    class(mlf_ode45_obj), intent(inout), target :: this
    info = mlf_step_obj_reinit(this)
    this%nAccept = 0; this%nReject = 0; this%nFun = 0
    this%Cont = 0; this%lastT = -HUGE(this%lastT)
    this%lastErr = -1d0; this%lastTheta = -1d0
    this%facMin = 0.2d0; this%facMax = 10d0
    this%hMax = HUGE(this%hMax); this%fac = 0.9
    this%X0 = 0; this%X = 0; this%t0 = 0; this%t = 0
    this%tMax = HUGE(this%tMax); this%alpha = 1.5d0
    this%nStiff = 1000_8; this%nonStiff = 0; this%iStiff = 0
  End Function mlf_ode45_reinit

  integer Function mlf_ode45_reinitT(this, X0, t0, tMax, atoli, rtoli, fac, facMin, &
      facMax, hMax, nStiff) result(info)
    class(mlf_ode45_obj), intent(inout), target :: this
    real(c_double), intent(in), optional :: X0(:), t0, atoli, rtoli, fac, facMin, facMax, hMax, tMax
    integer(c_int64_t), intent(in), optional :: nStiff
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
    if(present(nStiff)) this%nStiff = nStiff
  End Function mlf_ode45_reinitT

  integer Function mlf_ode45_init(this, fun, X0, t0, tMax, atoli, rtoli, fac, facMin, &
      facMax, hMax, nStiff, data_handler) result(info)
    class(mlf_ode45_obj), intent(inout), target :: this
    class(mlf_ode_fun), intent(inout), target :: fun
    class(mlf_data_handler), intent(inout), optional :: data_handler
    real(c_double), intent(in), optional :: X0(:), t0, tMax, atoli, rtoli
    real(c_double), intent(in), optional :: fac, facMin, facMax, hMax
    integer(c_int64_t), intent(in), optional :: nStiff
    type(mlf_step_numFields) :: numFields
    integer(c_int64_t) :: N, nK(2)
    this%fun => fun; N = -1
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
    if(.NOT. ALLOCATED(this%Cont)) ALLOCATE(this%Cont(N,4))
    if(.NOT. present(data_handler)) then
      info = this%reinitT(X0, t0, tMax, atoli, rtoli, fac, facMin, facMax, hMax, nStiff)
    endif
  End Function mlf_ode45_init

  Subroutine mlf_ode45_cancelStep(this)
    class(mlf_ode45_obj), intent(inout) :: this
    this%X = this%X0
    this%t = this%t0
  End Subroutine mlf_ode45_cancelStep

  ! Stiff detection algorithm used by DOPRI5 
  Logical Function mlf_ode45_stiffDetect(this, h, X, Xsti, K7, K6) result(is_stiff)
    class(mlf_ode45_obj), intent(inout) :: this
    real(c_double), intent(in) :: X(:), Xsti(:), K7(:), K6(:), h
    real(c_double) :: Xdist
    is_stiff = .FALSE.
    Xdist = norm2(X-Xsti)
    if(h*norm2(K7-K6) > Xdist*3.25d0) then
      this%nonStiff = 0
      this%iStiff = this%iStiff + 1
      if(this%iStiff == 15) is_stiff = .TRUE.
    else if(this%nonStiff < 6) then
      this%nonStiff = this%nonStiff + 1
      if(this%nonStiff == 6) this%iStiff = 0
    endif
  End Function mlf_ode45_stiffDetect

  real(c_double) Function mlf_ode45_errorFun(this, E, X0, X, theta) result(err)
    class(mlf_ode45_obj), intent(inout) :: this
    real(c_double), intent(in) :: E(:), X0(:), X(:), theta
    real(c_double) :: U(size(E))
    ASSOCIATE( at => this%atoli, rt => this%rtoli)
      U = theta*E/(at+rt*MAX(ABS(X0), ABS(X)))
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

  Integer Function mlf_ode45_findRoot(this, fun, hMax) result(info)
    class(mlf_ode45_obj), intent(inout) :: this
    class(mlf_ode_funCstr), intent(inout) :: fun
    real(c_double), intent(inout) :: hMax
    real(c_double) :: C(size(fun%cstrTmp)), U, dt, h
    integer :: ids(size(fun%cstrTmp)), i, j, id
    info = 0
    CALL fun%getCstr(this%X, C)
    j = 0
    Do i =1,size(C)
      if(C(i) /= 0) then
        U = C(i)*fun%cstrTmp(i)
        if(U > 0) CYCLE
        if(U == 0) then
          U = dot_product(this%K(:,1), fun%cstrVect(:,i))*C(i)
          if(U >= 0) CYCLE
        endif
      endif
      j = j+1
      ids(j) = i
    End Do
    if(j == 0) then
      h = fun%updateCstr(this%t, this%X, C, this%K(:,7))
      hMax = MIN(h, hMax)
      RETURN
    endif
    dt = this%t-this%t0
    hMax = ODE45FindRoot(this%rtoli, this%atoli, ids(1:j), fun%cstrVect, this%K, &
      fun%cstrTmp, C, dt, id)
    call this%denseEvaluation(this%t0+hMax, this%X)
    this%t = this%t0 + hMax
    info = fun%reachCstr(this%t, id, this%X, this%K(:,1))
    h = this%t-this%t0
    hMax = MIN(h, hMax)
  End Function mlf_ode45_findRoot

  real(c_double) Function ODE45FindRoot(rtol, atol, ids, A, K, C0, C, dt, id) result(hMax)
    real(c_double), intent(in) :: A(:,:), K(:, :), C0(:), C(:), rtol, atol, dt
    integer, intent(in) :: ids(:)
    integer, intent(out) :: id
    hMax = dt*Dense45FindRoot(rtol, atol, dt*MATMUL(transpose(A(:,ids)), K), C0(ids), C(ids), id)
    id = ids(id)
  End Function ODE45FindRoot

  ! Find root of the constraints using dense output
  ! CORNER CASE: the case where C0=C(t0)=0 and dC/dt(t0)=Q(:,1)=0 shall be avoided
  real(c_double) Function Dense45FindRoot(rtol, atol, Q, C0, C, id) result(th)
    real(c_double), intent(in) :: Q(:, :), C0(:), C(:), rtol, atol
    integer, intent(out) :: id
    integer :: i
    real(c_double) :: A(size(C),4), C12(size(C)), thI
    real(c_double) :: th1, Y, F, V, W, A5, A6
    A(:,1) = C-C0; A(:,2) = Q(:,1)-A(:,1)
    A(:,3) = -Q(:,7)+A(:,1)-A(:,2)
    A(:,4) = matmul(Q,DOPRI5_DC)
    ! Do a bissection step and then a secant step
    ! Evaluate Y(0.5)
    C12 = C0+0.5d0*(A(:,1)+0.5d0*(A(:,2)+0.5d0*(A(:,3)+0.5d0*A(:,4))))
    th = 2d0
    Do i=1,size(C)
      if(C12(i)*C0(i) >= 0) then
        if(th <= 0.5d0) CYCLE
        if(C(i) == 0) then
          thI = 1d0
        else
          thI = 1d0/(C(i)*(1d0/C(i)-1d0/C12(i)))
        endif
      else
        thI = 1d0/(C12(i)*(1d0/C12(i)-1d0/C0(i)))
      endif
      if(thI < th) then
        th = thI
        id = i
      endif
    End Do
    ! Use Newton-Ralphson to polish the root th
    V = atol+rtol*ABS(C0(id)+C(id))
    ASSOCIATE(A0 => C0(id), A1 => A(id,1), A2 => A(id,2), A3 => A(id,3), A4 => A(id,4))
      A5 = -4d0*A4-2d0*A3
      A6 = A4+A3-A2
      Do i=1,10
        th1 = 1d0-th
        W = A1+th1*(A2+th*(A3+th1*A4))
        F = W+th*(th*(3d0*A4*th+A5)+A6)
        Y = A0+th*W
        th = th-Y/F
        if(abs(Y)<V) EXIT
      End Do
    END ASSOCIATE
  End Function Dense45FindRoot

  Subroutine mlf_ode45_updateDense(this)
    class(mlf_ode45_obj), intent(inout) :: this
    real(c_double) :: dt
    dt = this%t-this%t0
    ASSOCIATE(A => this%Cont, K => this%K)
      A(:,1) = this%X-this%X0
      A(:,2) = dt*K(:,1)-A(:,1)
      A(:,3) = -dt*K(:,7)+A(:,1)-A(:,2)
      A(:,4) = dt*matmul(K,DOPRI5_DC)
    END ASSOCIATE
    this%lastT = this%t
  End Subroutine mlf_ode45_updateDense

  Subroutine mlf_ode45_denseEvaluation(this, t, Y)
    class(mlf_ode45_obj), intent(inout) :: this
    real(c_double), intent(in) :: t
    real(c_double), intent(out) :: Y(:)
    real(c_double) :: th, th1 
    if(this%lastT < this%t) call this%updateDense()
    th = (t-this%t0)/(this%lastT-this%t0)
    th1 = 1d0 - th
    ASSOCIATE(X0 => this%X0, X => this%X, A => this%Cont)
      Y = X0+th*(A(:,1)+th1*(A(:,2)+th*(A(:,3)+th1*A(:,4))))
    END ASSOCIATE
  End Subroutine mlf_ode45_denseEvaluation

  real(c_double) Function mlf_ode45_deltaFun(this, hMax) result(hopt)
    class(mlf_ode45_obj), intent(inout) :: this
    real(c_double), intent(in) :: hMax
    real(c_double) :: d0, d1, d2, dM, h0, h1, hCstr
    integer :: info, i
    hopt = -1
    hCstr = hMax
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
      h0 = MIN(h0, hCstr)
      ! Do an explicit euler step for evaluating the second derivative
      Do i=1,16
        info = this%fun%eval(this%t+h0, this%X+h0*this%K(:,1), this%K(:,2))
        if(info == 0) EXIT ! The point (t+h0,X+h0*F0) is valid
        if(info<0) RETURN
        h0 = 0.5d0*h0 ! The point (t+h0,X+h0*F0) is outside the constraints space
        hCstr = h0
      End Do
      if(info /= 0) RETURN
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
    hopt = min(hopt, hCstr)
    this%lastTheta = hopt
  End Function mlf_ode45_deltaFun

  Integer Function mlf_ode45_stepFun(this, niter) result(info)
    class(mlf_ode45_obj), intent(inout), target :: this
    integer(kind=8), intent(inout), optional :: niter
    integer(kind=8) :: i, niter0
    logical :: wasStopped
    real(c_double) :: h, hMax, err, Xsti(size(this%X))
    real(c_double) :: alphaH, t
    i=1; niter0=1; info = 0
    if(present(niter)) niter0 = niter
    info = mlf_ODE_FunError
    wasStopped = .FALSE.
    if(this%t == this%tMax) then
      info = 1
      niter = 0
      RETURN
    endif
    t = this%t
    ASSOCIATE(K => this%K, X0 =>this%X0, X => this%X, fun => this%fun)
      hMax = this%hMax; alphaH = HUGE(alphaH)
      X0 = X
      info = fun%eval(t, X, K(:,1))
      if(info /=0) RETURN
      this%nFun = this%nFun + 1
      SELECT TYPE(fun)
      class is (mlf_ode_funCstr)
        BLOCK
          real(c_double) :: C(size(fun%cstrTmp))
          CALL fun%getCstr(X0, C)
          hMax = fun%updateCstr(t, X, C, K(:,1))
        END BLOCK
      END SELECT
      do while(i <= niter0)
        this%t0 = t
        h = this%deltaFun(hMax)
        if(info > 0) wasStopped = .TRUE. ! was stopped by a constraint from eval
        if(t+h >= this%tMax) h = this%tMax-t
        if(h<0) RETURN
        info = fun%eval(t+DOPRI5_C(2)*h, X0+h*DOPRI5_A2*K(:,1), K(:,2))
        if(info<0) RETURN; if(info>0) then; hMax = 0.5*DOPRI5_C(2)*h; CYCLE; endif
        info = fun%eval(t+DOPRI5_C(3)*h, X0+h*MATMUL(K(:,1:2), DOPRI5_A3), K(:,3))
        if(info<0) RETURN; if(info>0) then; hMax = DOPRI5_C(2)*h; CYCLE; endif
        info = fun%eval(t+DOPRI5_C(4)*h, X0+h*MATMUL(K(:,1:3), DOPRI5_A4), K(:,4))
        if(info<0) RETURN; if(info>0) then; hMax = DOPRI5_C(3)*h; CYCLE; endif
        info = fun%eval(t+DOPRI5_C(5)*h, X0+h*MATMUL(K(:,1:4), DOPRI5_A5), K(:,5))
        if(info<0) RETURN; if(info>0) then; hMax = DOPRI5_C(4)*h; CYCLE; endif
        ! Ysti is used by DOPRI5 for stiffness detection
        Xsti = X0+h*MATMUL(K(:,1:5), DOPRI5_A6)
        info = fun%eval(t+DOPRI5_C(6)*h, Xsti, K(:,6))
        if(info<0) RETURN; if(info>0) then; hMax = DOPRI5_C(5)*h; CYCLE; endif
        ! Y Contains the value of X(t+h)
        X = X0+h*MATMUL(K(:,1:6), DOPRI5_A7)
        info = fun%eval(t+DOPRI5_C(7)*h, X, K(:,7))
        if(info<0) RETURN; if(info>0) then; hMax = 0.5*SUM(DOPRI5_C(5:6))*h; CYCLE; endif
        this%nFun = this%nFun + 6
        err = this%errorFun(MATMUL(K,DOPRI5_EC), X0, X, h)
        if(err > 1d0) CYCLE
        ! Update X and t
        t = t + h
        this%t = t
        if(MOD(this%nAccept, this%nStiff) == 0 .OR. this%iStiff > 0) then
          if(this%stiffDetect(h, X, Xsti, K(:,7), K(:,6))) then
            info = mlf_ODE_Stiff
            EXIT
          endif
        endif
        ! Check if the constraint is present
        SELECT TYPE(fun)
        class is (mlf_ode_funCstr)
          info = this%findRoot(fun, hMax)
          if(info > 0) EXIT
          if(info < 0) RETURN
        END SELECT
        if(wasStopped) then ! The evaluation constraint is reach
          info = mlf_ODE_HardCstr
          EXIT
        endif
        if(this%tMax <= t) then
            info = mlf_ODE_StopT
            EXIT
        endif
        if(i == niter0) EXIT
        K(:,1) = K(:,7)
        X0 = X
        i = i+1
      end do
    END ASSOCIATE
    if(present(niter)) niter = i
    if(info < 0) info = 0
  End Function mlf_ode45_stepFun
End Module mlf_ode45

